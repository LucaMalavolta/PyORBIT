from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_gaussian_processes import AbstractGaussianProcesses
from pyorbit.keywords_definitions import *

from scipy.linalg import cho_factor, cho_solve, lapack, LinAlgError
from scipy import matrix, spatial

class GP_Multidimensional_QuasiPeriodicActivity(AbstractModel, AbstractGaussianProcesses):
    ''' Three parameters out of four are the same for all the datasets, since they are related to
    the properties of the physical process rather than the observed effects on a dataset
     From Grunblatt+2015, Affer+2016
     - theta: is usually related to the rotation period of the star( or one of its harmonics);
     - lambda: is the correlation decay timescale, and it can be related to the lifetime of the active regions.
     - omega: is the length scale of the periodic component, and can be linked to the size evolution of the active regions;
     - h: represents the amplitude of the correlations '''

    default_common = 'activity'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

        self.model_class = 'multidimensional_gaussian_process'

        self.internal_likelihood = True
        self.delayed_lnlk_computation = True

        self.list_pams_common = OrderedSet([
            'Prot',  # Rotational period of the star
            'Pdec',  # Decay timescale of activity
            'Oamp',  # Granulation of activity
        ])
        self.list_pams_dataset = OrderedSet([
            'rot_amp', # Amplitude of the first derivative of the covariance matrix
            'con_amp' # Amplitude of the covariance matrix
        ])


        self.internal_parameter_values = None
        self._dist_t1 = None
        self._dist_t2 = None
        self._added_datasets = 0
        self.dataset_ordering = {}
        self.inds_cache = None

        self._dataset_x0 = []
        self._dataset_e2 = []
        self._dataset_names = {}

        self._dataset_nindex = []

        self.use_derivative_dict = {}

        self.internal_coefficients = []

        self._dataset_ej2 = []
        self._dataset_res = []

        self._added_datasets = 0
        self._n_cov_matrix = 0

        self.pi2 = np.pi * np.pi

        ### NEW Addded in PyORBIT v10.10
        self._dataset_rvflag_dict = {}

    def initialize_model(self, mc,  **kwargs):

        self._prepare_hyperparameter_conditions(mc, **kwargs)
        self._prepare_rotation_replacement(mc, **kwargs)
        self._prepare_decay_replacement(mc, **kwargs)

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        """ when reloading the .p files, the object is not reinitialized, so we have to skip the
        incremental addition of datasets if they are already present  """
        if dataset.name_ref in self._dataset_names:
            return

        ## NEW Addded in PyORBIT v10.10
        if dataset.kind == 'RV':
            self._dataset_rvflag_dict[dataset.name_ref] = True
        else:
            self._dataset_rvflag_dict[dataset.name_ref] = False

        self._dataset_nindex.append([self._n_cov_matrix,
                                    self._n_cov_matrix+dataset.n])
        self.internal_coefficients.append([0.00, 0.00])

        self._dataset_x0.append(dataset.x0)

        self._dataset_e2 = np.append(self._dataset_e2, dataset.e**2)

        self._dataset_names[dataset.name_ref] = self._added_datasets
        self._n_cov_matrix += dataset.n
        self._added_datasets += 1

        self._dataset_ej2 = self._dataset_e2 * 1.
        self._dataset_res = self._dataset_e2 * 0.
        self._nugget = self._dataset_e2 * 0. + 1e-7

        self._set_derivative_option(mc, dataset, **kwargs) 

        self.inds_cache = np.tri(self._n_cov_matrix, k=-1, dtype=bool)

        self._dataset_dists = [[0] * self._added_datasets for i in range(self._added_datasets)]
        for l_dataset in range(self._added_datasets):
            for m_dataset in range(self._added_datasets):

                dist_t1, dist_t2 = self._compute_distance(self._dataset_x0[l_dataset],
                                                          self._dataset_x0[m_dataset])
                self._dataset_dists[l_dataset][m_dataset] = [dist_t1, dist_t2]

        return


    def add_internal_dataset(self, parameter_values, dataset):

        self.update_parameter_values(parameter_values)

        self.internal_parameter_values = parameter_values

        d_ind = self._dataset_names[dataset.name_ref]
        d_nstart, d_nend = self._dataset_nindex[d_ind]

        self._dataset_ej2[d_nstart:d_nend] = self._dataset_e2[d_nstart:d_nend] + dataset.jitter**2.0
        self._dataset_res[d_nstart:d_nend] = dataset.residuals

        self.internal_coefficients[d_ind] = [parameter_values['con_amp'], parameter_values['rot_amp']]

    def _compute_distance(self, bjd0, bjd1):
        X0 = np.array([bjd0]).T
        X1 = np.array([bjd1]).T
        return spatial.distance.cdist(X0, X1, lambda u, v: u-v), \
            spatial.distance.cdist(X0, X1, 'sqeuclidean')
        # Change after Rajpaul, Barragàn,  priv. comm.
        #return spatial.distance.cdist(X0, X1, lambda u, v: v-u), \
        #    spatial.distance.cdist(X0, X1, 'sqeuclidean')

    def _compute_submatrix(self, dist_t1, dist_t2, Prot, Prot2, Pdec2, Oamp2):
        phi = 2. * np.pi * dist_t1 / Prot
        sin_phi = np.sin(phi)

        framework_GG = np.exp((-(np.sin(phi / 2.)) ** 2.) / (2.0 * Oamp2)) \
            * np.exp(- dist_t2 / (2 * Pdec2))

        framework_GdG = framework_GG * (- (np.pi * sin_phi / (2 * Prot * Oamp2))
                                        - dist_t1 / Pdec2)

        framework_dGdG = framework_GG * (- (self.pi2 * sin_phi ** 2) / (4. * Prot2 * Oamp2 * Oamp2)
                                         + (self.pi2 * np.cos(phi))
                                         / (Prot2 * Oamp2)
                                         - phi * sin_phi / (2 * Oamp2 * Pdec2)
                                         - dist_t2 / (Pdec2 * Pdec2) + 1. / Pdec2)

        return framework_GG, framework_GdG, framework_dGdG

    def _compute_cov_matrix(self):

        """ Notice the difference in the factor 2 of the decay time scale
            between Grunblatt+2015 (used for the "standard" GP model of PyORBIT) and Rajpaul+2015"""

        # this is faster than computing val**4 several times
        Prot = self.internal_parameter_values['Prot']
        Pdec2 = self.internal_parameter_values['Pdec']**2
        Prot2 = self.internal_parameter_values['Prot']**2
        Oamp2 = self.internal_parameter_values['Oamp']**2

        cov_matrix = np.empty([self._n_cov_matrix, self._n_cov_matrix])

        for l_dataset in range(0, self._added_datasets):
            for m_dataset in range(0, self._added_datasets):

                l_nstart, l_nend = self._dataset_nindex[l_dataset]
                m_nstart, m_nend = self._dataset_nindex[m_dataset]

                #dist_t1, dist_t2 = self._compute_distance(self._dataset_x0[l_dataset],
                #                                          self._dataset_x0[m_dataset])

                dist_t1, dist_t2 = self._dataset_dists[l_dataset][m_dataset]

                Al, Bl = self.internal_coefficients[l_dataset]
                Am, Bm = self.internal_coefficients[m_dataset]

                framework_GG, framework_GdG, framework_dGdG = \
                    self._compute_submatrix(dist_t1, dist_t2,
                                       Prot, Prot2, Pdec2, Oamp2)

                k_lm = Al * Am * framework_GG \
                        + Bl * Bm * framework_dGdG \
                        + Al * Bm * framework_GdG \
                        - Bl * Am * framework_GdG

                cov_matrix[l_nstart:l_nend, m_nstart:m_nend] = k_lm
                #cov_matrix[l_nstart:l_nend, m_nstart:m_nend] = matrix(k_lm)

        cov_matrix += np.diag(self._dataset_ej2)

        return cov_matrix

    def _compute_cov_diag(self, t_array):

        """ Notice the difference in the factor 2 of the decay time scale
            between Grunblatt+2015 (used for the "standard" GP model of PyORBIT) and Rajpaul+2015"""

        # this is faster than computing val**4 several times
        Prot = self.internal_parameter_values['Prot']
        Pdec2 = self.internal_parameter_values['Pdec']**2
        Prot2 = self.internal_parameter_values['Prot']**2
        Oamp2 = self.internal_parameter_values['Oamp']**2

        len_diag = len(t_array)
        cov_diag = np.empty(self._added_datasets*len_diag)
        dist_t1, dist_t2 = np.zeros(len_diag), np.zeros(len_diag)


        for l_dataset in range(0, self._added_datasets):
            l_nstart, l_nend = l_dataset*len_diag,  (l_dataset+1)*len_diag

            Al, Bl = self.internal_coefficients[l_dataset]

            framework_GG, framework_GdG, framework_dGdG = \
                self._compute_submatrix(dist_t1, dist_t2,
                                    Prot, Prot2, Pdec2, Oamp2)

            k_lm = Al * Al * framework_GG + Bl * Bl * framework_dGdG

            cov_diag[l_nstart:l_nend] = k_lm

        return cov_diag


    def _compute_cov_Ks(self, t_array):

        """ Notice the difference in the factor 2 of the decay time scale
            between Grunblatt+2015 (used for the "standard" GP model of PyORBIT) and Rajpaul+2015"""

        # this is faster than computing val**4 several times
        Prot = self.internal_parameter_values['Prot']
        Pdec2 = self.internal_parameter_values['Pdec']**2
        Prot2 = self.internal_parameter_values['Prot']**2
        Oamp2 = self.internal_parameter_values['Oamp']**2

        len_t_array = len(t_array)
        cov_matrix = np.empty([self._added_datasets*len_t_array, self._n_cov_matrix])

        for l_dataset in range(0, self._added_datasets):
            for m_dataset in range(0, self._added_datasets):

                l_nstart, l_nend = l_dataset*len_t_array,  (l_dataset+1)*len_t_array
                m_nstart, m_nend = self._dataset_nindex[m_dataset]


                dist_t1, dist_t2 = self._compute_distance(t_array,
                                                          self._dataset_x0[m_dataset])

                Al, Bl = self.internal_coefficients[l_dataset]
                Am, Bm = self.internal_coefficients[m_dataset]

                framework_GG, framework_GdG, framework_dGdG = \
                    self._compute_submatrix(dist_t1, dist_t2,
                                       Prot, Prot2, Pdec2, Oamp2)

                k_lm = Al * Am * framework_GG \
                        + Bl * Bm * framework_dGdG \
                        + Al * Bm * framework_GdG \
                        - Bl * Am * framework_GdG

                cov_matrix[l_nstart:l_nend, m_nstart:m_nend] = k_lm

        #cov_matrix += np.diag(self._nugget)

        return cov_matrix


    # https://stackoverflow.com/questions/40703042/more-efficient-way-to-invert-a-matrix-knowing-it-is-symmetric-and-positive-semi
    def fast_positive_definite_inverse(self, m):

        cholesky, info = lapack.dpotrf(m)

        if info != 0:
            return None, None, True

        detA = 2*np.sum(np.log(np.diagonal(cholesky)))

        inv, info = lapack.dpotri(cholesky)
        if info != 0:
            return None, None, True

        inv[self.inds_cache] = inv.T[self.inds_cache]
        return inv, detA, False

    def lnlk_compute(self):

        pass_conditions = self.check_hyperparameter_values(self.internal_parameter_values)
        if not pass_conditions:
            return pass_conditions
        
        cov_matrix = self._compute_cov_matrix()

        inv_M, det_A, failed = self.fast_positive_definite_inverse(cov_matrix)

        if failed:
            return -np.inf

        chi2 = np.dot(self._dataset_res, np.matmul(inv_M, self._dataset_res))
        log2_npi = self._n_cov_matrix * np.log(2 * np.pi)
        output = -0.5 * (log2_npi + chi2 + det_A)
        return output

        #cov_matrix = self._compute_cov_matrix(add_diagonal_errors=True)
        #chi2 = np.dot(_3res,np.matmul(inv_M,_3res))
        #
        # try:
        #    alpha = cho_solve(cho_factor(cov_matrix), self._3res)
        #    (s, d) = np.linalg.slogdet(cov_matrix)
        #
        #    return -0.5 * (self.n * np.log(2 * np.pi) +
        #                   np.dot(self._3res, alpha) + d)
        # except:
        #    return -np.inf

    ## NEW Addded in PyORBIT v10.10
    def lnlk_rvonly_compute(self):

        rv_loglike = 0.0

        cov_matrix = self._compute_cov_matrix()
        Ks = cov_matrix - np.diag(self._dataset_ej2)

        alpha = cho_solve(cho_factor(cov_matrix), self._dataset_res)
        mu = np.dot(Ks, alpha).flatten()


        for dataset_name, d_ind in self._dataset_names.items():

            l_nstart, l_nend = self._dataset_nindex[d_ind]

            env = 1. / self._dataset_ej2[l_nstart:l_nend]
            residuals = (self._dataset_res[l_nstart:l_nend] - mu[l_nstart:l_nend])
            n = len(residuals)

            rv_loglike += -0.5 * (n * np.log(2 * np.pi) + np.sum(residuals ** 2 * env - np.log(env)))

        return rv_loglike


    def sample_predict(self, dataset, x0_input=None, return_covariance=False, return_variance=False):

        dataset_index = self._dataset_names[dataset.name_ref]

        if x0_input is None:
            faster_computation = False
            t_predict = dataset.x0
            l_nstart, l_nend = self._dataset_nindex[dataset_index]
        else:
            faster_computation = True
            t_predict = x0_input
            l_nstart, l_nend = len(x0_input)*dataset_index, len(x0_input)*(dataset_index+1)

        cov_matrix = self._compute_cov_matrix()

        if faster_computation:
            Ks = self._compute_cov_Ks(t_predict)
        else:
            Ks = cov_matrix - np.diag(self._dataset_ej2)

        alpha = cho_solve(cho_factor(cov_matrix), self._dataset_res)
        mu = np.dot(Ks, alpha).flatten()
        (s, d) = np.linalg.slogdet(cov_matrix)

        B = cho_solve(cho_factor(cov_matrix), Ks.T)
        KsB_dot_diag = np.diag(np.dot(Ks, B))

        B = None
        Ks = None

        if faster_computation:
            Kss = self._compute_cov_diag(t_predict)
        else:
            Kss = np.diag(cov_matrix)
        cov_matrix = None

        std = np.sqrt(np.array(Kss - KsB_dot_diag).flatten())

        Kss = None

        if return_covariance:
            print('Covariance matrix output not implemented - ERROR')
            quit()

        if return_variance:
            return mu[l_nstart:l_nend], std[l_nstart:l_nend]
        else:
            return mu[l_nstart:l_nend]

    def sample_conditional(self, dataset, x0_input=None):
        val, std = self.sample_predict(dataset, x0_input)
        return val
