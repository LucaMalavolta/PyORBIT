from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *

from scipy.linalg import cho_factor, cho_solve
from scipy import matrix
from scipy import spatial
from scipy.linalg import lapack

class GP_Framework_QuasiPeriodicActivity(AbstractModel):
    ''' Three parameters out of four are the same for all the datasets, since they are related to
    the properties of the physical process rather than the observed effects on a dataset
     From Grunblatt+2015, Affer+2016
     - theta: is usually related to the rotation period of the star( or one of its harmonics);
     - lambda: is the correlation decay timescale, and it can be related to the lifetime of the active regions.
     - omega: is the length scale of the periodic component, and can be linked to the size evolution of the active regions;
     - h: represents the amplitude of the correlations '''

    internal_likelihood = True
    delayed_lnlk_computation = True

    model_class = 'gp_framework_semiperiod_activity'

    list_pams_common = {
        'Prot', # Rotational period of the star
        'Pdec', # Decay timescale of activity
        'Oamp', # Granulation of activity
        'Vc',
        'Vr',
        'Lc',
        'Bc',
        'Br'
    }

    list_pams_model = {}

    list_pams_dataset = {}

    recenter_pams_dataset = {}

    n_pams = 8

    """ Indexing is determined by the way the kernel is constructed, so it is specific of the Model and not of the 
    Common class"""
    gp_pams_index = {
        'Prot': 0,
        'Pdec': 1,
        'Oamp': 2,
        'Vc': 3,
        'Vr': 4,
        'Lc': 5,
        'Bc': 6,
        'Br': 7
    }



    def __init__(self, *args, **kwargs):
        super(GP_Framework_QuasiPeriodicActivity, self).__init__(*args, **kwargs)

        self.internal_gp_pams = None

        self._x0 = None
        self._nx0 = None
        self._3x0 = None
        self._3e = None
        self._3ej = None
        self._3res = None
        self._dist_t1 = None
        self._dist_t2 = None
        self._added_datasets = 0
        self.dataset_ordering = {}
        self.inds_cache = None


    def convert_val2gp(self, input_pams):
        """
        :param input_pam: dictonary with the 'physically meaningful' parameters of the GP kernel
        :return: dictonary with the parameters to be fed to 'george'
        WARNING: this subroutine is HIGHLY specific of your choice of the kernel! I reccomend to
        create a new Class with different transformations if you are planning of using a different
        kernel combination
        """
        output_pams = np.zeros(self.n_pams, dtype=np.double)

        """ You must check _george_ documentation (and possibily do a lot of testing) to know how to convert physical 
        values to the parameter vector accepted by george.set_parameter_vector() function. Note: these values may be 
        different from ones accepted by the kernel
        """
        output_pams[self.gp_pams_index['Prot']] = input_pams['Prot']
        output_pams[self.gp_pams_index['Pdec']] = input_pams['Pdec']
        output_pams[self.gp_pams_index['Oamp']] = input_pams['Oamp']
        output_pams[self.gp_pams_index['Vc']] = input_pams['Vc']
        output_pams[self.gp_pams_index['Vr']] = input_pams['Vr']
        output_pams[self.gp_pams_index['Lc']] = input_pams['Lc']
        output_pams[self.gp_pams_index['Bc']] = input_pams['Bc']
        output_pams[self.gp_pams_index['Br']] = input_pams['Br']

        return output_pams

    def convert_gp2val(self, input_pams):
        """
        :param input_pam: dictonary with the parameters to be fed to 'george'
        :return: dictonary with the 'physically meaningful' parameters of the GP kernel
        WARNING: this subroutine is HIGHLY specific of your choice of the kernel! I reccomend to
        create a new Class with different transformations if you are planning of using a different
        kernel combination
        """
        return {
            'Pdec': input_pams[self.gp_pams_index['Pdec']],
            'Oamp': input_pams[self.gp_pams_index['Oamp']],
            'Prot': input_pams[self.gp_pams_index['Prot']],
            'Vc': input_pams[self.gp_pams_index['Vc']],
            'Vr': input_pams[self.gp_pams_index['Vr']],
            'Lc': input_pams[self.gp_pams_index['Lc']],
            'Bc': input_pams[self.gp_pams_index['Bc']],
            'Br': input_pams[self.gp_pams_index['Br']],
        }

    def initialize_model(self, mc,  **kwargs):
        pass

    def setup_dataset(self, mc, dataset, **kwargs):

        if self._nx0:
            if not dataset.n == self._nx0:
                raise ValueError("GP framework error: inconsistent datasets")

        else:
            self._nx0 = dataset.n
            self._3x0 = np.zeros(3*dataset.n, dtype=np.double)
            self._3e = np.zeros(3*dataset.n, dtype=np.double)
            self._3ej = np.zeros(3*dataset.n, dtype=np.double)
            self._3res = np.zeros(3*dataset.n, dtype=np.double)

        if dataset.kind == 'RV':
            self._3x0[:self._nx0] = dataset.x0
            self._3e[:self._nx0] = dataset.e
            self._added_datasets += 1

        if dataset.kind == 'H-alpha' or \
                dataset.kind == 'S_index' or \
                dataset.kind == 'Ca_HK':
            self._3x0[self._nx0:2*self._nx0] = dataset.x0
            self._3e[self._nx0:2*self._nx0] = dataset.e
            self._added_datasets += 1

        if dataset.kind == 'BIS':
            self._3x0[2*self._nx0:] = dataset.x0
            self._3e[2*self._nx0:] = dataset.e
            self._added_datasets += 1

        if self._added_datasets == 3:

            self.inds_cache = np.tri(3*self._nx0, k=-1, dtype=bool)

            if self._3x0[:self._nx0].all() == self._3x0[self._nx0:2*self._nx0].all() == self._3x0[2*self._nx0:].all():
                self._x0 = self._3x0[:self._nx0]
                self._dist_t1, self._dist_t2 = self._compute_distance(self._x0, self._x0)
            else:
                raise ValueError("GP framework error: inconsistent datasets")

        return

    ### WHICH ONE SHOULD I KEEP???
    def common_initialization_with_dataset(self, dataset):
        self.define_kernel(dataset)
        return

    def define_kernel(self):

        # Prot, Pdec, Oamp
        return

    """
    def add_internal_dataset(self, variable_value, dataset, reset_status):

        self.internal_gp_pams = self.convert_val2gp(variable_value)

        if dataset.kind == 'RV':
            self._rv_x = dataset.x0
            self._rv_y = dataset.y - dataset.model
            self._rv_e = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
            scatter = np.std(self.gp_framwork.rv)
            self._rv /= scatter
            self._rverr /= scatter

        if dataset.kind == 'BIS':
            self.gp_framwork.bis = dataset.y - dataset.model
            self.gp_framwork.sig_bis = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
            scatter = np.std(self.gp_framwork.bis)
            self.gp_framwork.bis /= scatter
            self.gp_framwork.sig_bis /= scatter

        if dataset.kind == 'H-alpha' or \
                dataset.kind == 'S_index' or \
                dataset.kind == 'Ca_HK':
            self.gp_framwork.rhk = dataset.y - dataset.model
            self.gp_framwork.sig_rhk = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
            scatter = np.std(self.gp_framwork.rhk)
            self.gp_framwork.rhk /= scatter
            self.gp_framwork.sig_rhk /= scatter

    """

    def add_internal_dataset(self, variable_value, dataset, reset_status):
        if not reset_status:
            self.internal_gp_pams = self.convert_val2gp(variable_value)

        if dataset.kind == 'RV':
            self._3ej[:self._nx0] = np.sqrt(self._3e[:self._nx0]**2 + dataset.jitter**2.0)
            self._3res[:self._nx0] = dataset.residuals

        if dataset.kind == 'H-alpha' or \
                dataset.kind == 'S_index' or \
                dataset.kind == 'Ca_HK':
            self._3ej[self._nx0:2 * self._nx0] = \
                np.sqrt(self._3e[self._nx0:2 * self._nx0]**2 + dataset.jitter**2.0)
            self._3res[self._nx0:2 * self._nx0] = dataset.residuals

        if dataset.kind == 'BIS':
            self._3ej[2*self._nx0:] = np.sqrt(self._3e[2*self._nx0:]**2 + dataset.jitter**2.0)
            self._3res[2 * self._nx0:] = dataset.residuals

    def _compute_distance(self, bjd0, bjd1):
        X0 = np.array([bjd0]).T
        X1 = np.array([bjd1]).T
        return spatial.distance.cdist(X0, X1, lambda u, v: u-v), \
            spatial.distance.cdist(X0, X1, 'sqeuclidean')

    def _compute_cov_matrix(self, add_diagonal_errors=False, bjd0=None, bjd1=None):

        if bjd0 is None and bjd1 is None:
            dist_t1 = self._dist_t1
            dist_t2 = self._dist_t2
        else:
            dist_t1, dist_t2 = self._compute_distance(bjd0, bjd1)

        cov_matrix = np.empty([np.size(dist_t1, axis=0) * 3, np.size(dist_t1, axis=1) * 3])

        Prot = self.internal_gp_pams[0]
        Pdec = self.internal_gp_pams[1]
        Oamp = self.internal_gp_pams[2]
        Vc = self.internal_gp_pams[3]
        Vr = self.internal_gp_pams[4]
        Lc = self.internal_gp_pams[5]
        Bc = self.internal_gp_pams[6]
        Br = self.internal_gp_pams[7]

        #this is faster than computing val**4 several times
        Pdec2 = Pdec * Pdec
        Prot2 = Prot * Prot
        Oamp2 = Oamp * Oamp
        pi2 = np.pi * np.pi

        phi = 2. * np.pi * dist_t1 / Prot
        sin_phi = np.sin(phi)

        framework_GG = np.exp((-(np.sin(phi / 2.)) ** 2.) / (2.0 * Oamp2)) \
            * np.exp(- dist_t2 / (2 * Pdec2))

        framework_GdG = framework_GG * (- (np.pi * sin_phi / (2 * Prot * Oamp2)) \
                    - dist_t1 / Pdec2)

        framework_dGdG = framework_GG * ( - (pi2 * sin_phi ** 2) / (4. * Prot2 * Oamp2 * Oamp2) \
                    + (pi2 * np.cos(phi)) / (Prot2 * Oamp2) \
                    - phi * sin_phi / (2 * Oamp2 * Pdec2) \
                    - dist_t2 / (Pdec2 * Pdec2) + 1. / Pdec2)

        prod_11 = Vc ** 2 * framework_GG + Vr ** 2 * framework_dGdG
        prod_22 = Lc ** 2 * framework_GG
        prod_33 = Bc ** 2 * framework_GG + Br ** 2 * framework_dGdG
        prod_12 = Vc * Lc * framework_GG - Vr * Lc * framework_GdG
        prod_13 = Vc * Bc * framework_GG + Vr * Br * framework_dGdG + (Vc * Br - Vr * Bc) * framework_GdG
        prod_23 = Lc * Bc * framework_GG + Lc * Br * framework_GdG

        k_11 = matrix(prod_11)
        k_22 = matrix(prod_22)
        k_33 = matrix(prod_33)
        k_12 = matrix(prod_12)
        k_13 = matrix(prod_13)
        k_23 = matrix(prod_23)
        k_21 = k_12.T
        k_31 = k_13.T
        k_32 = k_23.T

        xs = np.size(dist_t1, axis=0)
        ys = np.size(dist_t1, axis=1)

        cov_matrix[:xs, :ys] = k_11
        cov_matrix[:xs, ys:2*ys] = k_12
        cov_matrix[:xs, 2*ys:] = k_13

        cov_matrix[xs:2*xs, :ys] = k_21
        cov_matrix[xs:2*xs, ys:2*ys] = k_22
        cov_matrix[xs:2*xs, 2*ys:] = k_23

        cov_matrix[2*xs:,:ys] = k_31
        cov_matrix[2*xs:,   ys:2*ys] = k_32
        cov_matrix[2 * xs:, 2*ys:] = k_33

        if add_diagonal_errors:
            cov_matrix += np.diag(self._3ej ** 2)

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

        cov_matrix = self._compute_cov_matrix(add_diagonal_errors=True)
        inv_M, det_A, failed = self.fast_positive_definite_inverse(cov_matrix)
        if failed:
            return -np.inf
        chi2 = np.dot(self._3res,np.matmul(inv_M,self._3res))
        log2_npi = self._nx0 * np.log(2 * np.pi)
        output = -0.5 * (log2_npi + chi2 + det_A)
        return output

        #cov_matrix = self._compute_cov_matrix(add_diagonal_errors=True)
        #chi2 = np.dot(_3res,np.matmul(inv_M,_3res))
        #
        #try:
        #    alpha = cho_solve(cho_factor(cov_matrix), self._3res)
        #    (s, d) = np.linalg.slogdet(cov_matrix)
        #
        #    return -0.5 * (self.n * np.log(2 * np.pi) +
        #                   np.dot(self._3res, alpha) + d)
        #except:
        #    return -np.inf


    def sample_predict(self, dataset, x0_input=None, return_covariance=False, return_variance=False):

        if x0_input is None:
            t_predict = self._3x0
            n_output = self._nx0
        else:
            t_predict = np.concatenate([x0_input, x0_input, x0_input])
            n_output = np.size(x0_input, axis=0)

        cov_matrix = self._compute_cov_matrix(add_diagonal_errors=True)
        Ks = self._compute_cov_matrix(add_diagonal_errors=False,
                                      bjd0=t_predict,
                                      bjd1=self._3x0)

        alpha = cho_solve(cho_factor(cov_matrix), self._3res)
        mu = np.dot(Ks, alpha).flatten()
        (s, d) = np.linalg.slogdet(cov_matrix)

        Kss = self._compute_cov_matrix(add_diagonal_errors=False,
                                 bjd0=t_predict,
                                 bjd1=t_predict)

        B = cho_solve(cho_factor(cov_matrix), Ks.T)
        std = np.sqrt(np.array(np.diag(Kss - np.dot(Ks, B))).flatten())

        if return_covariance:
            print('Covariance matrix output not implemented - ERROR')
            quit()

        if dataset.kind == 'RV':
            if return_variance:
                return mu[:n_output], std[:n_output]
            else:
                return mu[:n_output]

        if dataset.kind == 'H-alpha' or \
                dataset.kind == 'S_index' or \
                dataset.kind == 'Ca_HK':
            if return_variance:
                return mu[n_output:2*n_output], std[n_output:2*n_output]
            else:
                return mu[n_output:2*n_output]

        if dataset.kind == 'BIS':
            if return_variance:
                return mu[2*n_output:], std[2*n_output:]
            else:
                return mu[2*n_output:]

    def sample_conditional(self, dataset, x0_input=None):
        val, std = self.sample_predict(dataset, x0_input)
        return val

