from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *

from scipy.linalg import cho_factor, cho_solve, lapack, LinAlgError
from scipy import matrix, spatial


class GaussianProcess_QuasiPeriodicActivity_Derivative(AbstractModel):
    ''' Three parameters out of four are the same for all the datasets, since they are related to
    the properties of the physical process rather than the observed effects on a dataset
     From Grunblatt+2015, Affer+2016
     - theta: is usually related to the rotation period of the star( or one of its harmonics);
     - lambda: is the correlation decay timescale, and it can be related to the lifetime of the active regions.
     - omega: is the length scale of the periodic component, and can be linked to the size evolution of the active regions;
     - h: represents the amplitude of the correlations '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        try:
            import george
        except ImportError:
            print("ERROR: george not installed, this will not work")
            quit()

        self.model_class = 'gp_quasiperiodic_derivative'
        self.internal_likelihood = True

        self.list_pams_common = {
            'Prot',  # Rotational period of the star
            'Pdec',  # Decay timescale of activity
            'Oamp'  # Granulation of activity
        }

        self.list_pams_dataset = {
            'Hamp',  # Amplitude of the signal in the covariance matrix
            'Camp'  # Amplitude of the convective term of the covariance matrix
        }

        self._dist_t1 = {}
        self._dist_t2 = {}
        self.inds_cache = {}

    def _compute_distance(self, bjd0, bjd1):
        X0 = np.array([bjd0]).T
        X1 = np.array([bjd1]).T
        # Change after Rajpaul, priv. comm.
        #return spatial.distance.cdist(X0, X1, lambda u, v: u-v), \
        #    spatial.distance.cdist(X0, X1, 'sqeuclidean')
        return spatial.distance.cdist(X0, X1, lambda u, v: v-u), \
            spatial.distance.cdist(X0, X1, 'sqeuclidean')

    def _compute_cov_matrix(self, variable_value, dist_t1, dist_t2, diagonal_env=None):

        Prot2 = variable_value['Prot']**2
        Pdec2 = variable_value['Pdec']**2
        Oamp2 = variable_value['Oamp']**2
        phi = 2. * np.pi * dist_t1 / variable_value['Prot']

        """ The following equations are slitgly different from those in
            Rajpaul+2015 due to the differnece in the factor 2 of the decay time
            scale with respect to Grunblatt+2015 (uses for the "standard" GP
            model of PyORBIT)
        """

        framework_GG = np.exp((-np.sin(np.pi * dist_t1 / variable_value['Prot']) ** 2.)
                              / (2.0 * Oamp2)
                              - dist_t2 / (2*Pdec2))

        phi = 2. * np.pi * dist_t1 / variable_value['Prot']

        framework_dGdG = framework_GG * (- (np.pi ** 2 * np.sin(phi) ** 2) / (4. * Prot2 * Oamp2 ** 2)
                                         + (np.pi ** 2 * np.cos(phi)) /
                                         (Prot2 * Oamp2)
                                         - phi * np.sin(phi) / (Oamp2 * Pdec2)
                                         - 4 * dist_t2 / Pdec2 ** 2 + 2. / Pdec2)

        cov_matrix = variable_value['Hamp'] ** 2 * framework_GG + \
            variable_value['Camp'] ** 2 * framework_dGdG

        if diagonal_env is not None:
            cov_matrix += np.diag(diagonal_env)

        return cov_matrix

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        self._dist_t1[dataset.name_ref], self._dist_t2[dataset.name_ref] = \
            self._compute_distance(dataset.x0, dataset.x0)
        self.inds_cache[dataset.name_ref] = np.tri(dataset.n, k=-1, dtype=bool)

        return

    def lnlk_compute(self, variable_value, dataset):

        env = dataset.e ** 2.0 + dataset.jitter ** 2.0
        cov_matrix = self._compute_cov_matrix(variable_value,
                                              self._dist_t1[dataset.name_ref],
                                              self._dist_t2[dataset.name_ref],
                                              diagonal_env=env)

        # https://stackoverflow.com/questions/40703042/more-efficient-way-to-invert-a-matrix-knowing-it-is-symmetric-and-positive-semi
        cholesky, info = lapack.dpotrf(cov_matrix)
        if info != 0:
            return -np.inf

        det_A = 2*np.sum(np.log(np.diagonal(cholesky)))
        inv_M, info = lapack.dpotri(cholesky)
        if info != 0:
            return -np.inf

        inv_M[self.inds_cache[dataset.name_ref]
              ] = inv_M.T[self.inds_cache[dataset.name_ref]]

        chi2 = np.dot(dataset.residuals, np.matmul(inv_M, dataset.residuals))
        log2_npi = dataset.n * np.log(2 * np.pi)
        output = -0.5 * (log2_npi + chi2 + det_A)
        return output

    def sample_predict(self, variable_value, dataset, x0_input=None, return_covariance=False, return_variance=False):

        env = dataset.e ** 2.0 + dataset.jitter ** 2.0
        cov_matrix = self._compute_cov_matrix(variable_value,
                                              self._dist_t1[dataset.name_ref],
                                              self._dist_t2[dataset.name_ref],
                                              diagonal_env=env)

        if x0_input is not None:
            predict_t1, predict_t2 = self._compute_distance(x0_input, x0_input)
            crossed_t1, crossed_t2 = self._compute_distance(
                x0_input, dataset.x0)
        else:
            predict_t1 = self._dist_t1[dataset.name_ref]
            predict_t2 = self._dist_t2[dataset.name_ref]
            crossed_t1 = self._dist_t1[dataset.name_ref]
            crossed_t2 = self._dist_t2[dataset.name_ref]

        Ks = self._compute_cov_matrix(variable_value, crossed_t1, crossed_t2)

        alpha = cho_solve(cho_factor(cov_matrix), dataset.residuals)
        mu = np.dot(Ks, alpha).flatten()
        (s, d) = np.linalg.slogdet(cov_matrix)

        Kss = self._compute_cov_matrix(variable_value, predict_t1, predict_t2)

        B = cho_solve(cho_factor(cov_matrix), Ks.T)
        std = np.sqrt(np.array(np.diag(Kss - np.dot(Ks, B))).flatten())

        if return_covariance:
            print('Covariance matrix output not implemented - ERROR')
            quit()

        if return_variance:
            return mu, std
        else:
            return mu

    def sample_conditional(self, variable_value, dataset, x0_input=None):

        val, std = self.sample_predict(variable_value, dataset, x0_input)
        return val