from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *

from scipy.linalg import cho_factor, cho_solve, lapack
from scipy import spatial


class GaussianProcess_QuasiPeriodicCosineActivity(AbstractModel):
    '''
    This is an altervative version of GaussianProcess_QuasiPeriodicActivity class,
    the only practical difference is that we don't use the george package

    Three parameters out of four are the same for all the datasets, since they are related to
    the properties of the physical process rather than the observed effects on a dataset
     From Grunblatt+2015, Affer+2016
     - theta: is usually related to the rotation period of the star( or one of its harmonics);
     - lambda: is the correlation decay timescale, and it can be related to the lifetime of the active regions.
     - omega: is the length scale of the periodic component, and can be linked to the size evolution of the active regions;
     - h: represents the amplitude of the correlations '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'gp_quasiperiodic_cosine'
        self.internal_likelihood = True

        self.list_pams_common = {
            'Prot',  # Rotational period of the star
            'Pdec',  # Decay timescale of activity
            'Oamp'  # Granulation of activity
        }

        self.list_pams_dataset = {
            'Hamp',  # Amplitude of the signal in the covariance matrix
            'Camp',  # Amplitude of the signal in the covariance matrix
        }

        self._dist_t1 = {}
        self._dist_t2 = {}
        self.inds_cache = {}

    def _compute_distance(self, bjd0, bjd1):
        X0 = np.array([bjd0]).T
        X1 = np.array([bjd1]).T
        return spatial.distance.cdist(X0, X1, lambda u, v: u-v), \
            spatial.distance.cdist(X0, X1, 'sqeuclidean')
        # Change after Rajpaul, Barragàn,  priv. comm.
        #return spatial.distance.cdist(X0, X1, lambda u, v: v-u), \
        #    spatial.distance.cdist(X0, X1, 'sqeuclidean')

    def _compute_cov_matrix(self, variable_value, dist_t1, dist_t2, diagonal_env=None):

        cov_matrix = np.exp(- dist_t2 / (2*variable_value['Pdec']**2)) \
            * (variable_value['Hamp']**2 *
               np.exp(- np.sin(np.pi * dist_t1 / variable_value['Prot']) ** 2.
                      / (2.0 * variable_value['Oamp']**2))
                + variable_value['Camp']**2
                * np.cos(4 * np.pi * dist_t1 / variable_value['Prot']))

        if diagonal_env is not None:
            cov_matrix += np.diag(diagonal_env)

        return cov_matrix

    def initialize_model(self, mc,  **kwargs):

        if kwargs.get('hyperparameters_condition', False):
            self.hyper_condition = self._hypercond_01
        else:
            self.hyper_condition = self._hypercond_00

        if kwargs.get('rotation_decay_condition', False):
            self.rotdec_condition = self._hypercond_02
        else:
            self.rotdec_condition = self._hypercond_00

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        self._dist_t1[dataset.name_ref], self._dist_t2[dataset.name_ref] = self._compute_distance(
            dataset.x0, dataset.x0)
        self.inds_cache[dataset.name_ref] = np.tri(dataset.n, k=-1, dtype=bool)

        return

    def lnlk_compute(self, variable_value, dataset):

        if not self.hyper_condition(variable_value):
            return -np.inf
        if not self.rotdec_condition(variable_value):
            return -np.inf

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

    @staticmethod
    def _hypercond_00(variable_value):
        #Condition from Rajpaul 2017, Rajpaul+2021
        return True

    @staticmethod
    def _hypercond_01(variable_value):
        # Condition from Rajpaul 2017, Rajpaul+2021
        # Taking into account that Pdec^2 = 2*lambda_2^2
        return variable_value['Pdec']**2 > (3. / 4. / np.pi) * variable_value['Oamp']**2 * variable_value['Prot']**2 

    @staticmethod
    def _hypercond_02(variable_value):
        #Condition on Rotation period and decay timescale
        return variable_value['Pdec'] > 2. * variable_value['Prot']
