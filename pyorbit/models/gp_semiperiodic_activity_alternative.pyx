from pyorbit.classes.common import *
from pyorbit.models.abstract_common import *
from pyorbit.models.abstract_model import *

from scipy.linalg import cho_factor, cho_solve, LinAlgError
from scipy import spatial

class GaussianProcess_QuasiPeriodicActivity_Alternative(AbstractModel):
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

    internal_likelihood = True

    model_class = 'gp_quasiperiodic_alternative'

    list_pams_common = {
        'Prot', # Rotational period of the star
        'Pdec', # Decay timescale of activity
        'Oamp' # Granulation of activity
    }

    list_pams_dataset = {
        'Hamp',  # Amplitude of the signal in the covariance matrix
    }

    recenter_pams_dataset = {}


    def __init__(self, *args, **kwargs):
        super(GaussianProcess_QuasiPeriodicActivity_Alternative, self).__init__(*args, **kwargs)


        self._dist_t1 = {}
        self._dist_t2 = {}

    def _compute_distance(self, bjd0, bjd1):
        X0 = np.array([bjd0]).T
        X1 = np.array([bjd1]).T
        return spatial.distance.cdist(X0, X1, lambda u, v: u-v), \
               spatial.distance.cdist(X0, X1, 'sqeuclidean')

    def _compute_cov_matrix(self, variable_value, dist_t1, dist_t2, diagonal_env=None):

        cov_matrix = variable_value['Hamp'] ** 2 * \
                     np.exp( (-np.sin(np.pi * dist_t1 / variable_value['Prot']) ** 2.) /
                               (2.0 * variable_value['Oamp']**2)
                               - dist_t2 / variable_value['Pdec']**2)

        if diagonal_env is not None:
            cov_matrix += np.diag(diagonal_env)

        return cov_matrix

    def setup_dataset(self, mc, dataset, **kwargs):

        self._dist_t1[dataset.name_ref], self._dist_t2[dataset.name_ref] = \
            self._compute_distance(dataset.x0, dataset.x0)

        return

    def lnlk_compute(self, variable_value, dataset):
        """ 2 steps:
           1) theta parameters must be converted in physical units (e.g. from logarithmic to linear spaces)
           2) physical values must be converted to {\tt george} input parameters
        """

        env = dataset.e ** 2.0 + dataset.jitter ** 2.0
        cov_matrix = self._compute_cov_matrix(variable_value,
                                              self._dist_t1[dataset.name_ref],
                                              self._dist_t2[dataset.name_ref],
                                              diagonal_env=env)
        try:
            alpha = cho_solve(cho_factor(cov_matrix), dataset.residuals)
            (s, d) = np.linalg.slogdet(cov_matrix)
            return -0.5 * (dataset.n * np.log(2 * np.pi) + np.dot(dataset.residuals, alpha) + d)
        except LinAlgError:
            return -np.inf

    def sample_predict(self, variable_value, dataset, x0_input=None):

        env = dataset.e ** 2.0 + dataset.jitter ** 2.0
        cov_matrix = self._compute_cov_matrix(variable_value,
                                              self._dist_t1[dataset.name_ref],
                                              self._dist_t2[dataset.name_ref],
                                              diagonal_env=env)

        if x0_input is not None:
            predict_t1, predict_t2 = self._compute_distance(x0_input, x0_input)
            crossed_t1, crossed_t2 = self._compute_distance(x0_input, dataset.x0)
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
        return mu, std

    def sample_conditional(self, variable_value, dataset, x0_input=None):

        val, std = self.sample_predict(variable_value, dataset, x0_input)
        return val
