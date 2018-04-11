from common import *
from model_container_abstract import ModelContainer


class ModelContainerPolyChord(ModelContainer):
    # Default values, taken from the PyPolyChord wrapper in PolyChord official distribution, V1.9
    include_priors = False
    polychord_parameters = {'nlive_mult': 25,
                                 'num_repeats_mult': 5,
                                 'feedback': 1,
                                 'precision_criterion': 0.001,
                                 'max_ndead': -1,
                                 'boost_posterior': 0.0,
                                 'read_resume': True,
                                 'base_dir': 'polychord/',
                                 'shutdown_jitter': False}
    polychord_dir_output = None

    def polychord_priors(self, cube):
        theta = (self.bounds[:, 1] - self.bounds[:, 0]) * cube + self.bounds[:, 0]
        return theta.tolist()

    def polychord_call(self, theta1):
        #theta = np.empty(self.ndim)
        #for i in xrange(0, self.ndim):
        #    theta[i] = theta1[i]
        theta = [theta1[i] for i in xrange(0, self.ndim)]
        phi = [0.0] * 0
        chi_out = self(theta)
        if chi_out < -0.5e10:
            return -0.5e10, phi
        return chi_out, phi