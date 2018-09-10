from common import *
from model_container_abstract import ModelContainer


class ModelContainerPolyChord(ModelContainer):

    def __init__(self):
        super(self.__class__, self).__init__()

        # Default values, taken from the PyPolyChord wrapper in PolyChord official distribution, V1.9
        self.include_priors = True
        self.nested_sampling_parameters = {'shutdown_jitter': False,
                                     'include_priors': True}

        self.output_directory = None

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