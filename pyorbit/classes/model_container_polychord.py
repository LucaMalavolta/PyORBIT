from pyorbit.subroutines.common import *
from pyorbit.classes.model_container_abstract import ModelContainer


class ModelContainerPolyChord(ModelContainer):

    def __init__(self):
        super(self.__class__, self).__init__()

        # Default values, taken from the PyPolyChord wrapper in PolyChord official distribution, V1.9
        self.include_priors = False
        self.nested_sampling_parameters = {'shutdown_jitter': False,
                                     'include_priors': False}

        self.output_directory = None

    def polychord_priors(self, cube):
        #theta = []
        #
        #for i in range(0, len(cube)):
        #    theta.append(nested_sampling_prior_compute(cube[i], self.priors[i][0], self.priors[i][2], self.spaces[i]))

        theta = cube*0.
        for i in range(0, len(theta)):
            theta[i] = nested_sampling_prior_compute(cube[i], self.priors[i][0], self.priors[i][2], self.spaces[i])

        return theta

    def polychord_call(self, theta):

        #theta = [theta1[i] for i in range(0, self.ndim)]
        phi = [0.0] * 0 
        chi_out = self(theta, self.include_priors)
        if chi_out < -0.5e30:
            return -0.5e30, phi
        return chi_out, phi

