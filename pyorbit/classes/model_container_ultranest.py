from pyorbit.subroutines.common import np, nested_sampling_prior_compute
from pyorbit.classes.model_container_abstract import ModelContainer


class ModelContainerUltranest(ModelContainer):

    def __init__(self):
        super(self.__class__, self).__init__()

        # Default values, taken from the PyPolyChord wrapper in PolyChord official distribution, V1.9
        self.include_priors = False
        self.nested_sampling_parameters = {'shutdown_jitter': False,
                                           'nlive': 500,
                                           'desired_accuracy': 0.5,
                                           'update_interval_iter_fraction': 0.4,
                                           'min_ess':400,
                                           'improvement_loops': 3,
                                           'include_priors': False}

        self.output_directory = None

    def ultranest_transform(self, cube):
        theta = np.zeros(len(cube), dtype=np.double)

        for i in range(0, len(cube)):

            theta[i] = nested_sampling_prior_compute(
                cube[i], self.priors[i][0], self.priors[i][2], self.spaces[i])

        return theta

    def ultranest_call(self, theta):

        chi_out = self(theta, self.include_priors)

        if chi_out < -0.5e10:
            return -0.5e10
        return chi_out
