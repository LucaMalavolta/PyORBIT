from pyorbit.subroutines.common import np, nested_sampling_prior_compute
from pyorbit.classes.model_container_abstract import ModelContainer


class ModelContainerNestle(ModelContainer):

    def __init__(self):
        super(self.__class__, self).__init__()

        # Default values, taken from the PyPolyChord wrapper in PolyChord official distribution, V1.9
        self.include_priors = False
        self.nested_sampling_parameters = {'shutdown_jitter': False,
                                           'nthreads': 16,
                                           'include_priors': False,
                                           'dlogz':0.01,
                                           'pfrac': 1.0,
                                           'sample': 'auto',
                                           'bound': 'multi'}

        self.output_directory = None

    def nestle_priors(self, cube):
        theta = np.zeros(len(cube), dtype=np.double)

        for i in range(0, len(cube)):
            theta[i] = nested_sampling_prior_compute(
                cube[i], self.priors[i][0], self.priors[i][2], self.spaces[i])
        return theta

    def nestle_call(self, theta):

        chi_out = self(theta, self.include_priors)

        if chi_out < -0.5e10:
            return -0.5e10
        return chi_out
