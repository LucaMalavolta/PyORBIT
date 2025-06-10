from pyorbit.subroutines.common import np, nested_sampling_prior_compute
from pyorbit.classes.model_container_abstract import ModelContainer


class ModelContainerNautilus(ModelContainer):

    def __init__(self):
        super(self.__class__, self).__init__()

        # Default values, taken from the PyPolyChord wrapper in PolyChord official distribution, V1.9
        self.include_priors = False
        self.nested_sampling_parameters = {'shutdown_jitter': False,
                                           'nthreads': 16,
                                           'include_priors': False,
                                           'default': False,
                                           'use_default': False,
                                           'use_threading_pool':True }

        self.output_directory = None

    def nautilus_priors(self, cube):
        theta = np.zeros(len(cube), dtype=np.double)

        for i in range(0, len(cube)):
            theta[i] = nested_sampling_prior_compute(
                cube[i], self.priors[i][0], self.priors[i][2], self.spaces[i])
        return theta

    def nautilus_call(self, theta):

        chi_out = self(theta, self.include_priors)

        """ check for finite log_likelihood"""
        if  np.isnan(chi_out):
            chi_out = -np.inf

        return chi_out
