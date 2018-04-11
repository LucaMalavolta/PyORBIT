from common import *
from model_container_abstract import ModelContainer


class ModelContainerMultiNest(ModelContainer):

    include_priors = False
    polychord_parameters = {'nlive_mult': 25,
                            'num_repeats_mult': 5,
                            'feedback': 1,
                            'precision_criterion': 0.001,
                            'sampling_efficiency': 0.8,
                            'max_ndead': -1,
                            'boost_posterior': 0.0,
                            'read_resume': True,
                            'base_dir': './',
                            'shutdown_jitter': False}
    polychord_dir_output = None


    def multinest_priors(self, cube, ndim, nparams):
        # cube[:] = (self.bounds[:, 1] - self.bounds[:, 0]) * cube[:] + self.bounds[:, 0]
        for i in xrange(0, ndim):
            cube[i] = (self.bounds[i, 1] - self.bounds[i, 0]) * cube[i] + self.bounds[i, 0]

    def multinest_call(self, theta1, ndim, nparams):
        # Workaround for variable selection: if a variable as null index
        # (i.e. it has not been included in the model)
        # the numpy array will give back an empty list, the ctype will give back an error
        theta = np.empty(ndim)
        for i in xrange(0, ndim):
            theta[i] = theta1[i]
        chi_out = self(theta)
        if chi_out < -0.5e10:
            return -0.5e10
        return chi_out