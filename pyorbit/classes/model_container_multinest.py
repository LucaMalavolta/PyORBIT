from pyorbit.subroutines.common import *
from pyorbit.classes.model_container_abstract import ModelContainer


class ModelContainerMultiNest(ModelContainer):

    def __init__(self):
        super(self.__class__, self).__init__()

        self.include_priors = False
        self.nested_sampling_parameters = {'nlive_mult': 25,
                                           'base_dir': './',
                                           'verbose': True,
                                           'sampling_efficiency': 0.3,
                                           'shutdown_jitter': False,
                                           'include_priors': False}

        self.pymultinest_signature = [
            'n_params',
            'n_clustering_params',
            'wrapped_params',
            'importance_nested_sampling',
            'multimodal',
            'const_efficiency_mode',
            'n_live_points',
            'evidence_tolerance',
            'sampling_efficiency',
            'n_iter_before_update',
            'null_log_evidence',
            'max_modes',
            'mode_tolerance',
            'outputfiles_basename',
            'seed',
            'verbose',
            'resume',
            'context',
            'write_output',
            'log_zero',
            'max_iter',
            'init_MPI',
            'dump_callback']

        self.output_directory = None

    def multinest_priors(self, cube, ndim, nparams):

        for i in range(0, ndim):
            cube[i] = nested_sampling_prior_compute(cube[i], self.priors[i][0], self.priors[i][2], self.spaces[i])


    def multinest_call(self, theta1, ndim, nparams):
        # Workaround for parameter selection: if a parameter as null index
        # (i.e. it has not been included in the model)
        # the numpy array will give back an empty list, the ctype will give back an error
        theta = np.empty(ndim)
        for i in range(0, ndim):
            theta[i] = theta1[i]

        chi_out = self(theta, self.include_priors)

        if chi_out < -0.5e30:
            return -0.5e30

        return chi_out
