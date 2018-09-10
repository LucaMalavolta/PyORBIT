from common import *
from model_container_abstract import ModelContainer


class ModelContainerMultiNest(ModelContainer):

    def __init__(self):
        super(self.__class__, self).__init__()

        self.include_priors = True
        self.nested_sampling_parameters = {'nlive_mult': 25,
                                           'base_dir': './',
                                           'verbose': True,
                                           'sampling_efficiency': 0.3,
                                           'shutdown_jitter': False,
                                           'include_priors': True}

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
