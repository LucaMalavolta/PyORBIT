from pyorbit.subroutines.common import np, OrderedSet
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_gaussian_processes import AbstractGaussianProcesses
from pyorbit.keywords_definitions import *

try:
    import celerite2
except (ModuleNotFoundError,ImportError):
    pass


class Celerite2_Matern32(AbstractModel, AbstractGaussianProcesses):

    r"""A term that approximates a Matern-3/2 function

    Args:
       log_sigma (float) – The log of the parameter σ.
       log_rho (float) – The log of the parameter ρ.
       eps (Optional[float]) – The value of the parameter ϵ. (default: 0.01)
    """

    default_common = 'activity'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            import celerite2
        except (ModuleNotFoundError,ImportError):
            print("ERROR: celerite2 not installed, this will not work")
            quit()

        self.model_class = 'gaussian_process'
        self.internal_likelihood = True

        self.list_pams_common = OrderedSet()

        self.list_pams_dataset = OrderedSet([
            'matern32_scale',  # time scale of the Matern32
            'matern32_sigma', # amplitude of the covariance matrix
        ])

        self.n_pams = 2
        self.gp = {}

    def initialize_model(self, mc,  **kwargs):

        self._prepare_shared_hyperparameters(pam_scale='matern32_scale', pam_decay='matern32_scale', **kwargs)

        self._prepare_rotation_replacement(mc,
                                            parameter_name='matern32_scale',
                                            common_pam=self.use_shared_scale,
                                            check_common=False,
                                            **kwargs)
        self._prepare_decay_replacement(mc,
                                            parameter_name='matern32_scale',
                                            common_pam=self.use_shared_scale,
                                            check_common=False,
                                            **kwargs)

        self._check_extra_conditions(**kwargs)

    def initialize_model_dataset(self, mc, dataset, **kwargs):
        self.define_kernel(dataset)
        return

    def define_kernel(self, dataset):
        # random initialization

        kernel = celerite2.terms.Matern32Term(sigma=0.0, rho=0.0)
        self.gp[dataset.name_ref] = celerite2.GaussianProcess(kernel, mean=0.0)

        """ I've decided to add the jitter in quadrature instead of using a constant kernel to allow the use of
        different / selective jitter within the dataset
        """
        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag)
        return

    def lnlk_compute(self, parameter_values, dataset):

        """ 2 steps:
        In celerite2 the old function "set_parameter_vector" has been removed
        and the kernel has to be defined every time
        """

        if 'matern32_rho' in parameter_values:
            parameter_values['matern32_scale'] = parameter_values['matern32_rho']

        self.update_parameter_values(parameter_values,
                                        replace_rotation='matern32_scale',
                                        replace_decay='matern32_scale')

        self.gp[dataset.name_ref].mean = 0.
        self.gp[dataset.name_ref].kernel = celerite2.terms.Matern32Term(
            sigma=parameter_values['matern32_sigma'],
            rho=parameter_values['matern32_scale'])

        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        return self.gp[dataset.name_ref].log_likelihood(dataset.residuals)

    def sample_predict(self, parameter_values, dataset, x0_input=None, return_covariance=False, return_variance=False):

        if 'matern32_rho' in parameter_values:
            parameter_values['matern32_scale'] = parameter_values['matern32_rho']

        self.update_parameter_values(parameter_values,
                                        replace_rotation='matern32_scale',
                                        replace_decay='matern32_scale')

        self.gp[dataset.name_ref].mean = 0.
        self.gp[dataset.name_ref].kernel = celerite2.terms.Matern32Term(
            sigma=parameter_values['matern32_sigma'],
            rho=parameter_values['matern32_scale'])

        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        if x0_input is None:
            return self.gp[dataset.name_ref].predict(dataset.residuals, dataset.x0, return_cov=return_covariance, return_var=return_variance)
        else:
            return self.gp[dataset.name_ref].predict(dataset.residuals, x0_input, return_cov=return_covariance, return_var=return_variance)

    def sample_conditional(self, parameter_values, dataset,  x0_input=None):

        if 'matern32_rho' in parameter_values:
            parameter_values['matern32_scale'] = parameter_values['matern32_rho']

        self.update_parameter_values(parameter_values,
                                        replace_rotation='matern32_scale',
                                        replace_decay='matern32_scale')

        self.gp[dataset.name_ref].mean = 0.
        self.gp[dataset.name_ref].kernel = celerite2.terms.Matern32Term(
            sigma=parameter_values['matern32_sigma'],
            rho=parameter_values['matern32_scale'])

        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        if x0_input is None:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, dataset.x0)
        else:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, x0_input)
