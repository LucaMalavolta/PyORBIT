from pyorbit.subroutines.common import np, OrderedSet
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.keywords_definitions import *

try:
    import celerite2
except (ModuleNotFoundError,ImportError):
    pass


class Celerite2_Matern32(AbstractModel):

    r"""A term that approximates a Matern-3/2 function

    Args:
       log_sigma (float) – The log of the parameter σ.
       log_rho (float) – The log of the parameter ρ.
       eps (Optional[float]) – The value of the parameter ϵ. (default: 0.01)
    """

    default_common = 'activity'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        try:
            import celerite2
        except (ModuleNotFoundError,ImportError):
            print("ERROR: celerite2 not installed, this will not work")
            quit()

        self.model_class = 'celerite2_matern32'
        self.internal_likelihood = True

        self.list_pams_common = OrderedSet()

        self.list_pams_dataset = OrderedSet([
            'matern32_sigma',  # sigma
            'matern32_rho',  # rho
        ])

        self.n_pams = 2
        self.gp = {}
        self.use_stellar_rotation_period = False

    def initialize_model(self, mc,  **kwargs):

        try:
            for common_ref in self.common_ref:
                if mc.common_models[common_ref].model_class == 'activity':
                    self.use_stellar_rotation_period = getattr(mc.common_models[common_ref], 'use_stellar_rotation_period', False)
                    break
        except:
            self.use_stellar_rotation_period = False

        for keyword in keywords_stellar_rotation:
            self.use_stellar_rotation_period = kwargs.get(keyword, self.use_stellar_rotation_period)

        if self.use_stellar_rotation_period:
            self.list_pams_common.update(['rotation_period'])
            self.list_pams_dataset.discard('matern32_rho')

        self.use_shared_hyperparameters = False
        for keyword in keywords_shared_hyperparameters:
            self.use_shared_hyperparameters =  kwargs.get(keyword, self.use_shared_hyperparameters)
        if self.use_shared_hyperparameters:
            pams_copy = self.list_pams_dataset.copy()
            for pam in pams_copy:
                self.list_pams_common.update([pam])
                self.list_pams_dataset.discard(pam)

        self.use_shared_rho = False
        for keyword in keywords_shared_timescale:
            self.use_shared_rho =  kwargs.get(keyword, self.use_shared_rho)
        if self.use_shared_rho and not self.use_stellar_rotation_period:
            pam = 'matern32_rho'
            self.list_pams_common.update([pam])
            self.list_pams_dataset.discard(pam)

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
        if self.use_stellar_rotation_period:
            parameter_values['matern32_rho'] = parameter_values['rotation_period']

        self.gp[dataset.name_ref].mean = 0.
        self.gp[dataset.name_ref].kernel = celerite2.terms.Matern32Term(
            sigma=parameter_values['matern32_sigma'],
            rho=parameter_values['matern32_rho'])

        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        return self.gp[dataset.name_ref].log_likelihood(dataset.residuals)

    def sample_predict(self, parameter_values, dataset, x0_input=None, return_covariance=False, return_variance=False):

        if self.use_stellar_rotation_period:
            parameter_values['matern32_rho'] = parameter_values['rotation_period']

        self.gp[dataset.name_ref].mean = 0.
        self.gp[dataset.name_ref].kernel = celerite2.terms.Matern32Term(
            sigma=parameter_values['matern32_sigma'],
            rho=parameter_values['matern32_rho'])

        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        if x0_input is None:
            return self.gp[dataset.name_ref].predict(dataset.residuals, dataset.x0, return_cov=return_covariance, return_var=return_variance)
        else:
            return self.gp[dataset.name_ref].predict(dataset.residuals, x0_input, return_cov=return_covariance, return_var=return_variance)

    def sample_conditional(self, parameter_values, dataset,  x0_input=None):

        if self.use_stellar_rotation_period:
            parameter_values['matern32_rho'] = parameter_values['rotation_period']

        self.gp[dataset.name_ref].mean = 0.
        self.gp[dataset.name_ref].kernel = celerite2.terms.Matern32Term(
            sigma=parameter_values['matern32_sigma'],
            rho=parameter_values['matern32_rho'])

        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        if x0_input is None:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, dataset.x0)
        else:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, x0_input)
