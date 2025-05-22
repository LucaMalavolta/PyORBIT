from pyorbit.subroutines.common import np, OrderedSet
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_gaussian_processes import AbstractGaussianProcesses
from pyorbit.keywords_definitions import *

try:
    import celerite2
except (ModuleNotFoundError,ImportError):
    pass


class Celerite2_SHO(AbstractModel, AbstractGaussianProcesses):

    r"""A term representing a stochastically-driven, damped harmonic oscillator

    Args:
       sho_period (float) - (rho) the undamped period of the oscillator
       sho_tau (float) - the damping timescale of the process,
       sho_sigma (Optional[float]) - the standard deviation of the process
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

        self.list_pams_dataset = OrderedSet([
            'sho_scale',
            'sho_decay',
            'sho_sigma',  # sigma
        ])

        self.n_pams = 3
        self.gp = {}
        self.use_gp_notation = False


    def initialize_model(self, mc, **kwargs):

        self._prepare_hyperparameter_conditions(mc, **kwargs)

        self._prepare_shared_hyperparameters(pam_scale='sho_scale', pam_decay='sho_decay', **kwargs)

        self._prepare_rotation_replacement(mc,
                                            parameter_name='sho_scale',
                                            common_pam=self.use_shared_scale,
                                            check_common=False,
                                            **kwargs)
        self._prepare_decay_replacement(mc,
                                            parameter_name='sho_decay',
                                            common_pam=self.use_shared_scale,
                                            check_common=False,
                                            **kwargs)

    def initialize_model_dataset(self, mc, dataset, **kwargs):
        self.define_kernel(dataset)
        return

    def define_kernel(self, dataset):
        # random initialization

        kernel = celerite2.terms.SHOTerm(rho=1.0, tau=1.0, sigma=1.0)
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
        self.update_parameter_values(parameter_values,
                                        replace_rotation='sho_scale',
                                        replace_decay='sho_decay')

        pass_conditions = self.check_hyperparameter_values(parameter_values,
                                        pam_scale='sho_scale',
                                        pam_decay='sho_decay')

        if not pass_conditions:
            return pass_conditions

        self.gp[dataset.name_ref].mean = 0.
        self.gp[dataset.name_ref].kernel = celerite2.terms.SHOTerm(rho=parameter_values['sho_scale'],
                                                                    tau=parameter_values['sho_decay'],
                                                                    sigma=parameter_values['sho_sigma'])

        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        return self.gp[dataset.name_ref].log_likelihood(dataset.residuals)

    def sample_predict(self, parameter_values, dataset, x0_input=None, return_covariance=False, return_variance=False):

        self.update_parameter_values(parameter_values,
                                        replace_rotation='sho_scale',
                                        replace_decay='sho_decay')

        self.gp[dataset.name_ref].mean = 0.
        self.gp[dataset.name_ref].kernel = celerite2.terms.SHOTerm(rho=parameter_values['sho_scale'],
                                                                    tau=parameter_values['sho_decay'],
                                                                    sigma=parameter_values['sho_sigma'])

        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        if x0_input is None:
            return self.gp[dataset.name_ref].predict(dataset.residuals, dataset.x0, return_cov=return_covariance, return_var=return_variance)
        else:
            return self.gp[dataset.name_ref].predict(dataset.residuals, x0_input, return_cov=return_covariance, return_var=return_variance)

    def sample_conditional(self, parameter_values, dataset,  x0_input=None):

        self.update_parameter_values(parameter_values,
                                        replace_rotation='sho_scale',
                                        replace_decay='sho_decay')

        self.gp[dataset.name_ref].mean = 0.
        self.gp[dataset.name_ref].kernel = celerite2.terms.SHOTerm(rho=parameter_values['sho_scale'],
                                                                    tau=parameter_values['sho_decay'],
                                                                    sigma=parameter_values['sho_sigma'])

        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        if x0_input is None:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, dataset.x0)
        else:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, x0_input)

