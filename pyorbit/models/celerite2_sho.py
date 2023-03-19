from pyorbit.subroutines.common import np
from pyorbit.models.abstract_model import AbstractModel

try:
    import celerite2
except ImportError:
    pass


class Celerite2_SHO(AbstractModel):

    r"""A term representing a stochastically-driven, damped harmonic oscillator

    Args:
       sho_period (float) - (rho) the undamped period of the oscillator
       sho_tau (float) - the damping timescale of the process,
       sho_sigma (Optional[float]) - the standard deviation of the process
    """

    default_common = 'activity'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        try:
            import celerite2
        except:
            print("ERROR: celerite2 not installed, this will not work")
            quit()

        self.model_class = 'celerite2_sho'
        self.internal_likelihood = True

        self.list_pams_common = {
            'sho_period',
            'sho_tau'
        }

        self.list_pams_dataset = {
            'sho_sigma',  # sigma
        }

        self.n_pams = 3
        self.gp = {}


    def initialize_model(self, mc, **kwargs):

        self.retrieve_rho_tau = self._internal_transformation_mod00

        change_variable_names = [
            'use_GPquasiperiod_notation',
            'use_GPquasiperiodic_notation',
            'use_Prot_Pdec_notation',
            'use_Prot_Pdec'
        ]

        for dict_name in change_variable_names:
            if kwargs.get(dict_name, False):
                self.list_pams_common.update(['Prot'])
                self.list_pams_common.update(['Pdec'])
                self.list_pams_common.discard('sho_period')
                self.list_pams_common.discard('sho_tau')
                self.retrieve_rho_tau = self._internal_transformation_mod01


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
        rho, tau = self.retrieve_rho_tau(parameter_values)

        self.gp[dataset.name_ref].mean = 0.
        self.gp[dataset.name_ref].kernel = celerite2.terms.SHOTerm(rho=rho, tau=tau, sigma=parameter_values['sho_sigma'])

        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        return self.gp[dataset.name_ref].log_likelihood(dataset.residuals)

    def sample_predict(self, parameter_values, dataset, x0_input=None, return_covariance=False, return_variance=False):

        rho, tau = self.retrieve_rho_tau(parameter_values)

        self.gp[dataset.name_ref].mean = 0.
        self.gp[dataset.name_ref].kernel = celerite2.terms.SHOTerm(rho=rho, tau=tau, sigma=parameter_values['sho_sigma'])

        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        if x0_input is None:
            return self.gp[dataset.name_ref].predict(dataset.residuals, dataset.x0, return_cov=return_covariance, return_var=return_variance)
        else:
            return self.gp[dataset.name_ref].predict(dataset.residuals, x0_input, return_cov=return_covariance, return_var=return_variance)

    def sample_conditional(self, parameter_values, dataset,  x0_input=None):

        rho, tau = self.retrieve_rho_tau(parameter_values)

        self.gp[dataset.name_ref].mean = 0.
        self.gp[dataset.name_ref].kernel = celerite2.terms.SHOTerm(rho=rho, tau=tau, sigma=parameter_values['sho_sigma'])

        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        if x0_input is None:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, dataset.x0)
        else:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, x0_input)

    @staticmethod
    def _internal_transformation_mod00(parameter_values):
        return  parameter_values['sho_period'], parameter_values['sho_tau']

    @staticmethod
    def _internal_transformation_mod01(parameter_values):
        return  parameter_values['Prot'], parameter_values['Pdec']