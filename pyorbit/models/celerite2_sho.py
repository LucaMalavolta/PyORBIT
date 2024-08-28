from pyorbit.subroutines.common import np, OrderedSet
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.keywords_definitions import *

try:
    import celerite2
except (ModuleNotFoundError,ImportError):
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
        except (ModuleNotFoundError,ImportError):
            print("ERROR: celerite2 not installed, this will not work")
            quit()

        self.model_class = 'celerite2_sho'
        self.internal_likelihood = True

        self.list_pams_common = OrderedSet([
            'sho_period',
            'sho_tau'
        ])

        self.list_pams_dataset = OrderedSet([
            'sho_sigma',  # sigma
        ])

        self.n_pams = 3
        self.gp = {}


    def initialize_model(self, mc, **kwargs):

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == 'activity':
                self.use_stellar_rotation_period = getattr(mc.common_models[common_ref], 'use_stellar_rotation_period', False)
                break

        for keyword in keywords_stellar_rotation:
            self.use_stellar_rotation_period = kwargs.get(keyword, self.use_stellar_rotation_period)

        self.retrieve_rho_tau = self._internal_transformation_mod00


        if self.use_stellar_rotation_period:
            self.list_pams_common.update(['rotation_period'])
            self.list_pams_common.discard('sho_period')
            self.retrieve_rho_tau = self._internal_transformation_mod02

        for dict_name in keywords_change_variable_names:
            if kwargs.get(dict_name, False):

                self.list_pams_common.update(['Pdec'])
                self.list_pams_common.discard('sho_tau')

                if self.use_stellar_rotation_period:
                    self.list_pams_common.update(['rotation_period'])
                    self.retrieve_rho_tau = self._internal_transformation_mod04

                else:
                    self.list_pams_common.update(['Prot'])
                    self.list_pams_common.discard('sho_period')
                    self.retrieve_rho_tau = self._internal_transformation_mod01


        self.use_shared_hyperparameters = False
        for keyword in ['use_shared_hyperparameters',
                        'shared_hyperparameters',
                        'use_common_hyperparameters',
                        'common_hyperparameters']:
            self.use_shared_hyperparameters =  kwargs.get(keyword, self.use_shared_hyperparameters)
        if self.use_shared_hyperparameters:
            pams_copy = self.list_pams_dataset.copy()
            for pam in pams_copy:
                self.list_pams_common.update([pam])
                self.list_pams_dataset.discard(pam)


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

    @staticmethod
    def _internal_transformation_mod02(parameter_values):
        return  parameter_values['rotation_period'], parameter_values['sho_tau']

    @staticmethod
    def _internal_transformation_mod03(parameter_values):
        return  parameter_values['rotation_period'], parameter_values['Pdec']