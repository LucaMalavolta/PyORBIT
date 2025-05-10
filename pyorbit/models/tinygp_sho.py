from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_gaussian_processes import AbstractGaussianProcesses
from pyorbit.keywords_definitions import *
import sys

try:
    import jax
    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp
    from tinygp import kernels, GaussianProcess

    if sys.version_info[1] < 10:
        raise Warning("You should be using Python 3.10 - tinygp may not work")

    def _build_tinygp_sho(params):
        kernel = kernels.quasisep.SHO(omega=jnp.abs(params["omega"]),
                                      quality=jnp.abs(params["quality"]),
                                      sigma=jnp.abs(params["sigma"]))

        return GaussianProcess(
            kernel, params['x0'], diag=jnp.abs(params['diag']), mean=0.0
        )

    @jax.jit
    def _loss_tinygp_sho(params):
        gp = _build_tinygp_sho(params)
        return gp.log_probability(params['y'])

except:
    pass

__all__ = ['TinyGaussianProcess_SHO']



class TinyGaussianProcess_SHO(AbstractModel, AbstractGaussianProcesses):

    default_common = 'activity'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            import jax
            jax.config.update("jax_enable_x64", True)
        except (ModuleNotFoundError,ImportError):
            print("ERROR: tinygp or jax not installed, this will not work")
            quit()

        self.model_class = 'gaussian_process'
        self.internal_likelihood = True

        self.list_pams_common = OrderedSet([
            'sho_period',
            'sho_tau'
        ])

        self.list_pams_dataset = OrderedSet([
            'sho_sigma',  # sigma
        ])

        self.use_gp_notation = False

    def initialize_model(self, mc,  **kwargs):

        self._prepare_hyperparameter_conditions(mc, **kwargs)

        self.retrieve_rho = self._internal_transformation_period_mod00
        self.retrieve_tau = self._internal_transformation_decay_mod00

        for dict_name in keywords_change_variable_names:
            if kwargs.get(dict_name, False):

                self.use_gp_notation = True

                self.list_pams_common.update(['Pdec'])
                self.list_pams_common.discard('sho_tau')

                self.list_pams_common.update(['Prot'])
                self.list_pams_common.discard('sho_period')

                self.retrieve_rho = self._internal_transformation_period_mod01
                self.retrieve_tau = self._internal_transformation_decay_mod01

                self._prepare_rotation_replacement(mc, **kwargs)
                self._prepare_decay_replacement(mc, **kwargs)


        if not self.use_gp_notation:
            self._prepare_rotation_replacement(mc, **kwargs)
            self._prepare_decay_replacement(mc, **kwargs)

        if self.use_stellar_rotation_period:
            self.retrieve_rho = self._internal_transformation_period_mod02

        if self.use_stellar_activity_decay:
            self.retrieve_tau = self._internal_transformation_decay_mod02

        for dict_name in keywords_shared_hyperparameters:
            if kwargs.get(dict_name, False):
                pams_copy = self.list_pams_dataset.copy()
                for pam in pams_copy:
                    self.list_pams_common.update([pam])
                    self.list_pams_dataset.discard(pam)


    def lnlk_compute(self, parameter_values, dataset):

        parameter_values['Prot'] = self.retrieve_rho(parameter_values)
        parameter_values['Pdec'] = self.retrieve_tau(parameter_values)
        pass_conditions = self.check_hyperparameter_values(parameter_values)
        if not pass_conditions:
            return pass_conditions

        omega = 2 * np.pi / parameter_values['Prot']
        quality_factor =  omega * parameter_values['Pdec'] / 2.

        theta_dict =  dict(
            omega=omega,
            quality=quality_factor,
            sigma=parameter_values['sho_sigma'],
            diag=dataset.e ** 2.0 + dataset.jitter ** 2.0,
            x0=dataset.x0,
            y=dataset.residuals
        )
        return _loss_tinygp_sho(theta_dict)


    def sample_predict(self, parameter_values, dataset, x0_input=None, return_covariance=False, return_variance=False):

        parameter_values['Prot'] = self.retrieve_rho(parameter_values)
        parameter_values['Pdec'] = self.retrieve_tau(parameter_values)
        pass_conditions = self.check_hyperparameter_values(parameter_values)
        if not pass_conditions:
            return pass_conditions

        omega = 2 * np.pi / parameter_values['Prot']
        quality_factor =  omega * parameter_values['Pdec'] / 2.


        if x0_input is None:
            x0 = dataset.x0
        else:
            x0 = x0_input

        theta_dict =  dict(
            omega=omega,
            quality=quality_factor,
            sigma=parameter_values['sho_sigma'],
            diag=dataset.e ** 2.0 + dataset.jitter ** 2.0,
            x0=dataset.x0,
            y=dataset.residuals,
            x0_predict = x0
        )

        gp = _build_tinygp_sho(theta_dict)
        _, cond_gp = gp.condition(theta_dict['y'], theta_dict['x0_predict'])
        mu = cond_gp.loc # or cond_gp.mean?
        std = np.sqrt(cond_gp.variance)

        if return_variance:
            return mu, std
        else:
            return mu

    @staticmethod
    def _internal_transformation_period_mod00(parameter_values):
        return  parameter_values['sho_period']

    @staticmethod
    def _internal_transformation_period_mod01(parameter_values):
        return  parameter_values['Prot']

    @staticmethod
    def _internal_transformation_period_mod02(parameter_values):
        return  parameter_values['rotation_period']

    @staticmethod
    def _internal_transformation_decay_mod00(parameter_values):
        return  parameter_values['sho_tau']

    @staticmethod
    def _internal_transformation_decay_mod01(parameter_values):
        return  parameter_values['Pdec']

    @staticmethod
    def _internal_transformation_decay_mod02(parameter_values):
        return  parameter_values['activity_decay']