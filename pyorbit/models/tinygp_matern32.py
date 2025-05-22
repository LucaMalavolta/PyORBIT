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

    def _build_tinygp_matern32(params):
        kernel = jnp.power(params['amplitude'], 2.0) * kernels.Matern32(scale=jnp.abs(params["scale"]))

        return GaussianProcess(
            kernel, params['x0'], diag=jnp.abs(params['diag']), mean=0.0
        )

    @jax.jit
    def _loss_tinygp_matern32(params):
        gp = _build_tinygp_matern32(params)
        return gp.log_probability(params['y'])

except:
    pass

__all__ = ['TinyGaussianProcess_Matern32']



class TinyGaussianProcess_Matern32(AbstractModel, AbstractGaussianProcesses):
    '''
    - matern32_scale: the scale of the Matern32 kernel;
    - matern32_sigma: the amplitude of the correlations;
    '''

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

        self.list_pams_dataset = OrderedSet([
            'matern32_scale',  # time scale of the Matern32
            'matern32_sigma', # amplitude of the covariance matrix
        ])


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

    def lnlk_compute(self, parameter_values, dataset):

        if 'matern32_rho' in parameter_values:
            parameter_values['matern32_scale'] = parameter_values['matern32_rho']

        self.update_parameter_values(parameter_values,
                                        replace_rotation='matern32_scale',
                                        replace_decay='matern32_scale')

        theta_dict =  dict(
            scale=parameter_values['matern32_scale'],
            amplitude=parameter_values['matern32_sigma'],
            diag=dataset.e ** 2.0 + dataset.jitter ** 2.0,
            x0=dataset.x0,
            y=dataset.residuals
        )
        return _loss_tinygp_matern32(theta_dict)


    def sample_predict(self, parameter_values, dataset, x0_input=None, return_covariance=False, return_variance=False):

        if 'matern32_rho' in parameter_values:
            parameter_values['matern32_scale'] = parameter_values['matern32_rho']

        self.update_parameter_values(parameter_values, replace_rotation='matern32_scale')

        if x0_input is None:
            x0 = dataset.x0
        else:
            x0 = x0_input

        theta_dict =  dict(
            scale=parameter_values['matern32_scale'],
            amplitude=parameter_values['matern32_sigma'],
            diag=dataset.e ** 2.0 + dataset.jitter ** 2.0,
            x0=dataset.x0,
            y=dataset.residuals,
            x0_predict = x0
        )

        gp = _build_tinygp_matern32(theta_dict)
        _, cond_gp = gp.condition(theta_dict['y'], theta_dict['x0_predict'])
        mu = cond_gp.loc # or cond_gp.mean?
        std = np.sqrt(cond_gp.variance)

        if return_variance:
            return mu, std
        else:
            return mu

