from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
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

__all__ = ['TinyGaussianProcess_Matern32Activity']



class TinyGaussianProcess_Matern32Activity(AbstractModel):
    '''
    - matern32_rho: the scale of the Matern32 kernel;
    - matern32_sigma: the amplitude of the correlations;
    - matern32_sigma_deriv: amplitude of the first derivative
    '''

    default_common = 'activity'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        try:
            import jax
            jax.config.update("jax_enable_x64", True)
        except (ModuleNotFoundError,ImportError):
            print("ERROR: tinygp or jax not installed, this will not work")
            quit()

        self.model_class = 'tinygp_matern32'
        self.internal_likelihood = True

        self.list_pams_dataset = OrderedSet([
            'matern32_rho',  # time scale of the Matern32
            'matern32_sigma', # Amplitude of the covariance matrix
        ])
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


    def lnlk_compute(self, parameter_values, dataset):
        if self.use_stellar_rotation_period:
            parameter_values['matern32_rho'] = parameter_values['rotation_period']

        theta_dict =  dict(
            scale=parameter_values['matern32_rho'],
            amplitude=parameter_values['matern32_sigma'],
            diag=dataset.e ** 2.0 + dataset.jitter ** 2.0,
            x0=dataset.x0,
            y=dataset.residuals
        )
        return _loss_tinygp_matern32(theta_dict)


    def sample_predict(self, parameter_values, dataset, x0_input=None, return_covariance=False, return_variance=False):

        if self.use_stellar_rotation_period:
            parameter_values['matern32_rho'] = parameter_values['rotation_period']

        if x0_input is None:
            x0 = dataset.x0
        else:
            x0 = x0_input

        theta_dict =  dict(
            scale=parameter_values['matern32_rho'],
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

