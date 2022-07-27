from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *

try:
    import jax
    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp
    from tinygp import kernels, GaussianProcess

except:
    pass


__all__ = ['TinyGaussianProcess_QuasiPeriodicActivity']


def _build_tinygp_quasiperiodic(params):
    kernel = jnp.power(params['Hamp'], 2.0) \
        * kernels.ExpSquared(scale=jnp.abs(params["Pdec"])) \
            * kernels.ExpSineSquared(
            scale=jnp.abs(params["Prot"]),
            gamma=jnp.abs(params["gamma"]))

    return GaussianProcess(
        kernel, params['x0'], diag=jnp.abs(params['diag']), mean=0.0
    )

@jax.jit
def _loss_tinygp(params):
    gp = _build_tinygp_quasiperiodic(params)
    return gp.log_probability(params['y'])

@jax.jit
def _residuals_tinygp(params):
    gp = _build_tinygp_quasiperiodic(params)
    _, cond_gp = gp.condition(params['y'], params['x0'])
    return cond_gp

class TinyGaussianProcess_QuasiPeriodicActivity(AbstractModel):
    ''' Three parameters out of four are the same for all the datasets, since they are related to
    the properties of the physical process rather than the observed effects on a dataset
     From Grunblatt+2015, Affer+2016
     - theta: is usually related to the rotation period of the star( or one of its harmonics);
     - lambda: is the correlation decay timescale, and it can be related to the lifetime of the active regions.
     - omega: is the length scale of the periodic component, and can be linked to the size evolution of the active regions;
     - h: represents the amplitude of the correlations '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        try:
            import jax
            jax.config.update("jax_enable_x64", True)
        except ImportError:
            print("ERROR: tinygp or jax not installed, this will not work")
            quit()

        self.model_class = 'gp_quasiperiodic'
        self.internal_likelihood = True

        self.list_pams_common = {
            'Prot',  # Rotational period of the star
            'Pdec',  # Decay timescale of activity
            'Oamp'  # Granulation of activity
        }

        self.list_pams_dataset = {
            'Hamp'  # Amplitude of the signal in the covariance matrix
        }

    def initialize_model(self, mc,  **kwargs):

        if kwargs.get('hyperparameters_condition', False):
            self.hyper_condition = self._hypercond_01
        else:
            self.hyper_condition = self._hypercond_00

        if kwargs.get('rotation_decay_condition', False):
            self.rotdec_condition = self._hypercond_02
        else:
            self.rotdec_condition = self._hypercond_00

    def lnlk_compute(self, variable_value, dataset):
        if not self.hyper_condition(variable_value):
            return -np.inf
        if not self.rotdec_condition(variable_value):
            return -np.inf

        theta_dict =  dict(
            gamma=1. / (2.*variable_value['Oamp'] ** 2),
            Hamp=variable_value['Hamp'],
            Pdec=variable_value['Pdec'],
            Prot=variable_value['Prot'],
            diag=dataset.e ** 2.0 + dataset.jitter ** 2.0,
            x0=dataset.x0,
            y=dataset.residuals
        )
        return _loss_tinygp(theta_dict)


    def sample_predict(self, variable_value, dataset, x0_input=None, return_covariance=False, return_variance=False):

        if x0_input is None:
            x0 = dataset.x0
        else:
            x0 = x0_input

        theta_dict =  dict(
            gamma=1. / (2.*variable_value['Oamp'] ** 2),
            Hamp=variable_value['Hamp'],
            Pdec=variable_value['Pdec'],
            Prot=variable_value['Prot'],
            diag=dataset.e ** 2.0 + dataset.jitter ** 2.0,
            x0=x0,
            y=dataset.residuals
        )


        cond_gp = _residuals_tinygp(theta_dict)
        mu = cond_gp.mean
        std = np.sqrt(cond_gp.variance)
        if return_variance:
            return mu, std
        else:
            return mu

    @staticmethod
    def _hypercond_00(variable_value):
        #Condition from Rajpaul 2017, Rajpaul+2021
        return True

    @staticmethod
    def _hypercond_01(variable_value):
        # Condition from Rajpaul 2017, Rajpaul+2021
        # Taking into account that Pdec^2 = 2*lambda_2^2
        return variable_value['Pdec']**2 > (3. / 4. / np.pi) * variable_value['Oamp']**2 * variable_value['Prot']**2 

    @staticmethod
    def _hypercond_02(variable_value):
        #Condition on Rotation period and decay timescale
        return variable_value['Pdec'] > 2. * variable_value['Prot']
