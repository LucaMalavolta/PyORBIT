from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *

try:
    import jax
    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp
    from tinygp import kernels, GaussianProcess

except:
    pass

def _build_tinygp_quasiperiodic(params):
    kernel = jnp.power(params['Hamp'], 2.0) \
        * kernels.ExpSquared(scale=jnp.abs(params["Pdec"])) \
            * kernels.ExpSineSquared(
            scale=jnp.abs(params["Prot"]),
            gamma=jnp.abs(params["gamma"]))

    return GaussianProcess(
        kernel, params['x0'], diag=jnp.abs(params['diag']), mean=0.0
    )


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


    @staticmethod
    @jax.jit
    def _loss_tinygp(params):
        gp = _build_tinygp_quasiperiodic(params)
        return gp.log_probability(params['y'])

    def _residuals_tinygp(self, params):
        gp = _build_tinygp_quasiperiodic(params)
        _, cond_gp = gp.condition(params['y'], params['x0'])
        return cond_gp




    def lnlk_compute(self, variable_value, dataset):
        theta_dict =  dict(
            gamma=1. / (2.*variable_value['Oamp'] ** 2),
            Hamp=variable_value['Oamp'],
            Pdec=variable_value['Oamp'],
            Prot=variable_value['Oamp'],
            diag=dataset.e ** 2.0 + dataset.jitter ** 2.0,
            x0=dataset.x0,
            y=dataset.residuals
        )
        return self._loss_tinygp(theta_dict)


    def sample_predict(self, variable_value, dataset, x0_input=None, return_covariance=False, return_variance=False):

        if x0_input is None:
            x0 = dataset.x0
        else:
            x0 = x0_input
        theta_dict =  dict(
            gamma=1. / (2.*variable_value['Oamp'] ** 2),
            Hamp=variable_value['Oamp'],
            Pdec=variable_value['Oamp'],
            Prot=variable_value['Oamp'],
            diag=dataset.e ** 2.0 + dataset.jitter ** 2.0,
            x0=x0,
            y=dataset.residuals
        )


        cond_gp = self._residuals_tinygp(theta_dict)
        mu = cond_gp.mean
        std = np.sqrt(cond_gp.variance)
        if return_variance:
            return mu, std
        else:
            return mu