from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *

try:
    import jax
    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp
    import tinygp
except:
    pass


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
            import tinygp
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

        self.n_pams = 4

        self.gp_pams_index = {
            'Hamp': 0,  # amp2
            'Pdec': 1,  # metric
            'Oamp': 2,  # gamma
            'Prot': 3  # ln_P
        }

        self.gp = {}


    @staticmethod
    def _build_tinygp_quasiperiodic(variable_dict, diag, X):

        gamma = jnp.abs(1. / (2.*variable_dict['Oamp'] ** 2))
        kernel = jnp.power(variable_dict['Hamp'], 2.0) \
            * tinygp.kernels.ExpSquared(scale=jnp.abs(variable_dict["Pdec"])) \
                * tinygp.kernels.ExpSineSquared(
                scale=jnp.abs(variable_dict["Prot"]),
                gamma=jnp.abs(1. / (2.*variable_dict['Oamp'] ** 2))),

        return tinygp.GaussianProcess(
            kernel, X, diag=jnp.abs(diag), mean=0.0
        )


    def lnlk_compute(self, variable_value, dataset):
        """ 2 steps:
           1) theta parameters must be converted in physical units (e.g. from logarithmic to linear spaces)
           2) physical values must be converted to {\tt george} input parameters
        """
        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref] = self._build_tinygp_quasiperiodic(variable_value, diag, dataset.x0)
        return self.gp[dataset.name_ref].log_likelihood(dataset.residuals)

    def sample_predict(self, variable_value, dataset, x0_input=None, return_covariance=False, return_variance=False):

        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref] = self._build_tinygp_quasiperiodic(variable_value, diag, dataset.x0)

        if x0_input is None:
            x0 = dataset.x0
        else:
            x0 = x0_input
        _, cond_gp = self.gp[dataset].name_ref.condition(dataset.residuals, x0)
        mu = cond_gp.mean
        std = np.sqrt(cond_gp.variance)
        if return_variance:
            return mu, std
        else:
            return mu
