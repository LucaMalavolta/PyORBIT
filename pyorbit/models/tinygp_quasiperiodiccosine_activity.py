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

    def _build_tinygp_quasiperiodic_cosine(params):

        kernel = kernels.ExpSquared(scale=jnp.abs(params["Pdec"])) * (
            jnp.power(params['Hamp'], 2.0) * kernels.ExpSineSquared(
                scale=jnp.abs(params["Prot"]),
                gamma=jnp.abs(params["gamma"])) +  \
                jnp.power(params['Camp'], 2.0) * kernels.Cosine(scale= jnp.abs(params["Prot"]/2.)))

        return GaussianProcess(
            kernel, params['x0'], diag=jnp.abs(params['diag']), mean=0.0
        )

    @jax.jit
    def _loss_tinygp_QPcosine(params):
        gp = _build_tinygp_quasiperiodic_cosine(params)
        return gp.log_probability(params['y'])

except:
    pass


__all__ = ['TinyGaussianProcess_QuasiPeriodicCosineActivity']



class TinyGaussianProcess_QuasiPeriodicCosineActivity(AbstractModel, AbstractGaussianProcesses):
    ''' Three parameters out of four are the same for all the datasets, since they are related to
    the properties of the physical process rather than the observed effects on a dataset
     From Grunblatt+2015, Affer+2016
     - theta: is usually related to the rotation period of the star( or one of its harmonics);
     - lambda: is the correlation decay timescale, and it can be related to the lifetime of the active regions.
     - omega: is the length scale of the periodic component, and can be linked to the size evolution of the active regions;
     - h: represents the amplitude of the correlations '''

    default_common = 'activity'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            import tinygp
            import jax
            jax.config.update("jax_enable_x64", True)
        except (ModuleNotFoundError,ImportError):
            print("ERROR: tinygp or jax not installed, this will not work")
            quit()

        self.model_class = 'gaussian_process'
        self.internal_likelihood = True

        self.list_pams_common = OrderedSet([
            'Prot',  # Rotational period of the star
            'Pdec',  # Decay timescale of activity
            'Oamp'  # Granulation of activity
        ])

        self.list_pams_dataset = OrderedSet([
            'Hamp',  # Amplitude of the signal in the covariance matrix
            'Camp'  # Amplitude of the cycle
        ])

    def initialize_model(self, mc,  **kwargs):

        self._prepare_hyperparameter_conditions(mc, **kwargs)
        self._prepare_rotation_replacement(mc, **kwargs)
        self._prepare_decay_replacement(mc, **kwargs)

    def lnlk_compute(self, parameter_values, dataset):

        self.update_parameter_values(parameter_values)
        
        pass_conditions = self.check_hyperparameter_values(parameter_values)
        if not pass_conditions:
            return pass_conditions
        
        theta_dict =  dict(
            gamma=1. / (2.*parameter_values['Oamp'] ** 2),
            Hamp=parameter_values['Hamp'],
            Camp=parameter_values['Camp'],
            Pdec=parameter_values['Pdec'],
            Prot=parameter_values['Prot'],
            diag=dataset.e ** 2.0 + dataset.jitter ** 2.0,
            x0=dataset.x0,
            y=dataset.residuals
        )
        return _loss_tinygp_QPcosine(theta_dict)


    def sample_predict(self, parameter_values, dataset, x0_input=None, return_covariance=False, return_variance=False):

        self.update_parameter_values(parameter_values)
        
        if x0_input is None:
            x0 = dataset.x0
        else:
            x0 = x0_input

        theta_dict =  dict(
            gamma=1. / (2.*parameter_values['Oamp'] ** 2),
            Hamp=parameter_values['Hamp'],
            Camp=parameter_values['Camp'],
            Pdec=parameter_values['Pdec'],
            Prot=parameter_values['Prot'],
            diag=dataset.e ** 2.0 + dataset.jitter ** 2.0,
            x0=dataset.x0,
            y=dataset.residuals,
            x0_predict = x0
        )

        gp = _build_tinygp_quasiperiodic_cosine(theta_dict)
        _, cond_gp = gp.condition(theta_dict['y'], theta_dict['x0_predict'])
        mu = cond_gp.loc # or cond_gp.mean?
        std = np.sqrt(cond_gp.variance)

        if return_variance:
            return mu, std
        else:
            return mu

