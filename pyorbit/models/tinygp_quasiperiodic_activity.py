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


    @jax.jit
    def _loss_tinygp(params):
        gp = _build_tinygp_quasiperiodic(params)
        return gp.log_probability(params['y'])

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


class TinyGaussianProcess_QuasiPeriodicActivity(AbstractModel):
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

        try:
            import tinygp
            import jax
            jax.config.update("jax_enable_x64", True)
        except (ModuleNotFoundError,ImportError):
            print("ERROR: tinygp or jax not installed, this will not work")
            quit()

        self.model_class = 'gp_quasiperiodic'
        self.internal_likelihood = True

        self.list_pams_common = OrderedSet([
            'Prot',  # Rotational period of the star
            'Pdec',  # Decay timescale of activity
            'Oamp'  # Granulation of activity
        ])

        self.list_pams_dataset = OrderedSet([
            'Hamp'  # Amplitude of the signal in the covariance matrix
        ])

    def initialize_model(self, mc,  **kwargs):

        if kwargs.get('hyperparameters_condition', False):
            self.hyper_condition = self._hypercond_01
        else:
            self.hyper_condition = self._hypercond_00

        if kwargs.get('rotation_decay_condition', False):
            self.rotdec_condition = self._hypercond_02
        else:
            self.rotdec_condition = self._hypercond_00

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == 'activity':
                self.use_stellar_rotation_period = getattr(mc.common_models[common_ref], 'use_stellar_rotation_period', False)
                break

        for keyword in keywords_stellar_rotation:
            self.use_stellar_rotation_period = kwargs.get(keyword, self.use_stellar_rotation_period)

        if self.use_stellar_rotation_period:
            self.list_pams_common.update(['rotation_period'])
            self.list_pams_common.discard('Prot')

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == 'activity':
                self.use_stellar_activity_decay = getattr(mc.common_models[common_ref], 'use_stellar_activity_decay', False)
                break

        for keyword in keywords_stellar_activity_decay:
            self.use_stellar_activity_decay = kwargs.get(keyword, self.use_stellar_activity_decay)

        if self.use_stellar_activity_decay:
            self.list_pams_common.update(['activity_decay'])
            self.list_pams_common.discard('Pdec')

    def lnlk_compute(self, parameter_values, dataset):

        if self.use_stellar_rotation_period:
            parameter_values['Prot'] = parameter_values['rotation_period']

        if self.use_stellar_activity_decay:
            parameter_values['Pdec'] = parameter_values['activity_decay']

        if not self.hyper_condition(parameter_values):
            return -np.inf
        if not self.rotdec_condition(parameter_values):
            return -np.inf

        theta_dict =  dict(
            gamma=1. / (2.*parameter_values['Oamp'] ** 2),
            Hamp=parameter_values['Hamp'],
            Pdec=parameter_values['Pdec'],
            Prot=parameter_values['Prot'],
            diag=dataset.e ** 2.0 + dataset.jitter ** 2.0,
            x0=dataset.x0,
            y=dataset.residuals
        )
        return _loss_tinygp(theta_dict)


    def sample_predict(self, parameter_values, dataset, x0_input=None, return_covariance=False, return_variance=False):

        if self.use_stellar_rotation_period:
            parameter_values['Prot'] = parameter_values['rotation_period']

        if self.use_stellar_activity_decay:
            parameter_values['Pdec'] = parameter_values['activity_decay']


        if x0_input is None:
            x0 = dataset.x0
        else:
            x0 = x0_input

        theta_dict =  dict(
            gamma=1. / (2.*parameter_values['Oamp'] ** 2),
            Hamp=parameter_values['Hamp'],
            Pdec=parameter_values['Pdec'],
            Prot=parameter_values['Prot'],
            diag=dataset.e ** 2.0 + dataset.jitter ** 2.0,
            x0=dataset.x0,
            y=dataset.residuals,
            x0_predict = x0
        )

        gp = _build_tinygp_quasiperiodic(theta_dict)
        _, cond_gp = gp.condition(theta_dict['y'], theta_dict['x0_predict'])
        mu = cond_gp.loc # or cond_gp.mean?
        std = np.sqrt(cond_gp.variance)

        if return_variance:
            return mu, std
        else:
            return mu

    @staticmethod
    def _hypercond_00(parameter_values):
        #Condition from Rajpaul 2017, Rajpaul+2021
        return True

    @staticmethod
    def _hypercond_01(parameter_values):
        # Condition from Rajpaul 2017, Rajpaul+2021
        # Taking into account that Pdec^2 = 2*lambda_2^2
        return parameter_values['Pdec']**2 > (3. / 2. / np.pi) * parameter_values['Oamp']**2 * parameter_values['Prot']**2

    @staticmethod
    def _hypercond_02(parameter_values):
        #Condition on Rotation period and decay timescale
        return parameter_values['Pdec'] > 2. * parameter_values['Prot']

