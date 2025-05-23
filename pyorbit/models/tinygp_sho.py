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

        self.list_pams_dataset = OrderedSet([
            'sho_scale',
            'sho_decay',
            'sho_sigma',  # sigma
        ])

        self.use_gp_notation = False

        """ semi-separable matrix requires the input array to be ordered
            in increasing order of the x0 values.
        """
        self._sorting_mask = {}

    def initialize_model(self, mc,  **kwargs):

        self._prepare_hyperparameter_conditions(mc, **kwargs)

        self._prepare_shared_hyperparameters(pam_scale='sho_scale', pam_decay='sho_decay', **kwargs)

        self._prepare_rotation_replacement(mc,
                                            parameter_name='sho_scale',
                                            common_pam=self.use_shared_scale,
                                            check_common=False,
                                            **kwargs)
        self._prepare_decay_replacement(mc,
                                            parameter_name='sho_decay',
                                            common_pam=self.use_shared_scale,
                                            check_common=False,
                                            **kwargs)

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        """ when reloading the .p files, the object is not reinitialized, so we have to skip the
        incremental addition of datasets if they are already present  """
        if dataset.name_ref in self._sorting_mask:
            return

        self._sorting_mask[dataset.name_ref] = np.argsort(dataset.x0)

    def lnlk_compute(self, parameter_values, dataset):

        self.update_parameter_values(parameter_values,
                                        replace_rotation='sho_scale',
                                        replace_decay='sho_decay')

        pass_conditions = self.check_hyperparameter_values(parameter_values,
                                        pam_scale='sho_scale',
                                        pam_decay='sho_decay')
        if not pass_conditions:
            return pass_conditions

        omega = 2 * np.pi / parameter_values['sho_scale']
        quality_factor =  omega * parameter_values['sho_decay'] / 2.

        dataset_ej2 = dataset.e ** 2.0 + dataset.jitter ** 2.0

        theta_dict =  dict(
            omega=omega,
            quality=quality_factor,
            sigma=parameter_values['sho_sigma'],
            diag=dataset_ej2[self._sorting_mask[dataset.name_ref]],
            x0=dataset.x0[self._sorting_mask[dataset.name_ref]],
            y=dataset.residuals[self._sorting_mask[dataset.name_ref]]
        )
        return _loss_tinygp_sho(theta_dict)


    def sample_predict(self, parameter_values, dataset, x0_input=None, return_covariance=False, return_variance=False):

        self.update_parameter_values(parameter_values,
                                        replace_rotation='sho_scale',
                                        replace_decay='sho_decay')

        omega = 2 * np.pi / parameter_values['sho_scale']
        quality_factor =  omega * parameter_values['sho_decay'] / 2.

        if x0_input is None:
            x0 = dataset.x0
            sorting_predict = self._sorting_mask[dataset.name_ref]
        else:
            x0 = x0_input
            sorting_predict = np.argsort(x0)

        dataset_ej2 = dataset.e ** 2.0 + dataset.jitter ** 2.0

        theta_dict =  dict(
            omega=omega,
            quality=quality_factor,
            sigma=parameter_values['sho_sigma'],
            diag=dataset_ej2[self._sorting_mask[dataset.name_ref]],
            x0=dataset.x0[self._sorting_mask[dataset.name_ref]],
            y=dataset.residuals[self._sorting_mask[dataset.name_ref]],
            x0_predict = x0[sorting_predict]
        )

        mu = np.empty_like(sorting_predict)
        std = np.empty_like(sorting_predict)

        gp = _build_tinygp_sho(theta_dict)
        _, cond_gp = gp.condition(theta_dict['y'], theta_dict['x0_predict'])
        mu[sorting_predict] = cond_gp.loc # or cond_gp.mean?
        std[sorting_predict] = np.sqrt(cond_gp.variance)

        if return_variance:
            return mu, std
        else:
            return mu
