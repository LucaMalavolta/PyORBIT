from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_gaussian_processes import AbstractGaussianProcesses
from pyorbit.keywords_definitions import *

from scipy.linalg import cho_factor, cho_solve, lapack, LinAlgError
from scipy import matrix, spatial
import sys


__all__ = ['TinyGP_Multidimensional_Matern32']

try:
    import jax
    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp
    from tinygp import kernels, GaussianProcess

    if sys.version_info[1] < 10:
        raise Warning("You should be using Python 3.10 - tinygp may not work")


    class LatentKernel_Multi_Matern32(kernels.Kernel):
        """A custom kernel based on Matern32

        Args:
            kernel: The kernel function describing the latent process. This can be any other
                ``tinygp`` kernel.
            coeff_prim: The primal coefficients for each class. This can be thought of as how
                much the latent process itself projects into the observations for that class.
                This should be an array with an entry for each class of observation.
            coeff_deriv: The derivative coefficients for each class. This should have the same
                shape as ``coeff_prim``.
        """

        try:
            kernel : kernels.Kernel
            coeff_prim: jax.Array | float
            coeff_deriv: jax.Array | float
        except:
            pass

        def __init__(self, kernel, coeff_prim, coeff_deriv):
            self.kernel = kernel
            self.coeff_prim, self.coeff_deriv = jnp.broadcast_arrays(
                jnp.asarray(coeff_prim), jnp.asarray(coeff_deriv)
            )

        def evaluate(self, X1, X2):
            t1, label1 = X1
            t2, label2 = X2

            # Differentiate the kernel function: the first derivative wrt x1
            Kp = jax.grad(self.kernel.evaluate, argnums=0)

            # ... and the second derivative
            Kpp = jax.grad(Kp, argnums=1)

            # Evaluate the kernel matrix and all of its relevant derivatives
            K = self.kernel.evaluate(t1, t2)
            d2K_dx1dx2 = Kpp(t1, t2)

            # For stationary kernels, these are related just by a minus sign, but we'll
            # evaluate them both separately for generality's sake
            dK_dx2 = jax.grad(self.kernel.evaluate, argnums=1)(t1, t2)
            dK_dx1 = Kp(t1, t2)

            # Extract the coefficients
            a1 = self.coeff_prim[label1]
            a2 = self.coeff_prim[label2]
            b1 = self.coeff_deriv[label1]
            b2 = self.coeff_deriv[label2]

            # Construct the matrix element
            return (
                a1 * a2 * K
                + a1 * b2 * dK_dx2
                + b1 * a2 * dK_dx1
                + b1 * b2 * d2K_dx1dx2
            )


    def _build_tinygp_multidimensional_matern32(params):

        base_kernel = kernels.Matern32(scale=jnp.abs(params["scale"]))

        kernel = LatentKernel_Multi_Matern32(base_kernel, params['coeff_prime'], params['coeff_deriv'])
        return GaussianProcess(
            kernel, params['X'], diag=jnp.abs(params['diag']), mean=0.0
        )

    @jax.jit
    def _loss_tinygp_multi_matern32(params):
        gp = _build_tinygp_multidimensional_matern32(params)
        return gp.log_probability(params['y'])

    ### NEW Addded in PyORBIT v10.10
    @jax.jit
    def _condition_tinygp_matern32(params):
        gp = _build_tinygp_multidimensional_matern32(params)
        return gp.condition(params['y'], params['X'])

except:
    pass





class TinyGP_Multidimensional_Matern32(AbstractModel, AbstractGaussianProcesses):
    '''
    - matern32_scale: the scale of the Matern32 kernel;
    - matern32_multigp_sigma: the amplitude of the correlations;
    - matern32_multigp_sigma_deriv: amplitude of the first derivative
    '''

    default_common = 'activity'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

        self.model_class = 'multidimensional_gaussian_process'

        self.internal_likelihood = True
        self.delayed_lnlk_computation = True

        self.list_pams_common = OrderedSet([
            'matern32_scale',  # time scale of the Matern32
        ])
        self.list_pams_dataset = OrderedSet([
            'matern32_multigp_sigma', # Amplitude of the covariance matrix
            'matern32_multigp_sigma_deriv' # Amplitude of the first derivative of the covariance matrix
        ])


        self.internal_parameter_values = None
        self._dist_t1 = None
        self._dist_t2 = None
        self._added_datasets = 0
        self.dataset_ordering = {}
        self.inds_cache = None

        self._dataset_x0 = []
        self._dataset_label = []
        self._dataset_e2 = []
        self._dataset_names = {}

        self._dataset_nindex = []

        self.use_derivative_dict = {}

        self.internal_coeff_prime = []
        self.internal_coeff_deriv = []

        self._dataset_ej2 = []
        self._dataset_res = []

        self._added_datasets = 0
        self._n_cov_matrix = 0

        self.pi2 = np.pi * np.pi

        self.use_shared_scale = True

        ### NEW Addded in PyORBIT v10.10
        self._rv_dataset_flag = []

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

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        """ when reloading the .p files, the object is not reinitialized, so we have to skip the
        incremental addition of datasets if they are already present  """
        if dataset.name_ref in self._dataset_names:
            return

        ## NEW Addded in PyORBIT v10.10
        if dataset.kind == 'RV':
            rv_flag = np.ones_like(dataset.x0, dtype=bool)
        else:
            rv_flag = np.zeros_like(dataset.x0, dtype=bool)
        self._rv_dataset_flag = np.append(self._rv_dataset_flag, rv_flag).astype(bool)

        self._dataset_nindex.append([self._n_cov_matrix,
                                    self._n_cov_matrix+dataset.n])

        self._dataset_x0 = np.append(self._dataset_x0, dataset.x0)
        self._dataset_label = np.append(self._dataset_label, np.zeros_like(dataset.x0, dtype=int) + self._added_datasets)
        self._dataset_e2 = np.append(self._dataset_e2, dataset.e**2)

        self._dataset_names[dataset.name_ref] = self._added_datasets
        self._n_cov_matrix += dataset.n
        self._added_datasets += 1

        self._dataset_ej2 = self._dataset_e2 * 1.
        self._dataset_res = self._dataset_e2 * 0.

        self.internal_coeff_prime = np.empty(self._added_datasets)
        self.internal_coeff_deriv = np.empty(self._added_datasets)
        self._X = (self._dataset_x0, self._dataset_label.astype(int))

        if 'derivative'in kwargs:
            use_derivative = kwargs['derivative'].get(dataset.name_ref, False)
        elif dataset.name_ref in kwargs:
            use_derivative = kwargs[dataset.name_ref].get('derivative', False)
        else:
            use_derivative = True

        if not use_derivative:
            self.fix_list[dataset.name_ref] = {'matern32_multigp_sigma_deriv': [0., 0.]}


        return

    def add_internal_dataset(self, parameter_values, dataset):

        if 'matern32_rho' in parameter_values:
            parameter_values['matern32_scale'] = parameter_values['matern32_rho']

        self.update_parameter_values(parameter_values,
                                        replace_rotation='matern32_scale',
                                        replace_decay='matern32_scale')

        self.internal_parameter_values = parameter_values

        d_ind = self._dataset_names[dataset.name_ref]
        d_nstart, d_nend = self._dataset_nindex[d_ind]

        self._dataset_ej2[d_nstart:d_nend] = self._dataset_e2[d_nstart:d_nend] + dataset.jitter**2.0
        self._dataset_res[d_nstart:d_nend] = dataset.residuals

        self.internal_coeff_prime[d_ind] = parameter_values['matern32_multigp_sigma']
        self.internal_coeff_deriv[d_ind] = parameter_values['matern32_multigp_sigma_deriv']

    def lnlk_compute(self):

        theta_dict =  dict(
            scale=self.internal_parameter_values['matern32_scale'],
            diag=self._dataset_ej2,
            X=self._X,
            y=self._dataset_res,
            coeff_prime=self.internal_coeff_prime,
            coeff_deriv=self.internal_coeff_deriv
        )

        return _loss_tinygp_multi_matern32(theta_dict)

    ### NEW Addded in PyORBIT v10.10
    def lnlk_rvonly_compute(self):

        theta_dict =  dict(
            scale=self.internal_parameter_values['matern32_scale'],
            diag=self._dataset_ej2,
            X=self._X,
            y=self._dataset_res,
            coeff_prime=self.internal_coeff_prime,
            coeff_deriv=self.internal_coeff_deriv
        )

        _, cond_gp = _condition_tinygp_matern32(theta_dict)
        gp_model = cond_gp.mean
        residuals = (self._dataset_res - gp_model)[self._rv_dataset_flag]
        env = 1.0 / self._dataset_ej2[self._rv_dataset_flag]
        n = np.sum(self._rv_dataset_flag)

        return -0.5 * (n * np.log(2 * np.pi) +
                       np.sum(residuals ** 2 * env - np.log(env)))


    def sample_predict(self, dataset, x0_input=None, return_covariance=False, return_variance=False):

        dataset_index = self._dataset_names[dataset.name_ref]
        if x0_input is None:

            l_nstart, l_nend = self._dataset_nindex[dataset_index]
            X_input = self._X

        else:

            l_nstart, l_nend = len(x0_input)*dataset_index, len(x0_input)*(dataset_index+1)

            temp_input = []
            temp_label = []

            for ii in range(0, self._added_datasets):
                temp_input = np.append(temp_input, x0_input)
                temp_label = np.append(temp_label, np.zeros_like(x0_input, dtype=int) + ii)

            X_input = (temp_input, temp_label.astype(int))

        theta_dict =  dict(
            scale=self.internal_parameter_values['matern32_scale'],
            diag=self._dataset_ej2,
            X=self._X,
            y=self._dataset_res,
            coeff_prime=self.internal_coeff_prime,
            coeff_deriv=self.internal_coeff_deriv,
            x0_predict = X_input
        )

        gp = _build_tinygp_multidimensional_matern32(theta_dict)
        _, cond_gp = gp.condition(theta_dict['y'], theta_dict['x0_predict'])

        #mu = cond_gp.mean
        #std = np.sqrt(cond_gp.variance)
        mu_full = cond_gp.loc # or cond_gp.mean?
        mu = mu_full[l_nstart:l_nend]
        std = np.sqrt(cond_gp.variance)[l_nstart:l_nend]

        if return_variance:
            return mu, std
        else:
            return mu
