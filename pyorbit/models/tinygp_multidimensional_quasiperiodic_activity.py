from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *

from scipy.linalg import cho_factor, cho_solve, lapack, LinAlgError
from scipy import matrix, spatial

try:
    import jax
    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp
    from tinygp import kernels, GaussianProcess

except:
    pass


__all__ = ['TinyGP_Multidimensional_QuasiPeriodicActivity']


class LatentKernel(kernels.Kernel):
    """A custom kernel based on Rajpaul et al. (2015)

    Args:
        kernel: The kernel function describing the latent process. This can be any other
            ``tinygp`` kernel.
        coeff_prim: The primal coefficients for each class. This can be thought of as how
            much the latent process itself projects into the observations for that class.
            This should be an array with an entry for each class of observation.
        coeff_deriv: The derivative coefficients for each class. This should have the same
            shape as ``coeff_prim``.
    """

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



def _build_tinygp_multidimensional(params):

    base_kernel = kernels.ExpSquared(scale=jnp.abs(params["Pdec"])) \
            * kernels.ExpSineSquared(
            scale=jnp.abs(params["Prot"]),
            gamma=jnp.abs(params["gamma"]))

    kernel = LatentKernel(base_kernel, params['coeff_prime'], params['coeff_deriv'])
    return GaussianProcess(
        kernel, params['X'], diag=jnp.abs(params['diag']), mean=0.0
    )

@jax.jit
def _loss_tinygp(params):
    gp = _build_tinygp_multidimensional(params)
    return gp.log_probability(params['y'])

@jax.jit
def _residuals_tinygp(params):
    gp = _build_tinygp_multidimensional(params)
    _, cond_gp = gp.condition(params['y'], params['X'])
    return cond_gp


class TinyGP_Multidimensional_QuasiPeriodicActivity(AbstractModel):
    ''' Three parameters out of four are the same for all the datasets, since they are related to
    the properties of the physical process rather than the observed effects on a dataset
     From Grunblatt+2015, Affer+2016
     - theta: is usually related to the rotation period of the star( or one of its harmonics);
     - lambda: is the correlation decay timescale, and it can be related to the lifetime of the active regions.
     - omega: is the length scale of the periodic component, and can be linked to the size evolution of the active regions;
     - h: represents the amplitude of the correlations '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'gp_multidimensional_quasiperiodic_activity'

        self.internal_likelihood = True
        self.delayed_lnlk_computation = True

        self.list_pams_common = {
            'Prot',  # Rotational period of the star
            'Pdec',  # Decay timescale of activity
            'Oamp',  # Granulation of activity
        }
        self.list_pams_dataset = {
            'rot_amp', # Amplitude of the covariance matrix
            'con_amp' # Amplitude of the first derivative of the covariance matrix
        }


        self.internal_variable_value = None
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


    def initialize_model(self, mc,  **kwargs):

        if kwargs.get('hyperparameters_condition', False):
            self.hyper_condition = self._hypercond_01
        else:
            self.hyper_condition = self._hypercond_00

        if kwargs.get('rotation_decay_condition', False):
            self.rotdec_condition = self._hypercond_02
        else:
            self.rotdec_condition = self._hypercond_00

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        """ when reloading the .p files, the object is not reinitialized, so we have to skip the
        incremental addition of datasets if they are already present  """
        if dataset.name_ref in self._dataset_names:
            return

        self._dataset_nindex.append([self._n_cov_matrix,
                                    self._n_cov_matrix+dataset.n])

        self._dataset_x0 = np.append(self._dataset_x0, dataset.x0)
        self._dataset_label = np.append(self._dataset_label, np.zeros_like(dataset.x0)+ self._added_datasets)
        self._dataset_e2 = np.append(self._dataset_e2, dataset.e**2)

        self._dataset_names[dataset.name_ref] = self._added_datasets
        self._n_cov_matrix += dataset.n
        self._added_datasets += 1

        self._dataset_ej2 = self._dataset_e2 * 1.
        self._dataset_res = self._dataset_e2 * 0.

        self.internal_coeff_prime = np.empty(self._added_datasets)
        self.internal_coeff_deriv = np.empty(self._added_datasets)
        self._X = (self._dataset_x0, self._dataset_label)

        if 'derivative'in kwargs:
            use_derivative = kwargs['derivative'].get(dataset.name_ref, False)
        elif dataset.name_ref in kwargs:
            use_derivative = kwargs[dataset.name_ref].get('derivative', False)
        else:
            if dataset.kind == 'H-alpha' or \
                dataset.kind == 'S_index' or \
                dataset.kind == 'Ca_HK' or \
                dataset.kind == 'FWHM':
                    use_derivative = False
            else:
                use_derivative = True

        if not use_derivative:
            self.fix_list[dataset.name_ref] = {'rot_amp': [0., 0.]}


        return

    def add_internal_dataset(self, variable_value, dataset):

        self.internal_variable_value = variable_value

        d_ind = self._dataset_names[dataset.name_ref]
        d_nstart, d_nend = self._dataset_nindex[d_ind]

        self._dataset_ej2[d_nstart:d_nend] = self._dataset_e2[d_nstart:d_nend] + dataset.jitter**2.0
        self._dataset_res[d_nstart:d_nend] = dataset.residuals

        self.internal_coeff_prime[d_ind] = variable_value['con_amp']
        self.internal_coeff_deriv[d_ind] = variable_value['rot_amp']

    def lnlk_compute(self):
        if not self.hyper_condition(self.internal_variable_value):
            return -np.inf
        if not self.rotdec_condition(self.internal_variable_value):
            return -np.inf

        theta_dict =  dict(
            gamma=1. / (2.*self.internal_variable_value['Oamp'] ** 2),
            Pdec=self.internal_variable_value['Pdec'],
            Prot=self.internal_variable_value['Prot'],
            diag=self._dataset_ej2,
            X=self._X,
            y=self._dataset_res,
            coeff_prime=self.internal_coeff_prime,
            coeff_deriv=self.internal_coeff_deriv
        )

        return _loss_tinygp(theta_dict)


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
                temp_label = np.append(temp_label, np.zeros_like(x0_input, dtype=int)+ii)
            X_input = (temp_input, temp_label)

        theta_dict =  dict(
            gamma=1. / (2.*self.internal_variable_value['Oamp'] ** 2),
            Pdec=self.internal_variable_value['Pdec'],
            Prot=self.internal_variable_value['Prot'],
            diag=self._dataset_ej2,
            X=X_input,
            y=self._dataset_res,
            coeff_prime=self.internal_coeff_prime,
            coeff_deriv=self.internal_coeff_deriv
        )

        _, cond_gp = _residuals_tinygp(theta_dict)
        mu = cond_gp.mean[l_nstart:l_nend]
        std = np.sqrt(cond_gp.variance[l_nstart:l_nend])
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
