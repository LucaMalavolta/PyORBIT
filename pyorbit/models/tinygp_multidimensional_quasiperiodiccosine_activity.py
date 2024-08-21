from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
from pyorbit.keywords_definitions import *

from scipy.linalg import cho_factor, cho_solve, lapack, LinAlgError
from scipy import matrix, spatial
import sys

__all__ = ['TinyGP_Multidimensional_QuasiPeriodicCosineActivity']


try:
    import jax
    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp
    from tinygp import kernels, GaussianProcess
    #from tinygp.helpers import JAXArray

    if sys.version_info[1] < 10:
        raise Warning("You should be using Python 3.10 - tinygp may not work")

    class LatentKernel_Multi_QuasiPeriodicCosine(kernels.Kernel):
        """A custom kernel based on

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
            kernel_QP : kernels.Kernel
            kernel_CO : kernels.Kernel
            coeff_QP_prim: jax.Array | float
            coeff_QP_deriv: jax.Array | float
            coeff_CO_prim: jax.Array | float
            coeff_CO_deriv: jax.Array | float
        except:
            pass

        def __init__(self, kernel_QP, kernel_CO, coeff_QP_prim, coeff_QP_deriv, coeff_CO_prim, coeff_CO_deriv):
            self.kernel_QP = kernel_QP
            self.kernel_CO = kernel_CO
            self.coeff_QP_prim, self.coeff_QP_deriv = jnp.broadcast_arrays(
                jnp.asarray(coeff_QP_prim), jnp.asarray(coeff_QP_deriv)
            )
            self.coeff_CO_prim, self.coeff_CO_deriv = jnp.broadcast_arrays(
                jnp.asarray(coeff_CO_prim), jnp.asarray(coeff_CO_deriv)
            )
        def evaluate(self, X1, X2):
            t1, label1 = X1
            t2, label2 = X2

            # Differentiate the kernel function: the first derivative wrt x1
            QP_Kp = jax.grad(self.kernel_QP.evaluate, argnums=0)
            CO_Kp = jax.grad(self.kernel_CO.evaluate, argnums=0)

            # ... and the second derivative
            QP_Kpp = jax.grad(QP_Kp, argnums=1)
            CO_Kpp = jax.grad(CO_Kp, argnums=1)

            # Evaluate the kernel matrix and all of its relevant derivatives
            QP_K = self.kernel_QP.evaluate(t1, t2)
            QP_d2K_dx1dx2 = QP_Kpp(t1, t2)

            CO_K = self.kernel_CO.evaluate(t1, t2)
            CO_d2K_dx1dx2 = CO_Kpp(t1, t2)

            # For stationary kernels, these are related just by a minus sign, but we'll
            # evaluate them both separately for generality's sake
            QP_dK_dx2 = jax.grad(self.kernel_QP.evaluate, argnums=1)(t1, t2)
            QP_dK_dx1 = QP_Kp(t1, t2)

            CO_dK_dx2 = jax.grad(self.kernel_CO.evaluate, argnums=1)(t1, t2)
            CO_dK_dx1 = CO_Kp(t1, t2)

            # Extract the coefficients
            a1 = self.coeff_QP_prim[label1]
            a2 = self.coeff_QP_prim[label2]
            b1 = self.coeff_QP_deriv[label1]
            b2 = self.coeff_QP_deriv[label2]

            c1 = self.coeff_CO_prim[label1]
            c2 = self.coeff_CO_prim[label2]
            d1 = self.coeff_CO_deriv[label1]
            d2 = self.coeff_CO_deriv[label2]

            # Construct the matrix element
            return (
                a1 * a2 * QP_K
                + a1 * b2 * QP_dK_dx2
                + b1 * a2 * QP_dK_dx1
                + b1 * b2 * QP_d2K_dx1dx2
                + c1 * c2 * CO_K
                + c1 * d2 * CO_dK_dx2
                + d1 * c2 * CO_dK_dx1
                + d1 * d2 * CO_d2K_dx1dx2
            )


    def _build_tinygp_multidimensional_QPCos(params):

        base_kernel_QP = kernels.ExpSquared(scale=jnp.abs(params["Pdec"])) \
                * kernels.ExpSineSquared(
                scale=jnp.abs(params["Prot"]),
                gamma=jnp.abs(params["gamma"]))

        base_kernel_CO = kernels.ExpSquared(scale=jnp.abs(params["Pdec"])) *  kernels.Cosine(scale=jnp.abs(params["Prot"]/2.))

        kernel = LatentKernel_Multi_QuasiPeriodicCosine(base_kernel_QP, base_kernel_CO,
                                        params['coeff_QP_prime'], params['coeff_QP_deriv'],
                                        params['coeff_CO_prime'], params['coeff_CO_deriv'])
        return GaussianProcess(
            kernel, params['X'], diag=jnp.abs(params['diag']), mean=0.0
        )

    @jax.jit
    def _loss_tinygp_MultiQPCos(params):
        gp = _build_tinygp_multidimensional_QPCos(params)
        return gp.log_probability(params['y'])

except:
    pass




class TinyGP_Multidimensional_QuasiPeriodicCosineActivity(AbstractModel):
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

        self.model_class = 'gp_multidimensional_quasiperiodicsquaredexponential_activity'

        self.internal_likelihood = True
        self.delayed_lnlk_computation = True

        self.list_pams_common = OrderedSet([
            'Prot',  # Rotational period of the star
            'Pdec',  # Decay timescale of activity
            'Oamp',  # Granulation of activity
        ])


        self.list_pams_dataset = OrderedSet([
            'rot_amp', # Amplitude of the covariance matrix
            'con_amp', # Amplitude of the first derivative of the covariance matrix
            'cos_amp', # Amplitude of the covariance matrix
            'cos_der' # Amplitude of the first derivative of the covariance matrix
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

        self.internal_coeff_QP_prime = []
        self.internal_coeff_QP_deriv = []
        self.internal_coeff_CO_prime = []
        self.internal_coeff_CO_deriv = []

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

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        """ when reloading the .p files, the object is not reinitialized, so we have to skip the
        incremental addition of datasets if they are already present  """
        if dataset.name_ref in self._dataset_names:
            return

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

        self.internal_coeff_QP_prime = np.empty(self._added_datasets)
        self.internal_coeff_QP_deriv = np.empty(self._added_datasets)
        self.internal_coeff_CO_prime = np.empty(self._added_datasets)
        self.internal_coeff_CO_deriv = np.empty(self._added_datasets)
        self._X = (self._dataset_x0, self._dataset_label.astype(int))

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

        if 'derivative_quasiperiodic'in kwargs:
            use_derivative_QP = kwargs['derivative_quasiperiodic'].get(dataset.name_ref, False)
        elif dataset.name_ref in kwargs:
            use_derivative_QP = kwargs[dataset.name_ref].get('derivative_quasiperiodic', False)
        else:
            use_derivative_QP = True

        if 'derivative_squaredexponential'in kwargs:
            use_derivative_SE = kwargs['derivative_squaredexponential'].get(dataset.name_ref, False)
        elif dataset.name_ref in kwargs:
            use_derivative_SE= kwargs[dataset.name_ref].get('derivative_squaredexponential', False)
        else:
            use_derivative_SE = True

        if not use_derivative or not use_derivative_QP:
            self.fix_list[dataset.name_ref] = {'rot_amp': [0., 0.]}
        if not use_derivative or not use_derivative_SE:
            self.fix_list[dataset.name_ref] = {'cos_der': [0., 0.]}


        return

    def add_internal_dataset(self, parameter_values, dataset):

        if self.use_stellar_rotation_period:
            parameter_values['Prot'] = parameter_values['rotation_period']

        if self.use_stellar_activity_decay:
            parameter_values['Pdec'] = parameter_values['activity_decay']

        self.internal_parameter_values = parameter_values

        d_ind = self._dataset_names[dataset.name_ref]
        d_nstart, d_nend = self._dataset_nindex[d_ind]

        self._dataset_ej2[d_nstart:d_nend] = self._dataset_e2[d_nstart:d_nend] + dataset.jitter**2.0
        self._dataset_res[d_nstart:d_nend] = dataset.residuals

        self.internal_coeff_QP_prime[d_ind] = parameter_values['con_amp']
        self.internal_coeff_QP_deriv[d_ind] = parameter_values['rot_amp']
        self.internal_coeff_CO_prime[d_ind] = parameter_values['cos_amp']
        self.internal_coeff_CO_deriv[d_ind] = parameter_values['cos_der']

    def lnlk_compute(self):
        if not self.hyper_condition(self.internal_parameter_values):
            return -np.inf
        if not self.rotdec_condition(self.internal_parameter_values):
            return -np.inf

        theta_dict =  dict(
            gamma=1. / (2.*self.internal_parameter_values['Oamp'] ** 2),
            Pdec=self.internal_parameter_values['Pdec'],
            Prot=self.internal_parameter_values['Prot'],
            diag=self._dataset_ej2,
            X=self._X,
            y=self._dataset_res,
            coeff_QP_prime=self.internal_coeff_QP_prime,
            coeff_QP_deriv=self.internal_coeff_QP_deriv,
            coeff_CO_prime=self.internal_coeff_CO_prime,
            coeff_CO_deriv=self.internal_coeff_CO_deriv
        )

        return _loss_tinygp_MultiQPCos(theta_dict)


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
            gamma=1. / (2.*self.internal_parameter_values['Oamp'] ** 2),
            Pdec=self.internal_parameter_values['Pdec'],
            Prot=self.internal_parameter_values['Prot'],
            diag=self._dataset_ej2,
            X=self._X,
            y=self._dataset_res,
            coeff_QP_prime=self.internal_coeff_QP_prime,
            coeff_QP_deriv=self.internal_coeff_QP_deriv,
            coeff_CO_prime=self.internal_coeff_CO_prime,
            coeff_CO_deriv=self.internal_coeff_CO_deriv,
            x0_predict = X_input
        )

        gp = _build_tinygp_multidimensional_QPCos(theta_dict)
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
