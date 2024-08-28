from pyorbit.subroutines.common import np, OrderedSet
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.keywords_definitions import *

try:
    import celerite2
    from celerite2 import terms
except (ModuleNotFoundError,ImportError):
    pass


class Celerite2_Granulation_Oscillation_Rotation(AbstractModel):

    r"""A


    Args:
        SHOTerm (granulation):
        grn_period:  the undamped period of the oscillator
        grn_sigma:  the standard deviation of the process
        Q_factor = 1 / np.sqrt(2)

        SHOTerm (oscillations):
        osc_period:  the undamped period of the oscillator
        osc_sigma:  the standard deviation of the process
        osc_qfactor > 1.

        RotationTerm (rotation):
        rot_sigma: The standard deviation of the process.
        Prot - rot_period: The primary period of variability,
            named as the quasi-period kernel parameter to preserve compatiblity
        rot_Q0: The quality factor (or really the quality factor
            minus one half) for the secondary oscillation.
        rot_deltaQ: The difference between the quality factors of the first
            and the second modes. This parameterization (if ``deltaQ > 0``)
            ensures that the primary mode alway has higher quality.
        rot_fmix: The fractional amplitude of the secondary mode compared to the
            primary. This should probably always be ``0 < mix < 1``.
    """

    default_common = 'activity'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        try:
            import celerite2
        except (ModuleNotFoundError,ImportError):
            print("ERROR: celerite2 not installed, this will not work")
            quit()

        self.model_class = 'celerite2_granulation_rotation'
        self.internal_likelihood = True

        self.list_pams_common = OrderedSet()
        self.list_pams_dataset = OrderedSet()

        self.n_pams = 7
        self.Q_granulation = 1./np.sqrt(2.)
        self.gp = {}


        self.rotation_kernels = 1
        self.granulation_kernels = 2
        self.oscillation_kernels = 1

    def initialize_model(self, mc, **kwargs):

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == 'activity':
                self.use_stellar_rotation_period = getattr(mc.common_models[common_ref], 'use_stellar_rotation_period', False)
                break

        for keyword in keywords_stellar_rotation:
            self.use_stellar_rotation_period = kwargs.get(keyword, self.use_stellar_rotation_period)

        self.rotation_kernels = kwargs.get('rotation_kernels', 1)
        self.granulation_kernels = kwargs.get('granulation_kernels', 2)
        self.oscillation_kernels = kwargs.get('oscillation_kernels', 1)

        """ Only one rotation kernel is allowed """
        if self.rotation_kernels == 1 :

            if self.use_stellar_rotation_period:
                self.list_pams_common.update(['rotation_period'])
            else:
                self.list_pams_common.update(['Prot'])

            self.list_pams_common.update(['rot_Q0'])
            self.list_pams_common.update(['rot_deltaQ'])
            self.list_pams_common.update(['rot_fmix'])
            if kwargs.get('common_amplitudes', False):
                self.list_pams_common.update(['rot_sigma'])
            else:
                self.list_pams_dataset.update(['rot_sigma'])

        elif self.rotation_kernels > 1 :
            print(" Celerite2_Granulation_Oscillation_Rotation ERROR: only one rotational kernel supported")
            quit()

        for i_k in range(0, self.granulation_kernels):
            if i_k >= 10:
                print(" Celerite2_Granulation_Oscillation_Rotation WARNING: only up to 10 kernel supported ")
            self.list_pams_common.update(['grn_k'+repr(i_k) + '_period'])
            if kwargs.get('common_amplitudes', False):
                self.list_pams_common.update(['grn_k'+repr(i_k) + '_sigma'])
            else:
                self.list_pams_dataset.update(['grn_k'+repr(i_k) + '_sigma'])

        for i_k in range(0, self.oscillation_kernels):
            if i_k >= 10:
                print(" Celerite2_Granulation_Oscillation_Rotation WARNING: only up to 10 kernel supported ")
            self.list_pams_common.update(['osc_k'+repr(i_k) + '_period'])
            self.list_pams_common.update(['osc_k'+repr(i_k) + '_Q0'])
            if kwargs.get('common_amplitudes', False):
                self.list_pams_common.update(['osc_k'+repr(i_k) + '_sigma'])
            else:
                self.list_pams_dataset.update(['osc_k'+repr(i_k) + '_sigma'])

    def initialize_model_dataset(self, mc, dataset, **kwargs):
        self.define_kernel(dataset)
        return

    def define_kernel(self, dataset):
        """ Kernel initialized with fake values... don't worry, they'll be overwritten soon"""

        i_kernels = 0

        if self.rotation_kernels == 1:
            kernel = terms.RotationTerm(sigma=1.0, period=10., Q0=1.0, dQ=0.5, f=0.5)
            i_kernels += 1

        for i_k in range(0, self.granulation_kernels):
            if i_kernels > 0:
                kernel += terms.SHOTerm(sigma=1.0, rho=1.0, Q=self.Q_granulation)
            else:
                kernel = terms.SHOTerm(sigma=1.0, rho=1.0, Q=self.Q_granulation)
            i_kernels += 1

        for i_k in range(0, self.oscillation_kernels):
            if i_kernels > 0:
                kernel += terms.SHOTerm(sigma=1.0, rho=1.0, Q=2.0)
            else:
                kernel = terms.SHOTerm(sigma=1.0, rho=1.0, Q=2.0)
            i_kernels += 1

        # Setup the GP
        self.gp[dataset.name_ref] = celerite2.GaussianProcess(kernel, mean=0.0)

        """ I've decided to add the jitter in quadrature instead of using a constant kernel to allow the use of
        different / selective jitter within the dataset
        """
        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0

        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag)
        return

    def lnlk_compute(self, parameter_values, dataset):
        """
        In celerite2 the old function "set_parameter_vector" has been removed
        and the kernel has to be defined every time
        """
        if self.use_stellar_rotation_period:
            parameter_values['Prot'] = parameter_values['rotation_period']

        self.gp[dataset.name_ref].mean = 0.
        i_kernels = 0
        if self.rotation_kernels == 1:
            kernel = terms.RotationTerm(sigma=parameter_values['rot_sigma'],
                                period=parameter_values['Prot'],
                                Q0=parameter_values['rot_Q0'],
                                dQ=parameter_values['rot_deltaQ'],
                                f=parameter_values['rot_fmix'])
            i_kernels += 1

        for i_k in range(0, self.granulation_kernels):
            if i_kernels > 0:
                kernel += terms.SHOTerm(sigma=parameter_values['grn_k'+repr(i_k) + '_sigma'],
                                                         rho=parameter_values['grn_k'+repr(i_k) + '_period'],
                                                         Q=self.Q_granulation)
            else:
                kernel = terms.SHOTerm(sigma=parameter_values['grn_k'+repr(i_k) + '_sigma'],
                                                         rho=parameter_values['grn_k'+repr(i_k) + '_period'],
                                                         Q=self.Q_granulation)
            i_kernels += 1

        for i_k in range(0, self.oscillation_kernels):
            if i_kernels > 0:
                kernel += terms.SHOTerm(sigma=parameter_values['osc_k'+repr(i_k) + '_sigma'],
                                                         rho=parameter_values['osc_k'+repr(i_k) + '_period'],
                                                         Q=parameter_values['osc_k'+repr(i_k) + '_Q0'])
            else:
                kernel = terms.SHOTerm(sigma=parameter_values['osc_k'+repr(i_k) + '_sigma'],
                                                         rho=parameter_values['osc_k'+repr(i_k) + '_period'],
                                                         Q=parameter_values['osc_k'+repr(i_k) + '_Q0'])
            i_kernels += 1


        self.gp[dataset.name_ref].kernel = kernel
        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        return self.gp[dataset.name_ref].log_likelihood(dataset.residuals)

    def sample_predict(self, parameter_values, dataset, x0_input=None, return_covariance=False, return_variance=False):

        if self.use_stellar_rotation_period:
            parameter_values['Prot'] = parameter_values['rotation_period']

        self.gp[dataset.name_ref].mean = 0.

        i_kernels = 0
        if self.rotation_kernels == 1:
            kernel = terms.RotationTerm(sigma=parameter_values['rot_sigma'],
                                period=parameter_values['Prot'],
                                Q0=parameter_values['rot_Q0'],
                                dQ=parameter_values['rot_deltaQ'],
                                f=parameter_values['rot_fmix'])
            i_kernels += 1

        for i_k in range(0, self.granulation_kernels):
            if i_kernels > 0:
                kernel += terms.SHOTerm(sigma=parameter_values['grn_k'+repr(i_k) + '_sigma'],
                                                         rho=parameter_values['grn_k'+repr(i_k) + '_period'],
                                                         Q=self.Q_granulation)
            else:
                kernel = terms.SHOTerm(sigma=parameter_values['grn_k'+repr(i_k) + '_sigma'],
                                                         rho=parameter_values['grn_k'+repr(i_k) + '_period'],
                                                         Q=self.Q_granulation)
            i_kernels += 1

        for i_k in range(0, self.oscillation_kernels):
            if i_kernels > 0:
                kernel += terms.SHOTerm(sigma=parameter_values['osc_k'+repr(i_k) + '_sigma'],
                                                         rho=parameter_values['osc_k'+repr(i_k) + '_period'],
                                                         Q=parameter_values['osc_k'+repr(i_k) + '_Q0'])
            else:
                kernel = terms.SHOTerm(sigma=parameter_values['osc_k'+repr(i_k) + '_sigma'],
                                                         rho=parameter_values['osc_k'+repr(i_k) + '_period'],
                                                         Q=parameter_values['osc_k'+repr(i_k) + '_Q0'])
            i_kernels += 1


        self.gp[dataset.name_ref].kernel = kernel
        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        if x0_input is None:
            return self.gp[dataset.name_ref].predict(dataset.residuals, dataset.x0, return_cov=return_covariance, return_var=return_variance)
        else:
            return self.gp[dataset.name_ref].predict(dataset.residuals, x0_input, return_cov=return_covariance, return_var=return_variance)

    def sample_conditional(self, parameter_values, dataset,  x0_input=None):

        if self.use_stellar_rotation_period:
            parameter_values['Prot'] = parameter_values['rotation_period']

        self.gp[dataset.name_ref].mean = 0.
        i_kernels = 0
        if self.rotation_kernels == 1:
            kernel = terms.RotationTerm(sigma=parameter_values['rot_sigma'],
                                period=parameter_values['Prot'],
                                Q0=parameter_values['rot_Q0'],
                                dQ=parameter_values['rot_deltaQ'],
                                f=parameter_values['rot_fmix'])
            i_kernels += 1

        for i_k in range(0, self.granulation_kernels):
            if i_kernels > 0:
                kernel += terms.SHOTerm(sigma=parameter_values['grn_k'+repr(i_k) + '_sigma'],
                                                         rho=parameter_values['grn_k'+repr(i_k) + '_period'],
                                                         Q=self.Q_granulation)
            else:
                kernel = terms.SHOTerm(sigma=parameter_values['grn_k'+repr(i_k) + '_sigma'],
                                                         rho=parameter_values['grn_k'+repr(i_k) + '_period'],
                                                         Q=self.Q_granulation)
            i_kernels += 1

        for i_k in range(0, self.oscillation_kernels):
            if i_kernels > 0:
                kernel += terms.SHOTerm(sigma=parameter_values['osc_k'+repr(i_k) + '_sigma'],
                                                         rho=parameter_values['osc_k'+repr(i_k) + '_period'],
                                                         Q=parameter_values['osc_k'+repr(i_k) + '_Q0'])
            else:
                kernel = terms.SHOTerm(sigma=parameter_values['osc_k'+repr(i_k) + '_sigma'],
                                                         rho=parameter_values['osc_k'+repr(i_k) + '_period'],
                                                         Q=parameter_values['osc_k'+repr(i_k) + '_Q0'])
            i_kernels += 1


        self.gp[dataset.name_ref].kernel = kernel
        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        if x0_input is None:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, dataset.x0)
        else:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, x0_input)
