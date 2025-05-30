from pyorbit.subroutines.common import np, OrderedSet
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_gaussian_processes import AbstractGaussianProcesses
from pyorbit.keywords_definitions import *

try:
    import celerite2
    from celerite2 import terms
except (ModuleNotFoundError,ImportError):
    pass


class Celerite2_Granulation_Rotation(AbstractModel, AbstractGaussianProcesses):

    r"""A mixture of two SHO terms that can be used to model stellar rotation
    This term has two modes in Fourier space: one at ``period`` and one at
    ``0.5 * period``. This can be a good descriptive model for a wide range of
    stochastic variability in stellar time series from rotation to pulsations.
    from Foreman-Mackey+2017 and exoplanet, but keeping the notation of
    the semi-periodic goerge kernel used in PyORBIT
    differently from the example provided in the paper, here the terms are passed in the linear space already. It will
    the job of the sampler to convert from Logarithmic to Linear space for those parameters that the user has decided
    to explore in logarithmic space

    Args:
        SHOTerm (granulation):
        grn_period:  the undamped period of the oscillator
        grn_sigma:  the standard deviation of the process

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
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            import celerite2
        except (ModuleNotFoundError,ImportError):
            print("ERROR: celerite2 not installed, this will not work")
            quit()

        self.model_class = 'gaussian_process'
        self.internal_likelihood = True

        self.list_pams_common = OrderedSet([
            'grn_period',
            'Prot',
            'rot_Q0',
            'rot_deltaQ',
            'rot_fmix',
        ])

        self.list_pams_dataset = OrderedSet([
            'grn_sigma',
            'rot_sigma'
        ])

        self.n_pams = 7
        self.Q_granulation = 1./np.sqrt(2.)
        self.gp = {}

    def initialize_model(self, mc,  **kwargs):

        self._prepare_rotation_replacement(mc, **kwargs)

    def initialize_model_dataset(self, mc, dataset, **kwargs):
        self.define_kernel(dataset)
        return

    def define_kernel(self, dataset):
        """ Kernel initialized with fake values... don't worry, they'll be overwritten soon"""
        kernel = terms.SHOTerm(sigma=1.0, rho=1.0, Q=self.Q_granulation) \
            + terms.RotationTerm(sigma=1.0, period=10., Q0=1.0, dQ=0.5, f=0.5)

        # Setup the GP
        self.gp[dataset.name_ref] = celerite2.GaussianProcess(kernel, mean=0.0)

        """ I've decided to add the jitter in quadrature instead of using a constant kernel to allow the use of
        different / selective jitter within the dataset
        """
        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0

        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag)
        return

    def lnlk_compute(self, parameter_values, dataset):

        self.update_parameter_values(parameter_values)

        """
        In celerite2 the old function "set_parameter_vector" has been removed
        and the kernel has to be defined every time
        """
        self.gp[dataset.name_ref].mean = 0.
        self.gp[dataset.name_ref].kernel = terms.SHOTerm(sigma=parameter_values['grn_sigma'],
                                                        rho=parameter_values['grn_period'],
                                                        Q=self.Q_granulation) \
            + terms.RotationTerm(sigma=parameter_values['rot_sigma'],
                                 period=parameter_values['Prot'],
                                 Q0=parameter_values['rot_Q0'],
                                 dQ=parameter_values['rot_deltaQ'],
                                 f=parameter_values['rot_fmix'])

        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        return self.gp[dataset.name_ref].log_likelihood(dataset.residuals)

    def sample_predict(self, parameter_values, dataset, x0_input=None, return_covariance=False, return_variance=False):

        self.update_parameter_values(parameter_values)

        self.gp[dataset.name_ref].mean = 0.
        self.gp[dataset.name_ref].kernel = terms.SHOTerm(sigma=parameter_values['grn_sigma'],
                                                         rho=parameter_values['grn_period'],
                                                         Q=self.Q_granulation) \
            + terms.RotationTerm(sigma=parameter_values['rot_sigma'],
                                 period=parameter_values['Prot'],
                                 Q0=parameter_values['rot_Q0'],
                                 dQ=parameter_values['rot_deltaQ'],
                                 f=parameter_values['rot_fmix'])

        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        if x0_input is None:
            return self.gp[dataset.name_ref].predict(dataset.residuals, dataset.x0, return_cov=return_covariance, return_var=return_variance)
        else:
            return self.gp[dataset.name_ref].predict(dataset.residuals, x0_input, return_cov=return_covariance, return_var=return_variance)

    def sample_conditional(self, parameter_values, dataset,  x0_input=None):

        self.update_parameter_values(parameter_values)

        self.gp[dataset.name_ref].mean = 0.
        self.gp[dataset.name_ref].kernel = terms.SHOTerm(sigma=parameter_values['grn_sigma'],
                                                         rho=parameter_values['grn_period'],
                                                         Q=self.Q_granulation) \
            + terms.RotationTerm(sigma=parameter_values['rot_sigma'],
                                 period=parameter_values['Prot'],
                                 Q0=parameter_values['rot_Q0'],
                                 dQ=parameter_values['rot_deltaQ'],
                                 f=parameter_values['rot_fmix'])

        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        if x0_input is None:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, dataset.x0)
        else:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, x0_input)
