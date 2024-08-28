from pyorbit.subroutines.common import np, OrderedSet
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.keywords_definitions import *

try:
    from pyorbit.models.celerite_term import celerite, SHOTerm
except (ModuleNotFoundError,ImportError):
    pass


class Celerite_Rotation(AbstractModel):

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
        amp: The amplitude of the variability.
        period: The primary period of variability.
        Q0: The quality factor (or really the quality factor
            minus one half) for the secondary oscillation.
        deltaQ: The difference between the quality factors of the first
            and the second modes. This parameterization (if ``deltaQ > 0``)
            ensures that the primary mode alway has higher quality.
        mix: The fractional amplitude of the secondary mode compared to the
            primary. This should probably always be ``0 < mix < 1``.
    """

    default_common = 'activity'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        try:
            from pyorbit.models.celerite_term import celerite, SHOTerm
        except (ModuleNotFoundError,ImportError):
            print("ERROR: celerite not installed, this will not work")
            quit()

        self.model_class = 'celerite_rotation'
        self.internal_likelihood = True

        self.list_pams_common = OrderedSet([
            'Prot',  # Rotational period of the star
            'Q0',
            'deltaQ',
            'mix'
        ])

        self.list_pams_dataset = OrderedSet([
            'amp',
        ])

        self.n_pams = 6
        self.gp = {}

    def convert_val2gp(self, input_pams):
        """
        :param input_pams: dictionary with the 'physically meaningful'
                           parameters of the GP kernel
        :return: array with the parameters to be fed to 'celerite'

        WARNING: this subroutine is HIGHLY specific of your choice of the
        kernel! I recommend to create a new Class with different
        transformations if you are planning of using a different kernel
        """

        output_pams = np.zeros(self.n_pams, dtype=np.double)

        """ You must check _celerite_ documentation (and possibily do a lot of
            testing) to know how to convert physical values to the parameter
            vector accepted by celerite.set_parameter_vector() function. Note:
            these values may be different from ones accepted by the kernel
        """
        # S0, Q, w0 = output_pams[:3] for SHOterm 1
        output_pams[1] = 0.5 + input_pams['Q0'] + input_pams['deltaQ']
        output_pams[2] = 4 * np.pi * output_pams[1] \
            / (input_pams['Prot'] * np.sqrt(4 * output_pams[1] ** 2 - 1))
        output_pams[0] = input_pams['amp'] \
            / (output_pams[2] * output_pams[1])

        # Another term at half the period
        output_pams[4] = 0.5 + input_pams['Q0']
        output_pams[5] = 8 * np.pi * output_pams[4] \
            / (input_pams['Prot'] * np.sqrt(4 * output_pams[4] ** 2 - 1))
        output_pams[3] = input_pams['mix'] * input_pams['amp'] \
            / (output_pams[5] * output_pams[4])

        return output_pams

    def initialize_model(self, mc,  **kwargs):

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == 'activity':
                self.use_stellar_rotation_period = getattr(mc.common_models[common_ref], 'use_stellar_rotation_period', False)
                break

        for keyword in keywords_stellar_rotation:
            self.use_stellar_rotation_period = kwargs.get(keyword, self.use_stellar_rotation_period)

        if self.use_stellar_rotation_period:
            self.list_pams_common.update(['rotation_period'])
            self.list_pams_common.discard('Prot')

    def initialize_model_dataset(self, mc, dataset, **kwargs):
        self.define_kernel(dataset)
        return

    def define_kernel(self, dataset):
        input_pams = {
            'Prot': 10.0,
            'Q0': 1.0,
            'deltaQ': 0.5,
            'mix': 0.5,
            'amp': 10.0
        }
        gp_pams = self.convert_val2gp(input_pams)

        kernel = SHOTerm(S0=gp_pams[0], Q=gp_pams[1], w0=gp_pams[2]) \
            + SHOTerm(S0=gp_pams[3], Q=gp_pams[4], w0=gp_pams[5])
        self.gp[dataset.name_ref] = celerite.GP(kernel)

        """ I've decided to add the jitter in quadrature instead of using a constant kernel to allow the use of
        different / selective jitter within the dataset
        """
        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].compute(dataset.x0, env)
        return

    def lnlk_compute(self, parameter_values, dataset):
        """ 2 steps:
           1) theta parameters must be converted in physical units (e.g. from logarithmic to linear spaces)
           2) physical values must be converted to {\tt george} input parameters
        """
        if self.use_stellar_rotation_period:
            parameter_values['Prot'] = parameter_values['rotation_period']

        gp_pams = self.convert_val2gp(parameter_values)

        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(dataset.x0, env)

        return self.gp[dataset.name_ref].log_likelihood(dataset.residuals)

    def sample_predict(self, parameter_values, dataset, x0_input=None, return_covariance=False, return_variance=False):

        if self.use_stellar_rotation_period:
            parameter_values['Prot'] = parameter_values['rotation_period']

        gp_pams = self.convert_val2gp(parameter_values)

        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(dataset.x0, env)

        if x0_input is None:
            return self.gp[dataset.name_ref].predict(dataset.residuals, dataset.x0, return_cov=return_covariance, return_var=return_variance)
        else:
            return self.gp[dataset.name_ref].predict(dataset.residuals, x0_input, return_cov=return_covariance, return_var=return_variance)

    def sample_conditional(self, parameter_values, dataset,  x0_input=None):

        if self.use_stellar_rotation_period:
            parameter_values['Prot'] = parameter_values['rotation_period']

        gp_pams = self.convert_val2gp(parameter_values)

        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(dataset.x0, env)

        if x0_input is None:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, dataset.x0)
        else:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, x0_input)
