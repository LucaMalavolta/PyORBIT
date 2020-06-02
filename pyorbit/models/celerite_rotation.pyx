from pyorbit.classes.common import *
from pyorbit.models.abstract_common import *
from pyorbit.models.abstract_model import *
from celerite.terms import SHOTerm, TermSum

class Celerite_Rotation_Term(TermSum):
    """

    """

    """A mixture of two SHO terms that can be used to model stellar rotation
    This term has two modes in Fourier space: one at ``period`` and one at
    ``0.5 * period``. This can be a good descriptive model for a wide range of
    stochastic variability in stellar time series from rotation to pulsations.
    Args:
        amp: The amplitude of the variability.
        tensor perio: The primary period of variability.
        tensor Q0 or log_Q0: The quality factor (or really the quality factor
            minus one half) for the secondary oscillation.
        tensor deltaQ or log_deltaQ: The difference between the quality factors
            of the first and the second modes. This parameterization (if
            ``deltaQ > 0``) ensures that the primary mode alway has higher
            quality.
        mix: The fractional amplitude of the secondary mode compared to the
            primary. This should probably always be ``0 < mix < 1``.
    """

    # from Foreman-Mackey+2017 and exoplanet, but keeping the notation of the semi-periodic goerge kernel used in PyORBIT
    # differently from the example provided in the paper, here the terms are passed in the linear space already. It will
    # the job of the sampler to convert from Logarithmic to Linear space for those variables that the user has decided
    # to explore in logarithmic space

    parameter_names = ("period", "Q0", "deltaQ", "mix", "amp")
    def __init__(self, **kwargs):
        super(Celerite_Rotation_Term, self).__init__(**kwargs)

        # One term with a period of period
        Q1 = 0.5 + self.Q0 + self.deltaQ
        w1 = 4 * np.pi * Q1 / (self.period * np.sqrt(4 * Q1 ** 2 - 1))
        S1 = self.amp / (w1 * Q1)

        # Another term at half the period
        Q2 = 0.5 + self.Q0
        w2 = 8 * np.pi * Q2 / (self.period * np.sqrt(4 * Q2 ** 2 - 1))
        S2 = self.mix * self.amp / (w2 * Q2)

        self.terms = (SHOTerm(S0=S1, w0=w1, Q=Q1), SHOTerm(S0=S2, w0=w2, Q=Q2))
        self.coefficients = self.get_coefficients()


class Celerite_Rotation(AbstractModel):

    internal_likelihood = True

    model_class = 'celerite_quasiperiodic'

    list_pams_common = {
        'Prot',  # Rotational period of the star
        'Q0',
        'deltaQ',
        'mix'
    }

    list_pams_dataset = {
        'amp',
    }

    recenter_pams_dataset = {}

    n_pams = 5

    def __init__(self, *args, **kwargs):
        super(Celerite_Rotation, self).__init__(*args, **kwargs)
        self.gp = {}

    def convert_val2gp(self, input_pams):
        """
        :param input_pams: dictionary with the 'physically meaningful' parameters of the GP kernel
        :return: array with the parameters to be fed to 'celerite'
        WARNING: this subroutine is HIGHLY specific of your choice of the kernel! I recommend to
        create a new Class with different transformations if you are planning of using a different
        kernel combination
        """

        output_pams = np.zeros(self.n_pams, dtype=np.double)

        """ You must check _celerite_ documentation (and possibily do a lot of testing) to know how to convert physical 
        values to the parameter vector accepted by celerite.set_parameter_vector() function. Note: these values may be 
        different from ones accepted by the kernel
        """

        output_pams[0] = input_pams['Prot']
        output_pams[1] = input_pams['Q0']
        output_pams[2] = input_pams['deltaQ']
        output_pams[3] = input_pams['mix']
        output_pams[4] = input_pams['amp']

        return output_pams

    def setup_dataset(self, mc, dataset, **kwargs):
        self.define_kernel(dataset)
        return

    def define_kernel(self, dataset):
        gp_pams = np.ones(self.n_pams)

        kernel = celerite.terms.Celerite_Rotation_Term(gp_pams)
        self.gp[dataset.name_ref] = celerite.GP(kernel)

        """ I've decided to add the jitter in quadrature instead of using a constant kernel to allow the use of 
        different / selective jitter within the dataset
        """
        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].compute(dataset.x0, env)
        return

    def lnlk_compute(self, variable_value, dataset):
        """ 2 steps:
           1) theta parameters must be converted in physical units (e.g. from logarithmic to linear spaces)
           2) physical values must be converted to {\tt george} input parameters
        """
        gp_pams = self.convert_val2gp(variable_value)

        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(dataset.x0, env)

        return self.gp[dataset.name_ref].log_likelihood(dataset.residuals)

    def sample_predict(self, variable_value, dataset, x0_input=None):

        gp_pams = self.convert_val2gp(variable_value)

        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(dataset.x0, env)

        if x0_input is None:
            return self.gp[dataset.name_ref].predict(dataset.residuals, dataset.x0, return_var=True)
        else:
            return self.gp[dataset.name_ref].predict(dataset.residuals, x0_input, return_var=True)

    def sample_conditional(self, variable_value, dataset,  x0_input=None):

        gp_pams = self.convert_val2gp(variable_value)

        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(dataset.x0, env)

        if x0_input is None:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, dataset.x0)
        else:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, x0_input)
