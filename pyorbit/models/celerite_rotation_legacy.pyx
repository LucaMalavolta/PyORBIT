from pyorbit.models.abstract_model import AbstractModel

# This class was written by Daniel Foreman-Mackey for his paper:
# https://github.com/dfm/celerite/blob/master/paper/figures/rotation/rotation.ipynb


class RotationTerm(Term):
    parameter_names = ("period", "timescale", "factor", "amp")

    def get_real_coefficients(self, params):
        # in DFM formulation: P, c, b, a
        period, timescale, factor, amp = params
        return (
            amp * (1.0 + factor) / (2.0 + factor),
            1./timescale,
        )

    def get_complex_coefficients(self, params):
        period, timescale, factor, amp = params
        return (
            amp / (2.0 + factor),
            0.0,
            1./timescale,
            2 * np.pi * 1./period,
        )


class Celerite_Rotation_Legacy(AbstractModel):

    r"""A mixture of two SHO terms that can be used to model stellar rotation
    This term has two modes in Fourier space: one at ``period`` and one at
    ``0.5 * period``. This can be a good descriptive model for a wide range of
    stochastic variability in stellar time series from rotation to pulsations.
    from Foreman-Mackey+2017 and exoplanet, but keeping the notation of 
    the semi-periodic goerge kernel used in PyORBIT
    differently from the example provided in the paper, here the terms are passed in the linear space already. It will
    the job of the sampler to convert from Logarithmic to Linear space for those variables that the user has decided
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

    model_class = 'celerite_rotation_legacy'
    internal_likelihood = True

    def __init__(self, *args, **kwargs):
        super(Celerite_Rotation_Legacy, self).__init__(*args, **kwargs)

        try:
            import celerite
        except:
            print("ERROR: celerite not installed, this will not work")
            quit()


        self.list_pams_common = {
            'Prot',  # Rotational period of the star
            'Pdec',
            'cel_factor',
        }

        self.list_pams_dataset = {
            'Hamp',
        }

        self.n_pams = 4
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
        output_pams[0] = input_pams['Prot']
        output_pams[1] = input_pams['Pdec']
        output_pams[2] = input_pams['cel_factor']
        output_pams[3] = input_pams['Hamp']

        return output_pams

    def setup_dataset(self, mc, dataset, **kwargs):
        self.define_kernel(dataset)
        return

    def define_kernel(self, dataset):
        input_pams = {
            'Prot': 10.0,
            'Pdec': 1000.,
            'cel_factor': 0.5,
            'Hamp': 10.
        }
        gp_pams = self.convert_val2gp(input_pams)

        kernel = RotationTerm(period=gp_pams[0],
                              timescale=gp_pams[1], 
                              factor=gp_pams[2], 
                              amp=gp_pams[3])
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

    def sample_predict(self, variable_value, dataset, x0_input=None, return_covariance=False, return_variance=False):

        gp_pams = self.convert_val2gp(variable_value)

        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(dataset.x0, env)

        if x0_input is None:
            return self.gp[dataset.name_ref].predict(dataset.residuals, dataset.x0, return_cov=return_covariance, return_var=return_variance)
        else:
            return self.gp[dataset.name_ref].predict(dataset.residuals, x0_input, return_cov=return_covariance, return_var=return_variance)

    def sample_conditional(self, variable_value, dataset,  x0_input=None):

        gp_pams = self.convert_val2gp(variable_value)

        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(dataset.x0, env)

        if x0_input is None:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, dataset.x0)
        else:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, x0_input)
