from pyorbit.classes.common import *
from pyorbit.models.abstract_common import *
from pyorbit.models.abstract_model import *


""" Straight copy from Celerite exmple"""
class Celerite_SemiPeriodic_Term(celerite.terms.Term):
    """

    """

    # from Foreman-Mackey+2017, but keeping the notation of the semi-periodic goerge kernel used in PyORBIT
    # differently from the example provided in the paper, here the terms are passed in the linear space already. It will
    # the job of the sampler to convert from Logarithmic to Linear space for those variables that the user has decided
    # to explore in logarithmic space
    parameter_names = ("Hamp", "Pdec", "Prot", "cel_factor")
    #parameter_names = ("Hamp", "cel_b", "cel_c", "Prot")

    def get_real_coefficients(self, params):
        Hamp, Pdec, Prot, cel_factor = params
        return (
            Hamp * (1.0 + cel_factor) / (2.0 + cel_factor), 1./Pdec,
        )

    def get_complex_coefficients(self, params):
        Hamp, Pdec, Prot, cel_factor = params
        return (
            Hamp / (2.0 + cel_factor),
            0.0,
            1. / Pdec,
            2 * np.pi * (1./Prot),
        )


class Celerite_QuasiPeriodicActivity(AbstractModel):

    internal_likelihood = True

    model_class = 'celerite_quasiperiodic'

    list_pams_common = {
        'Prot', # Rotational period of the star
        'Pdec',
    }

    list_pams_dataset = {
        'Hamp',
        'cel_factor'
    }

    recenter_pams_dataset = {}

    n_pams = 4

    """ Indexing is determined by the way the kernel is constructed, so it is specific of the Model and not of the 
    Common class"""
    gp_pams_index = {
        'Hamp': 0,
        'Pdec': 1,
        'Prot': 2,
        'cel_factor': 3
    }

    def __init__(self, *args, **kwargs):
        super(Celerite_QuasiPeriodicActivity, self).__init__(*args, **kwargs)
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

        """ You must check _george_ documentation (and possibily do a lot of testing) to know how to convert physical 
        values to the parameter vector accepted by celerite.set_parameter_vector() function. Note: these values may be 
        different from ones accepted by the kernel
        """
        output_pams[self.gp_pams_index['Hamp']] = input_pams['Hamp']
        output_pams[self.gp_pams_index['Pdec']] = input_pams['Pdec']
        output_pams[self.gp_pams_index['Prot']] = input_pams['Prot']
        output_pams[self.gp_pams_index['cel_factor']] = input_pams['cel_factor']

        return output_pams

    def convert_gp2val(self, input_pams):
        """
        :param input_pam: array with the parameters to be fed to 'george'
        :return: dictonary with the 'physically meaningful' parameters of the GP kernel
        WARNING: this subroutine is HIGHLY specific of your choice of the kernel! I reccomend to
        create a new Class with different transformations if you are planning of using a different
        kernel combination
        """
        return {
            'Hamp': input_pams[self.gp_pams_index['Hamp']],
            'Pdec': input_pams[self.gp_pams_index['Pdec']],
            'Prot': input_pams[self.gp_pams_index['Prot']],
            'cel_factor': input_pams[self.gp_pams_index['cel_factor']]
        }

    def setup_dataset(self, dataset, **kwargs):
        self.define_kernel(dataset)
        return

    def define_kernel(self, dataset):
        gp_pams = np.ones(self.n_pams)
        kernel = Celerite_SemiPeriodic_Term(Hamp=gp_pams[0], Pdec=gp_pams[1], Prot=gp_pams[2], cel_factor=gp_pams[3])

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
