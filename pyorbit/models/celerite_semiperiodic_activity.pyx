from ..classes.common import *
from abstract_common import *
from abstract_model import *

""" Straight copy from Celerite exmple"""
class Celerite_SemiPeriodic_Term(celerite.terms.Term):
    parameter_names = ("cel_a", "cel_b", "cel_c", "Prot")

    def get_real_coefficients(self, params):
        a, b, c, P = params
        return (
            a * (1.0 + b) / (2.0 + b), c,
        )

    def get_complex_coefficients(self, params):
        a, b, c, P = params
        return (
            a / (2.0 + b), 0.0,
            c, 2 * np.pi * (1./P),
        )


class Celerite_QuasiPeriodicActivity(AbstractModel):
    ''' Three parameters out of four are the same for all the datasets, since they are related to
    the properties of the physical process rather than the observed effects on a dataset
     From Grunblatt+2015, Affer+2016
     - theta: is usually related to the rotation period of the star( or one of its harmonics);
     - lambda: is the correlation decay timescale, and it can be related to the lifetime of the active regions.
     - omega: is the length scale of the periodic component, and can be linked to the size evolution of the active regions;
     - h: represents the amplitude of the correlations '''

    internal_likelihood = True

    model_class = 'celerite_quasiperiodic'

    list_pams_common = {
        'Prot': 'U', # Rotational period of the star
    }

    list_pams_dataset = {
        'cel_a': 'LU',  # celerite term A
        'cel_b': 'LU',  # celerite term B
        'cel_c': 'LU',  # celerite term C
    }

    recenter_pams_dataset = {}

    n_pams = 4

    """ Indexing is determined by the way the kernel is constructed, so it is specific of the Model and not of the 
    Common class"""
    gp_pams_index = {
        'cel_a': 0, # celerite term A
        'cel_b': 1, # celerite term B
        'cel_c': 2, # celerite term C
        'Prot': 3  # ln_P
    }

    gp = {}

    def convert_val2gp(self, input_pams):
        """
        :param input_pam: dictonary with the 'physically meaningful' parameters of the GP kernel
        :return: array with the parameters to be fed to 'celerite'
        WARNING: this subroutine is HIGHLY specific of your choice of the kernel! I reccomend to
        create a new Class with different transformations if you are planning of using a different
        kernel combination
        """
        output_pams = np.zeros(self.n_pams, dtype=np.double)

        """ You must check _george_ documentation (and possibily do a lot of testing) to know how to convert physical 
        values to the parameter vector accepted by george.set_parameter_vector() function. Note: these values may be 
        different from ones accepted by the kernel
        """
        output_pams[self.gp_pams_index['cel_a']] = input_pams['cel_a']
        output_pams[self.gp_pams_index['cel_b']] = input_pams['cel_b']
        output_pams[self.gp_pams_index['cel_c']] = input_pams['cel_c']
        output_pams[self.gp_pams_index['Prot']] = input_pams['Prot']

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
            'cel_a': input_pams[self.gp_pams_index['cel_a']],
            'cel_b': input_pams[self.gp_pams_index['cel_b']],
            'cel_c': input_pams[self.gp_pams_index['cel_c']],
            'Prot': input_pams[self.gp_pams_index['Prot']]
        }

    def setup_dataset(self, dataset, **kwargs):
        self.define_kernel(dataset)
        return

    def define_kernel(self, dataset):
        gp_pams = np.ones(self.n_pams)
        kernel = Celerite_SemiPeriodic_Term(cel_a=gp_pams[0], cel_b=gp_pams[1], cel_c=gp_pams[2], Prot=gp_pams[3])

        """
         gp_pams[0] = h^2 -> h^2 * ExpSquaredKernel * ExpSine2Kernel
           -> set_parameter_vector() accepts the natural logarithm of this value
         gp_pams[1] = metric = r^2 = lambda**2  -> ExpSquaredKernel(metric=r^2)
           -> set_parameter_vector() accepts the natural logarithm of this value
         gp_pams[2] = Gamma =  1/ (2 omega**2) -> ExpSine2Kernel(gamma, ln_period)
         gp_pams[3] = ln_theta = ln_Period -> ExpSine2Kernel(gamma, ln_period)
         
        """

        self.gp[dataset.name_ref] = celerite.GP(kernel)

        """ I've decided to add the jitter in quadrature instead of using a constant kernel to allow the use of 
        different / selective jitter within the dataset
        """
        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].compute(dataset.x0, env)
        return

    def lnlk_compute(self, variable_value, dataset):
        """ 2 steps:
           1) theta parameters must be converted in physical units (e.g. from logarithmic to linear space)
           2) physical values must be converted to {\tt george} input parameters
        """
        gp_pams = self.convert_val2gp(variable_value)
        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(dataset.x0, env)
        #self.gp[dataset.name_ref].recompute()
        return self.gp[dataset.name_ref].log_likelihood(dataset.y - dataset.model)

    def sample_predict(self, variable_value, dataset, x0_input=None):

        gp_pams = self.convert_val2gp(variable_value)

        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(dataset.x0, env)

        if x0_input is None:
            return self.gp[dataset.name_ref].predict(dataset.y - dataset.model, dataset.x0)
        else:
            return self.gp[dataset.name_ref].predict(dataset.y - dataset.model, x0_input)

    def sample_conditional(self, variable_value, dataset,  x0_input=None):

        gp_pams = self.convert_val2gp(variable_value)

        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(dataset.x0, env)

        if x0_input is None:
            return self.gp[dataset.name_ref].sample_conditional(dataset.y - dataset.model, dataset.x0)
        else:
            return self.gp[dataset.name_ref].sample_conditional(dataset.y - dataset.model, x0_input)
