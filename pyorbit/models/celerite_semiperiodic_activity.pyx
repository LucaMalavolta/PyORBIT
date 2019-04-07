from pyorbit.classes.common import *
from pyorbit.models.abstract_common import *
from pyorbit.models.abstract_model import *


class Celerite_QuasiPeriodicActivity(AbstractModel):

    internal_likelihood = True

    model_class = 'celerite_quasiperiodic'

    list_pams_common = {
        'Prot', # Rotational period of the star
        'Pdec',
    }

    list_pams_dataset = {
        'cel_B',
        'cel_C'
    }

    recenter_pams_dataset = {}

    n_pams = 5

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

        """ You must check _celerite_ documentation (and possibily do a lot of testing) to know how to convert physical 
        values to the parameter vector accepted by celerite.set_parameter_vector() function. Note: these values may be 
        different from ones accepted by the kernel
        """
        B = input_pams['cel_B']
        Pdec = input_pams['Pdec']
        Prot = input_pams['Prot']
        C = input_pams['cel_C']

        output_pams[0] = np.log(B*(1+C)/(2+C))
        output_pams[1] = np.log(1/Pdec)
        output_pams[2] = np.log(B/(2+C))
        output_pams[3] = np.log(1/Pdec)
        output_pams[4] = np.log(2*np.pi/Prot)

        return output_pams

    def setup_dataset(self, dataset, **kwargs):
        self.define_kernel(dataset)
        return

    def define_kernel(self, dataset):
        gp_pams = np.ones(self.n_pams)

        kernel = celerite.terms.RealTerm(
            log_a=gp_pams[0],
            log_c=gp_pams[1]
        )

        kernel += celerite.terms.ComplexTerm(
            log_a=gp_pams[2],
            log_c=gp_pams[3],
            log_d=gp_pams[4]
        )

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
