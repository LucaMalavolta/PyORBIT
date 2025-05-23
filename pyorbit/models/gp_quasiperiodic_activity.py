from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_gaussian_processes import AbstractGaussianProcesses
from pyorbit.keywords_definitions import *

try:
    import george
except:
    pass


class GaussianProcess_QuasiPeriodicActivity(AbstractModel, AbstractGaussianProcesses):
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
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            import george
        except (ModuleNotFoundError,ImportError):
            print("ERROR: george not installed, this will not work")
            quit()

        self.model_class = 'gaussian_process'
        self.internal_likelihood = True

        self.list_pams_common = OrderedSet([
            'Prot',  # Rotational period of the star
            'Pdec',  # Decay timescale of activity
            'Oamp'  # Granulation of activity
        ])

        self.list_pams_dataset = OrderedSet([
            'Hamp'  # Amplitude of the signal in the covariance matrix
        ])

        self.n_pams = 4

        self.gp_pams_index = {
            'Hamp': 0,  # amp2
            'Pdec': 1,  # metric
            'Oamp': 2,  # gamma
            'Prot': 3  # ln_P
        }

        self.gp = {}
        self.use_HODLR = False

    def convert_val2gp(self, input_pams):
        """
        :param input_pam: dictonary with the 'physically meaningful' parameters of the GP kernel
        :return: dictonary with the parameters to be fed to 'george'
        WARNING: this subroutine is HIGHLY specific of your choice of the kernel! I reccomend to
        create a new Class with different transformations if you are planning of using a different
        kernel combination
        """
        output_pams = np.zeros(self.n_pams, dtype=np.double)

        """ You must check _george_ documentation (and possibily do a lot of testing) to know how to convert physical
        values to the parameter vector accepted by george.set_parameter_vector() function. Note: these values may be
        different from ones accepted by the kernel
        """
        output_pams[self.gp_pams_index['Hamp']] = np.log(input_pams['Hamp'])*2
        output_pams[self.gp_pams_index['Pdec']] = np.log(input_pams['Pdec'])*2
        output_pams[self.gp_pams_index['Oamp']] = 1. / \
            (2.*input_pams['Oamp'] ** 2)
        output_pams[self.gp_pams_index['Prot']] = np.log(input_pams['Prot'])

        return output_pams

    def convert_gp2val(self, input_pams):
        """
        :param input_pam: dictonary with the parameters to be fed to 'george'
        :return: dictonary with the 'physically meaningful' parameters of the GP kernel
        WARNING: this subroutine is HIGHLY specific of your choice of the kernel! I reccomend to
        create a new Class with different transformations if you are planning of using a different
        kernel combination
        """
        return {
            'Hamp': np.exp(input_pams[self.gp_pams_index['Hamp']]/2.0),
            'Pdec': np.exp(input_pams[self.gp_pams_index['Pdec']]/2.0),
            'Oamp': np.sqrt(1. / (2.*input_pams[self.gp_pams_index['Oamp']])),
            'Prot': np.exp(input_pams[self.gp_pams_index['Prot']])
        }

    def initialize_model(self, mc,  **kwargs):

        if 'use_HODLR' in kwargs:
            self.use_HODLR = kwargs['use_HODLR']

        self._prepare_hyperparameter_conditions(mc, **kwargs)
        self._prepare_rotation_replacement(mc, **kwargs)
        self._prepare_decay_replacement(mc, **kwargs)

    def initialize_model_dataset(self, mc, dataset, **kwargs):
        self.define_kernel(dataset)
        return

    def define_kernel(self, dataset):

        gp_pams = np.ones(self.n_pams)
        """ Kernel initialized with fake values... don't worry, they'll be overwritten soon"""
        kernel = np.exp(gp_pams[0]) * \
            george.kernels.ExpSquaredKernel(metric=np.exp(gp_pams[1])) * \
            george.kernels.ExpSine2Kernel(
                gamma=gp_pams[2], log_period=gp_pams[3])

        """
         gp_pams[0] = h^2 -> h^2 * ExpSquaredKernel * ExpSine2Kernel
           -> set_parameter_vector() accepts the natural logarithm of this value
         gp_pams[1] = metric = r^2 = lambda**2  -> ExpSquaredKernel(metric=r^2)
           -> set_parameter_vector() accepts the natural logarithm of this value
         gp_pams[2] = Gamma =  1/ (2 omega**2) -> ExpSine2Kernel(gamma, ln_period)
         gp_pams[3] = ln_theta = ln_Period -> ExpSine2Kernel(gamma, ln_period)

        """

        if self.use_HODLR:
            self.gp[dataset.name_ref] = george.GP(
                kernel, solver=george.HODLRSolver, mean=0.00)
            print(' *** USING HODLR *** ')
            print()

        else:
            self.gp[dataset.name_ref] = george.GP(kernel)

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

        self.update_parameter_values(parameter_values)
        
        pass_conditions = self.check_hyperparameter_values(parameter_values)
        if not pass_conditions:
            return pass_conditions
        
        gp_pams = self.convert_val2gp(parameter_values)
        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(dataset.x0, env)

        return self.gp[dataset.name_ref].log_likelihood(dataset.residuals, quiet=True)

    def sample_predict(self, parameter_values, dataset, x0_input=None, return_covariance=False, return_variance=False):

        self.update_parameter_values(parameter_values)
        
        gp_pams = self.convert_val2gp(parameter_values)

        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(dataset.x0, env)
        if x0_input is None:
            return self.gp[dataset.name_ref].predict(dataset.residuals, dataset.x0, return_cov=return_covariance, return_var=return_variance)
        else:
            return self.gp[dataset.name_ref].predict(dataset.residuals, x0_input, return_cov=return_covariance, return_var=return_variance)

    def sample_conditional(self, parameter_values, dataset, x0_input=None):

        self.update_parameter_values(parameter_values)

        gp_pams = self.convert_val2gp(parameter_values)

        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(dataset.x0, env)
        if x0_input is None:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, dataset.x0)
        else:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, x0_input)
