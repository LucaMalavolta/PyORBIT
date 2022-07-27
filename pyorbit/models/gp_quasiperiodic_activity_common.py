from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *

try:
    import george
except:
    pass


class GaussianProcess_QuasiPeriodicActivity_Common(AbstractModel):
    ''' Three parameters out of four are the same for all the datasets, since they are related to
    the properties of the physical process rather than the observed effects on a dataset
     From Grunblatt+2015, Affer+2016
     - theta: is usually related to the rotation period of the star( or one of its harmonics);
     - lambda: is the correlation decay timescale, and it can be related to the lifetime of the active regions.
     - omega: is the length scale of the periodic component, and can be linked to the size evolution of the active regions;
     - h: represents the amplitude of the correlations '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        try:
            import george
        except ImportError:
            print("ERROR: george not installed, this will not work")
            quit()

        self.model_class = 'gp_quasiperiodic_common'
        self.internal_likelihood = True
        self.delayed_lnlk_computation = True

        self.list_pams_common = {
            'Prot',  # Rotational period of the star
            'Pdec',  # Decay timescale of activity
            'Oamp',  # Granulation of activity
            'Hamp'  # Amplitude of the signal in the covariance matrix
        }

        self.n_pams = 4

        """ Indexing is determined by the way the kernel is constructed, so it
        is specific of the Model and not of the Common class"""
        self.gp_pams_index = {
            'Hamp': 0,  # amp2
            'Pdec': 1,  # metric
            'Oamp': 2,  # gamma
            'Prot': 3  # ln_P
        }

        self.gp = {}
        self.internal_dataset = {'x0': [], 'yr': [], 'ej': []}
        self.internal_variable_value = None
        self.internal_gp_pams = None
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
            'Pdec': np.exp(input_pams[self.gp_pams_index['Pdec']] / 2.0),
            'Oamp': np.sqrt(1. / (2.*input_pams[self.gp_pams_index['Oamp']])),
            'Prot': np.exp(input_pams[self.gp_pams_index['Prot']])
        }

    def initialize_model(self, mc, **kwargs):

        if 'use_HODLR' in kwargs:
            self.use_HODLR = kwargs['use_HODLR']

        if kwargs.get('hyperparameters_condition', False):
            self.hyper_condition = self._hypercond_01
        else:
            self.hyper_condition = self._hypercond_00

        if kwargs.get('rotation_decay_condition', False):
            self.rotdec_condition = self._hypercond_02
        else:
            self.rotdec_condition = self._hypercond_00

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
            self.gp = george.GP(kernel, solver=george.HODLRSolver, mean=0.00)
            print(' *** USING HODLR *** ')
            print()
        else:
            self.gp = george.GP(kernel)
        # self.gp = george.GP(kernel, solver=george.HODLRSolver, mean=0.00)

        """ I've decided to add the jitter in quadrature instead of using a constant kernel to allow the use of
        different / selective jitter within the dataset
        """
        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)

        self.gp.compute(dataset.x0, env)
        # Temporary patch

        return

    def add_internal_dataset(self, variable_value, dataset, reset_status):
        if not reset_status:
            self.internal_dataset['x0'] = []
            self.internal_dataset['yr'] = []
            self.internal_dataset['ej'] = []
            self.internal_variable_value = variable_value
            self.internal_gp_pams = self.convert_val2gp(variable_value)

        self.internal_dataset['x0'].extend(dataset.x0)
        self.internal_dataset['yr'].extend(dataset.residuals)
        self.internal_dataset['ej'].extend(
            np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0))

    def lnlk_compute(self):
        """ 2 steps:
           1) theta parameters must be converted in physical units (e.g. from logarithmic to linear spaces)
           2) physical values must be converted to {\tt george} input parameters
        """
        if not self.hyper_condition(self.internal_variable_value):
            return -np.inf
        if not self.rotdec_condition(self.internal_variable_value):
            return -np.inf

        self.gp.set_parameter_vector(self.internal_gp_pams)
        self.gp.compute(
            self.internal_dataset['x0'], self.internal_dataset['ej'])
        return self.gp.log_likelihood(self.internal_dataset['yr'], quiet=True)

    def sample_predict(self, dataset, x0_input=None, return_covariance=False, return_variance=False):

        self.gp.set_parameter_vector(self.internal_gp_pams)
        self.gp.compute(
            self.internal_dataset['x0'], self.internal_dataset['ej'])

        if x0_input is None:
            return self.gp.predict(self.internal_dataset['yr'], dataset.x0, return_cov=return_covariance, return_var=return_variance)
        else:
            return self.gp.predict(self.internal_dataset['yr'], x0_input, return_cov=return_covariance, return_var=return_variance)

    def sample_conditional(self, dataset, x0_input=None):

        self.gp.set_parameter_vector(self.internal_gp_pams)
        self.gp.compute(
            self.internal_dataset['x0'], self.internal_dataset['ej'])

        if x0_input is None:
            return self.gp.sample_conditional(self.internal_dataset['yr'], dataset.x0)
        else:
            return self.gp.sample_conditional(self.internal_dataset['yr'], x0_input)

    @staticmethod
    def _hypercond_00(variable_value):
        #Condition from Rajpaul 2017, Rajpaul+2021
        return True

    @staticmethod
    def _hypercond_01(variable_value):
        # Condition from Rajpaul 2017, Rajpaul+2021
        # Taking into account that Pdec^2 = 2*lambda_2^2
        return variable_value['Pdec']**2 > (3. / 4. / np.pi) * variable_value['Oamp']**2 * variable_value['Prot']**2 

    @staticmethod
    def _hypercond_02(variable_value):
        #Condition on Rotation period and decay timescale
        return variable_value['Pdec'] > 2. * variable_value['Prot']
