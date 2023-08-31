from pyorbit.subroutines.common import \
    get_var_val,\
    get_2var_c1,\
    get_2var_c2,\
    nested_sampling_prior_prepare

from pyorbit.common.abstract_common import AbstractCommon


class LimbDarkening_1Pam(AbstractCommon):
    ''' This class must be created for each planet in the system
            model_name is the way the planet is identified
    '''

    "Linear LD"

    parameters_dictionary = {
        'ld_c1': # Limb darkening linear coefficient
            {
                'bounds': [0.0, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'days',
            }
    }

    ld_ncoeff = 1

    default_fixed = {}
    recenter_pams = {}


class LimbDarkening_2Pam(AbstractCommon):
    ''' This class must be created for each planet in the system
            model_name is the way the planet is identified
    '''

    "2-coefficient LD"
    parameters_dictionary = {
        'ld_c1': # Limb darkening linear coefficient
            {
                'bounds': [0.0, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'adimensional',
            },
        'ld_c2':# Limb darkening quadratic coefficient
            {
                'bounds': [-1.00, 1.00],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'adimensional',
            },
        'ld_q1': # Limb darkening linear coefficient in Kipping parametrization
            {
                'bounds': [0.0, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'adimensional',
            },
        'ld_q2': # Limb darkening quadratic coefficient in Kipping parametrization
            {
                'bounds': [0.0, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'adimensional',
            },
    }

    ld_ncoeff = 2

    default_fixed = {}
    recenter_pams = {}

    parametrization_list = ['Standard', 'Kipping']

    def __init__(self, *args, **kwargs):
        super(LimbDarkening_2Pam, self).__init__(*args, **kwargs)

        self.parametrization = 'Standard'

    def initialize_model(self, mc, **kwargs):

        self.parametrization =  kwargs.get('parametrization', self.parametrization)
        if self.parametrization in self.parametrization_list:
            print('Using limb darkening coefficient parametrization: ', self.parametrization)
        else:
            print('ERROR in configuration file - limb darkening: parametrization not supported')
            quit()

    def define_derived_parameters(self):

        derived_list = []

        if 'ld_q1' in self.sampler_parameters and  \
            'ld_q2' in self.sampler_parameters:

            pam00_index = self.sampler_parameters['ld_q1']
            pam01_index = self.sampler_parameters['ld_q2']

            try:
                del self.parameter_index['ld_q1']
                del self.parameter_index['ld_q2']
            except:
                pass
            
            if 'ld_c1' not in self.parameter_index:
                self.transformation['ld_c1'] = get_2var_c1
                self.parameter_index['ld_c1'] = [pam00_index, pam01_index]
                derived_list.append('ld_c1')

            if 'ld_c2' not in self.parameter_index:
                self.transformation['ld_c2'] = get_2var_c2
                self.parameter_index['ld_c2'] = [pam00_index, pam01_index]
                derived_list.append('ld_c2')

        for pam in derived_list:
            if pam not in self.bounds:
                self.bounds[pam] = self.default_bounds[pam]

            if pam not in self.prior_pams:

                if pam in self.bounds:
                    self.prior_pams[pam] = self.bounds[pam]
                else:
                    self.prior_pams[pam] = self.default_bounds[pam]

                self.prior_kind[pam] = 'Uniform'

        return

    def define_starting_point_from_derived(self, starting_point, var_sampler):
        if var_sampler == 'ld_q1' or var_sampler == 'ld_q2':

            if 'ld_c1' in self.starts and 'ld_c2' in self.starts:
                starting_point[self.sampler_parameters['ld_q1']] = \
                    self.starts['ld_c1']**2 + self.starts['ld_c2']**2
                starting_point[self.sampler_parameters['ld_q2']] = \
                    self.starts['ld_c1'] / \
                    (self.starts['ld_c1'] + self.starts['ld_c2'])/2.0

            return True
        return False


class LimbDarkening_4Pam(AbstractCommon):
    ''' This class must be created for each planet in the system
            model_name is the way the planet is identified
    '''

    "4-coefficient LD"
    parameters_dictionary = {
        'ld_c1': # Limb darkening linear coefficient
            {
                'bounds': [0.0, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'adimensional',
            },
        'ld_c2':# Limb darkening quadratic coefficient
            {
                'bounds': [0.0, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'adimensional',
            },
        'ld_c3': # Limb darkening linear coefficient in Kipping parametrization
            {
                'bounds': [0.0, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'adimensional',
            },
        'ld_c4': # Limb darkening quadratic coefficient in Kipping parametrization
            {
                'bounds': [0.0, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'adimensional',
            },
    }

    ld_ncoeff = 4

    default_fixed = {}
    recenter_pams = {}


class LimbDarkening_Linear(LimbDarkening_1Pam):
    model_class = 'limb_darkening'
    ld_type = 'linear'


class LimbDarkening_Quadratic(LimbDarkening_2Pam):
    model_class = 'limb_darkening'
    ld_type = 'quadratic'


class LimbDarkening_SquareRoot(LimbDarkening_2Pam):
    model_class = 'limb_darkening'
    ld_type = 'square-root'


class LimbDarkening_Logarithmic(LimbDarkening_2Pam):
    model_class = 'limb_darkening'
    ld_type = 'logarithmic'


class LimbDarkening_Exponential(LimbDarkening_2Pam):
    model_class = 'limb_darkening'
    ld_type = 'exponential'


class LimbDarkening_Power2(LimbDarkening_1Pam):
    model_class = 'limb_darkening'
    ld_type = 'power2'


class LimbDarkening_NonLinear(LimbDarkening_4Pam):
    model_class = 'limb_darkening'
    ld_type = 'nonlinear'
