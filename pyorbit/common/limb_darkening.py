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

    def define_special_variable_properties(self, ndim, output_lists, var):
        """

        :param ndim:
        :param output_lists:
        :param var:
        :return:
        """

        if 'ld_c1' in self.fix_list or \
           'ld_c2' in self.fix_list:
            return ndim, output_lists, False

        for var_check in ['ld_c1', 'ld_c2', 'ld_q1', 'ld_q2']:
            if var_check in self.variable_sampler:
                return ndim, output_lists, False

        if self.parametrization[:8] == 'Standard':
            self.transformation['ld_c1'] = get_var_val
            self.variable_index['ld_c1'] = ndim
            self.transformation['ld_c2'] = get_var_val
            self.variable_index['ld_c2'] = ndim + 1
            variable_list = ['ld_c1', 'ld_c2']
        else:
            self.transformation['ld_c1'] = get_2var_c1
            self.variable_index['ld_c1'] = [ndim, ndim + 1]
            self.transformation['ld_c2'] = get_2var_c2
            self.variable_index['ld_c2'] = [ndim, ndim + 1]
            variable_list = ['ld_q1', 'ld_q2']

        for var in variable_list:
            if var not in self.bounds:
                self.bounds[var] = self.default_bounds[var]

            self.spaces[var] = self.default_spaces[var]

            output_lists['bounds'].append(self.bounds[var])

            if var not in self.prior_pams:
                self.prior_kind[var] = self.default_priors[var][0]
                self.prior_pams[var] = self.default_priors[var][1]

            nested_coeff = nested_sampling_prior_prepare(self.prior_kind[var],
                                                         output_lists['bounds'][-1],
                                                         self.prior_pams[var],
                                                         self.spaces[var])

            output_lists['spaces'].append(self.spaces[var])
            output_lists['priors'].append(
                [self.prior_kind[var], self.prior_pams[var], nested_coeff])

            self.variable_sampler[var] = ndim
            ndim += 1

        for var in ['ld_c1', 'ld_c2']:

            if var not in self.bounds:
                self.bounds[var] = self.default_bounds[var]

            if var not in self.prior_pams:

                if var in self.bounds:
                    self.prior_pams[var] = self.bounds[var]
                else:
                    self.prior_pams[var] = self.default_bounds[var]

                self.prior_kind[var] = 'Uniform'

        return ndim, output_lists, True

    def define_special_starting_point(self, starting_point, var_sampler):
        if var_sampler == 'ld_q1' or var_sampler == 'ld_q2':

            if 'ld_c1' in self.starts and 'ld_c2' in self.starts:
                starting_point[self.variable_sampler['ld_q1']] = \
                    self.starts['ld_c1']**2 + self.starts['ld_c2']**2
                starting_point[self.variable_sampler['ld_q2']] = \
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
    model_class = 'ld_linear'
    ld_type = 'linear'


class LimbDarkening_Quadratic(LimbDarkening_2Pam):
    model_class = 'ld_quadratic'
    ld_type = 'quadratic'


class LimbDarkening_SquareRoot(LimbDarkening_2Pam):
    model_class = 'ld_square-root'
    ld_type = 'square-root'


class LimbDarkening_Logarithmic(LimbDarkening_2Pam):
    model_class = 'ld_logarithmic'
    ld_type = 'logarithmic'


class LimbDarkening_Exponential(LimbDarkening_2Pam):
    model_class = 'ld_exponential'
    ld_type = 'exponential'


class LimbDarkening_Power2(LimbDarkening_1Pam):
    model_class = 'ld_power2'
    ld_type = 'power2'


class LimbDarkening_NonLinear(LimbDarkening_4Pam):
    model_class = 'ld_nonlinear'
    ld_type = 'nonlinear'
