from abstract_common import *
#from abstract_model import *

class LimbDarkening_1Pam(AbstractCommon):
    ''' This class must be created for each planet in the system
            model_name is the way the planet is identified
    '''

    "Linear LD"
    list_pams = {
        'ld_c1'
    }

    """These default boundaries are used when the user does not define them in the yaml file"""
    default_bounds = {
        'ld_c1': [0.00, 1.00]

    }

    default_spaces = {
        'ld_c1': 'Linear'
    }

    default_priors = {
        'ld_c1': ['Uniform', []]
    }

    ld_ncoeff = 1

    default_fixed = {}
    recenter_pams = {}


class LimbDarkening_2Pam(AbstractCommon):
    ''' This class must be created for each planet in the system
            model_name is the way the planet is identified
    '''

    "2-coefficient LD"
    list_pams = {
        'ld_c1',
        'ld_c2'
    }

    """These default boundaries are used when the user does not define them in the yaml file"""
    default_bounds = {
        'ld_c1': [0.00, 1.00],
        'ld_c2': [0.00, 1.00]
    }

    default_spaces = {
        'ld_c1': 'Linear',
        'ld_c2': 'Linear'
    }

    default_priors = {
        'ld_c1': ['Uniform', []],
        'ld_c2': ['Uniform', []]
    }

    ld_ncoeff = 2

    default_fixed = {}
    recenter_pams = {}


class LimbDarkening_4Pam(AbstractCommon):
    ''' This class must be created for each planet in the system
            model_name is the way the planet is identified
    '''

    "4-coefficient LD"
    list_pams = {
        'ld_c1',
        'ld_c2',
        'ld_c3',
        'ld_c4'
    }

    """These default boundaries are used when the user does not define them in the yaml file"""
    default_bounds = {
        'ld_c1': [0.00, 1.00],
        'ld_c2': [0.00, 1.00],
        'ld_c3': [0.00, 1.00],
        'ld_c4': [0.00, 1.00]
    }

    default_spaces = {
        'ld_c1': 'Linear',
        'ld_c2': 'Linear',
        'ld_c3': 'Linear',
        'ld_c4': 'Linear'
    }

    default_priors = {
        'ld_c1': ['Uniform', []],
        'ld_c2': ['Uniform', []],
        'ld_c3': ['Uniform', []],
        'ld_c4': ['Uniform', []]
    }

    ld_ncoeff = 4

    default_fixed = {}
    recenter_pams = {}


class Batman_LimbDarkening_Linear(LimbDarkening_1Pam):
    model_class = 'batman_ld_linear'
    ld_type = 'linear'


class Batman_LimbDarkening_Quadratic(LimbDarkening_2Pam):
    model_class = 'batman_ld_quadratic'
    ld_type = 'quadratic'


class Batman_LimbDarkening_SquareRoot(LimbDarkening_2Pam):
    model_class = 'batman_ld_square-root'
    ld_type = 'square-root'


class Batman_LimbDarkening_Logarithmic(LimbDarkening_2Pam):
    model_class = 'batman_ld_logarithmic'
    ld_type = 'logarithmic'


class Batman_LimbDarkening_Exponential(LimbDarkening_2Pam):
    model_class = 'batman_ld_exponential'
    ld_type = 'exponential'


class Batman_LimbDarkening_Power2(LimbDarkening_1Pam):
    model_class = 'batman_ld_power2'
    ld_type = 'power2'


class Batman_LimbDarkening_NonLinear(LimbDarkening_4Pam):
    model_class = 'batman_ld_nonlinear'
    ld_type = 'nonlinear'

