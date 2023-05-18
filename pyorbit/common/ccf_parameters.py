from pyorbit.subroutines.common import *
from pyorbit.common.abstract_common import *

class CommonCCFParameters(AbstractCommon):
    ''' This class contains CCF parameters employed in some models, such as the Rossiter McLaughlin Revolutions
    '''
    model_class = 'ccf_parameters'
    unitary_model = False
    normalization_model = False

    parameters_dictionary = {
        'contrast_m':
            {
                'bounds': [-10, 10.00],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.000,
                'unit': 'as input',
            },
        'contrast_q':
            {
                'bounds': [0.0, 1.00],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.000,
                'unit': 'as input',
            },
        'fwhm_m':
            {
                'bounds': [-10.0, 10.00],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.000,
                'unit': 'as input',
            },
        'fwhm_q':
            {
                'bounds': [0.0, 70.00],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.000,
                'unit': 'as input',
            },
        'rv_offset':
            {
                'bounds': [-10.0, 10.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.000,
                'unit': 'as input',
            },
    }
    recenter_pams = {}
