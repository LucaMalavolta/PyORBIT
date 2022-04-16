from pyorbit.subroutines.common import *
from pyorbit.common.abstract_common import *

class CommonNormalizationFactor(AbstractCommon):
    ''' This class must be created for each planet in the system
            model_name is the way the planet is identified
    '''
    model_class = 'normalization_factor'
    unitary_model = False
    normalization_model = True

    parameters_dictionary = {
        'n_factor':
            {
                'bounds': [1e-06, 1e06],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : 0.000,
                'unit': 'as input',
            },
    }
    recenter_pams = {}
