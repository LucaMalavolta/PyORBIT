from pyorbit.subroutines.common import *
from pyorbit.common.abstract_common import *

class CommonDilutionFactor(AbstractCommon):
    ''' This class must be created for each planet in the system
            model_name is the way the planet is identified
    '''
    model_class = 'dilution_factor'
    unitary_model = True


    parameters_dictionary = {
        'd_factor':
            {
                'bounds': [0.0000, 1.0000],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.000,
                'unit': 'adimensional',
            },
    }

    recenter_pams = {}
