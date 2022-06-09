from pyorbit.subroutines.common import *
from pyorbit.common.abstract_common import *

class CommonJitter(AbstractCommon):
    ''' Common offset for datasets in different files
    '''

    model_class = 'common_jitter'

    parameters_dictionary = {
        'jitter':
            {
                'bounds': {},
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.000,
                'unit': 'as input',
            },
        }

    recenter_pams = {}

