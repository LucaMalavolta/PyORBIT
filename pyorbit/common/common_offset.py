from pyorbit.common.abstract_common import *

class CommonOffset(AbstractCommon):
    ''' Common offset for datasets in different files
    '''

    model_class = 'common_offset'

    parameters_dictionary = {
        'offset':
            {
                'bounds': {},
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.000,
                'unit': 'as input',
            },
    }

    recenter_pams = {}
