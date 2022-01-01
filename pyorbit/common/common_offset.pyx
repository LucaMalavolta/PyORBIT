from pyorbit.common.abstract_common import *

class CommonOffset(AbstractCommon):
    ''' Common offset for datasets in different files
    '''

    model_class = 'common_offset'
    list_pams = {
        'offset'  # order 1
    }

    default_bounds = {}
    default_spaces = {'offset': 'Linear'}
    default_priors = {'offset': ['Uniform', []]}
    default_fixed = {}
    recenter_pams = {}
