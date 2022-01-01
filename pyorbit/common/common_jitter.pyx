from pyorbit.classes.common import *
from pyorbit.common.abstract_common import *

class CommonJitter(AbstractCommon):
    ''' Common offset for datasets in different files
    '''

    model_class = 'common_jitter'

    list_pams = {
        'jitter'
    }

    default_bounds = {}
    default_spaces = {'jitter': 'Linear'}
    default_priors = {'jitter': ['Uniform', []]}
    default_fixed = {}

    recenter_pams = {}

