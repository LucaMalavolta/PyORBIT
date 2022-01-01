from pyorbit.classes.common import *
from pyorbit.common.abstract_common import *

class CommonNormalizationFactor(AbstractCommon):
    ''' This class must be created for each planet in the system
            model_name is the way the planet is identified
    '''
    model_class = 'normalization_factor'
    unitary_model = False
    normalization_model = True

    list_pams = {
        'n_factor',  # normalization factor
    }

    default_bounds = {
        'n_factor': [0.0000010, 100000.0000]
    }

    """ Must be the same parameters as in list_pams, because priors are applied only to _physical_ parameters """
    default_priors = {
        'n_factor': ['Uniform', []]
    }

    default_spaces = {
        'n_factor': 'Logarithmic'
    }

    default_fixed = {
        'n_factor': 1.0000
    }

    recenter_pams = {}
