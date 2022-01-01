from pyorbit.classes.common import *
from pyorbit.common.abstract_common import *

class CommonDilutionFactor(AbstractCommon):
    ''' This class must be created for each planet in the system
            model_name is the way the planet is identified
    '''
    model_class = 'dilution_factor'
    unitary_model = True

    list_pams = {
        'd_factor',  # dilution factor
    }

    default_bounds = {
        'd_factor': [0.0000, 1.0000] # default boundaries for the parameter
    }

    default_priors = {
        'd_factor': ['Uniform', []] # Uniform prior (within the boundaries)
    }

    default_spaces = {
        'd_factor': 'Linear' # parameters is sampled in the linear spaced (i.e. not  Logarithmic)
    }

    default_fixed = {
        'd_factor': 0.0000 # default value when the variable is set to fixed
    }

    recenter_pams = {}
