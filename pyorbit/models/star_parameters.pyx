from ..classes.common import *
from abstract_common import *

class CommonStarParameters(AbstractCommon):
    ''' This class must be created for each planet in the system
            model_name is the way the planet is identified
    '''
    model_class = 'star_parameters'

    list_pams = {
        'radius',  # radius of the star, in Solar radii
        'mass',  # mass of the star, in Solar masses
        'rho'  # density of the star, in Solar density units
    }

    default_bounds = {
        'radius': [0.0000, 2.0000],
        'mass': [0.0000, 2.0000],
        'rho': [0.0000, 4.0000]
    }

    """ Must be the same parameters as in list_pams, because priors are applied only to _physical_ parameters """
    default_priors = {
        'radius': ['Uniform', []],
        'mass': ['Uniform', []],
        'rho': ['Uniform', []]
    }

    default_spaces = {
        'radius': 'Linear',
        'mass': 'Linear',
        'rho': 'Linear'
    }

    default_fixed = {
        'radius': 1.0000,
        'mass': 1.0000,
        'rho': 1.0000
    }

    recenter_pams = {}
