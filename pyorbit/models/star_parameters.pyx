from pyorbit.classes.common import *
from pyorbit.models.abstract_common import *
from pyorbit.models.abstract_model import *


class CommonStarParameters(AbstractCommon):
    ''' This class must be created for each planet in the system
            model_name is the way the planet is identified
    '''
    model_class = 'star_parameters'

    list_pams = {
        'radius',  # radius of the star, in Solar radii
        'mass',  # mass of the star, in Solar masses
        'rho',  # density of the star, in Solar density units
        'teff' #effective temperature of the star, in K 
    }

    default_bounds = {
        'radius': [0.0000, 2.0000],
        'mass': [0.0000, 2.0000],
        'rho': [0.0000, 5.0000],
        'teff': [3000., 11000.]
    }

    """ Must be the same parameters as in list_pams, because priors are applied only to _physical_ parameters """
    default_priors = {
        'radius': ['Uniform', []],
        'mass': ['Uniform', []],
        'rho': ['Uniform', []],
        'teff': ['Uniform', []]
    }

    default_spaces = {
        'radius': 'Linear',
        'mass': 'Linear',
        'rho': 'Linear',
        'teff': 'Linear'
    }

    default_fixed = {
        'radius': 1.0000,
        'mass': 1.0000,
        'rho': 1.0000,
        'teff': 5777
    }

    recenter_pams = {}
