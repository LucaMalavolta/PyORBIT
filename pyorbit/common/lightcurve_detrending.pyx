from pyorbit.subroutines.common import *
from pyorbit.common.abstract_common import *


class CommonLightcurveDetrending(AbstractCommon):
    """
    Inherited class from AbstractCommon

    Attributes:
        :model_class (string): identify the kind of class
        :list_pams: all the possible parameters that can be assigned to a planet are listed here
        :default_bounds: these default boundaries are used when the user does not define them in the yaml file
        :recenter_pams: circular parameters that may need a recentering around the most likely value after the global
            optimization run
        :period_average: variable used only by TRADES
    """

    model_class = 'lightcurve_detrending'

    parameters_dictionary = {
        'coeff_linear': # Orbital period of the planet
            {
                'bounds': [-10., 10.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'coeff_poly': # Orbital period of the planet
            {
                'bounds': [-10., 10.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'coeff_c0': # Orbital period of the planet
            {
                'bounds': [-100000., 100000.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'x_zero': # Orbital period of the planet
            {
                'bounds': [-1000., 1000.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
    }

    default_fixed = {}

    recenter_pams = {}
