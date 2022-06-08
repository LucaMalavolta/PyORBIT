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
        'coeff_linear':
            {
                'bounds': [-10., 10.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'coeff_poly':
            {
                'bounds': [-10., 10.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'coeff_c0':
            {
                'bounds': [-100000., 100000.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'x_zero':
            {
                'bounds': [-1e09, 1e09],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': 'as input',
            },
    }

    for i_pol in range(0,10):
        # Coefficient of the i_pol order of the polynomial
        parameters_dictionary['poly_c'+repr(i_pol)] = {
                'bounds': [-1e06, 1e06],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': 'poly order '+repr(i_pol),
        }

    default_fixed = {}

    recenter_pams = {}
