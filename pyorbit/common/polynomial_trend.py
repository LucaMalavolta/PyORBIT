from pyorbit.subroutines.common import *
from pyorbit.common.abstract_common import *

class CommonPolynomialTrend(AbstractCommon):

    model_class = 'polynomial_trend'

    parameters_dictionary = {
        'x_zero': # reference value of the polynomial
            {
                'bounds': [-1e09, 1e09],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': 'as input',
            },
        'x_offset': # reference value of the polynomial
            {
                'bounds': [-1e09, 1e09],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': 'as input',
            },
        'poly_factor': # reference value of the polynomial
            {
                'bounds': [-1e09, 1e09],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': 'as input',
            }
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
