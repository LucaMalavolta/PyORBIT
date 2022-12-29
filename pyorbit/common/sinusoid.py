from pyorbit.subroutines.common import *
from pyorbit.common.abstract_common import *

class CommonSinusoid(AbstractCommon):
    """
    Inherited class from AbstractCommon

    """

    model_class = 'sinusoid'

    parameters_dictionary = {
        'sine_period': # Orbital period of the planet
            {
                'bounds': [0.4, 100000.0],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : None,
                'unit': 'days',
            },
        'sine_amp': 
            {
                'bounds': [-1e09, 1e09],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
            },
        'sine_phase':
            {
                'bounds': [0.0, 360.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0000,
            },
        'sine_offset':
            {
                'bounds': [0.0, 360.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0000,
            },
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

    recenter_pams = {'sine_phase', 'sine_offset'}

