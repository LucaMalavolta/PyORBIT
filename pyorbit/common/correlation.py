from pyorbit.subroutines.common import *
from pyorbit.common.abstract_common import *


class CommonCorrelation(AbstractCommon):
    """
    Inherited class from AbstractCommon

    Attributes:
        :model_class (string): identify the kind of class
    """

    model_class = 'correlation'

    parameters_dictionary = {
        'x_zero': # Orbital period of the planet
            {
                'bounds': [-10 ** 9, 10 ** 9],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'corr_c0': # Orbital period of the planet
            {
                'bounds': [-10 ** 5, 10 ** 6],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'corr_c1': # Orbital period of the planet
            {
                'bounds': [-10 ** 5, 10 ** 6],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'corr_c2': # Orbital period of the planet
            {
                'bounds': [-10 ** 5, 10 ** 6],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'corr_c3': # Orbital period of the planet
            {
                'bounds': [-10 ** 5, 10 ** 6],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'corr_c4': # Orbital period of the planet
            {
                'bounds': [-10 ** 5, 10 ** 6],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'corr_c5': # Orbital period of the planet
            {
                'bounds': [-10 ** 5, 10 ** 6],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'corr_c6': # Orbital period of the planet
            {
                'bounds': [-10 ** 5, 10 ** 6],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'corr_c7': # Orbital period of the planet
            {
                'bounds': [-10 ** 5, 10 ** 6],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'corr_c8': # Orbital period of the planet
            {
                'bounds': [-10 ** 5, 10 ** 6],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'corr_c9': # Orbital period of the planet
            {
                'bounds': [-10 ** 5, 10 ** 6],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'corr_c10': # Orbital period of the planet
            {
                'bounds': [-10 ** 5, 10 ** 6],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
    }

    recenter_pams = {}




class CommonComplexCorrelation(CommonCorrelation):
    """
    Inherited class from AbstractCommon

    Attributes:
        :model_class (string): identify the kind of class
    """
    model_class = 'complex_correlation'
