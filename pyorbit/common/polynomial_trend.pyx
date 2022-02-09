from pyorbit.subroutines.common import *
from pyorbit.common.abstract_common import *

class CommonPolynomialTrend(AbstractCommon):

    model_class = 'polynomial_trend'

    "polynomial trend up to 10th order"
    list_pams = {
        'x_zero',
        'poly_c0',  # order 0
        'poly_c1',  # order 1
        'poly_c2',  # order 2
        'poly_c3',  # order 3
        'poly_c4',  # order 4
        'poly_c5',  # order 5
        'poly_c6',  # order 6
        'poly_c7',  # order 7
        'poly_c8',  # order 8
        'poly_c9',  # order 9
    }

    """These default boundaries are used when the user does not define them in the yaml file"""
    default_bounds = {
        'x_zero': [-10 ** 9, 10 ** 9],
        'poly_c0': [-10 ** 5, 10 ** 6],
        'poly_c1': [-10 ** 5, 10 ** 6],
        'poly_c2': [-10 ** 6, 10 ** 6],
        'poly_c3': [-10 ** 6, 10 ** 6],
        'poly_c4': [-10 ** 6, 10 ** 6],
        'poly_c5': [-10 ** 6, 10 ** 6],
        'poly_c6': [-10 ** 6, 10 ** 6],
        'poly_c7': [-10 ** 6, 10 ** 6],
        'poly_c8': [-10 ** 6, 10 ** 6],
        'poly_c9': [-10 ** 6, 10 ** 6]
    }

    default_spaces = {
        'x_zero': 'Linear',
        'poly_c0': 'Linear',  # order 1
        'poly_c1': 'Linear',  # order 1
        'poly_c2': 'Linear',  # order 2
        'poly_c3': 'Linear',  # order 3
        'poly_c4': 'Linear',  # order 4
        'poly_c5': 'Linear',  # order 5
        'poly_c6': 'Linear',  # order 6
        'poly_c7': 'Linear',  # order 7
        'poly_c8': 'Linear',  # order 8
        'poly_c9': 'Linear',  # order 9
    }

    default_priors = {
        'x_zero': ['Uniform', []],
        'poly_c0': ['Uniform', []],
        'poly_c1': ['Uniform', []],
        'poly_c2': ['Uniform', []],
        'poly_c3': ['Uniform', []],
        'poly_c4': ['Uniform', []],
        'poly_c5': ['Uniform', []],
        'poly_c6': ['Uniform', []],
        'poly_c7': ['Uniform', []],
        'poly_c8': ['Uniform', []],
        'poly_c9': ['Uniform', []]
    }

    default_fixed = {}

    recenter_pams = {}
