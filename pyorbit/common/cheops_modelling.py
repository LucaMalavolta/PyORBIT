from pyorbit.subroutines.common import *
from pyorbit.common.abstract_common import *


class CommonCheopsModelling(AbstractCommon):
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

    model_class = 'cheops_modelling'

    #detrending_params = [
    #    "dfdt",
    #    "d2fdt2",
    #    "dfdbg",
    #    "dfdcontam",
    #    "dfdsmear",
    #    "dfdx",
    #    "dfdy",
    #    "d2fdx2",
    #    "d2fdxdy",
    #    "d2fdy2",
    #    "dfdsinphi",
    #    "dfdcosphi",
    #    "dfdcos2phi",
    #    "dfdsin2phi",
    #    "dfdcos3phi",
    #    "dfdsin3phi",
    #]


    parameters_dictionary = {
        'scale_factor': 
            {
                'bounds': [-5., 20.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 1.00,
                'unit': None,
            },
        'dfdt': 
            {
                'bounds': [-1., 1.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'd2fdt2': 
            {
                'bounds': [-1., 1.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'dfdbg': 
            {
                'bounds': [-1., 1.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'dfdcontam': 
            {
                'bounds': [-1., 1.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'dfdsmear': 
            {
                'bounds': [-1., 1.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'dfdx': 
            {
                'bounds': [-1., 1.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'dfdy': 
            {
                'bounds': [-1., 1.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'd2fdx2': 
            {
                'bounds': [-1., 1.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'd2fdxdy': 
            {
                'bounds': [-1., 1.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'd2fdy2': 
            {
                'bounds': [-1., 1.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'dfdsinphi': 
            {
                'bounds': [-1., 1.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'dfdcosphi': 
            {
                'bounds': [-1., 1.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'dfdcos2phi': 
            {
                'bounds': [-1., 1.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'dfdsin2phi': 
            {
                'bounds': [-1., 1.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'dfdcos3phi': 
            {
                'bounds': [-1., 1.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'dfdsin3phi': 
            {
                'bounds': [-1., 1.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
    }

    default_fixed = {}

    recenter_pams = {}
