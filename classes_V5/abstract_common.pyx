from common import *


class AbstractCommon:
    ''' This class must be created for each planet in the system
        model_name is the way the planet is identified

    '''
    def __init__(self, common_ref):
        self.common_ref = common_ref
        self.variable_sampler = {}

        self.transformation = {}
        self.variable_index = {}
        self.bounds = {}
        self.variables = {}

        self.starts = {}

        self.fix_list = {}
        self.fixed = []
        self.nfix = 0

        self.prior_kind = {}
        self.prior_pams = {}

    def convert(self, theta):
        variable_value = {}
        for var in self.variable_index:
            variable_value[var] = self.common_model.transformation[var](
                theta, self.fixed, self.common_model.variable_sampler[var])

    def return_priors(self, theta):
        prior_out = 0.00
        variable_value = self.convert(theta)
        for var in self.value_model:
            if var in self.prior_pams:
                prior_out += giveback_priors(self.prior_kind[var],
                                             self.prior_pams[var],
                                             variable_value[var])
        return prior_out


class CommonPlanets(AbstractCommon):
    ''' This class must be created for each planet in the system
            model_name is the way the planet is identified
    '''

    model_class = 'planet'

    ''' all the possible parameters that can be assigned to a planet are listed here'''
    list_pams = {
        'P': 'LU',  # Period, log-uniform prior
        'K': 'LU',  # RV semi-amplitude, log-uniform prior
        'f': 'U',  # RV curve phase, log-uniform
        'e': 'U',  # eccentricity, uniform prior - to be fixed
        'o': 'U',  # argument of pericenter (in radians)
        'M': 'LU',  # Mass in Earth masses
        'i': 'U',  # orbital inclination (in degrees)
        'lN': 'U',  # longitude of ascending node
        'R': 'U',  # planet radius (in units of stellar radii)
        'a': 'U'  # semi-major axis (in units of stellar radii)
    }

    """These default boundaries are used when the user does not define them in the yaml file"""
    default_bounds = {
        'P': [0.4, 100000.0],
        'K': [0.5, 2000.0],
        'f': [0.0, 2 * np.pi],
        'ecoso': [-1.0, 1.0],
        'esino': [-1.0, 1.0],
        'e': [0.0, 1.0],
        'o': [0.0, 2 * np.pi],
        # Used by TTVfast/TRADES
        'M': [0.5, 10000],  # Fix the unit
        'i': [0.0, 180.0],
        'lN': [0.0, 2 * np.pi],
        # Used by BATMAN
        'R': [0.00001, 0.5],  # Fix the unit
        'a': [0.00001, 50.]  # Fix the unit
    }
