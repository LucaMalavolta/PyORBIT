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
            variable_value[var] = self.transformation[var](theta, self.fixed, self.variable_index[var])
        return variable_value

    def define_starting_point(self):
        starting_point = {}
        for var in list(set(self.starts) and set(self.variable_sampler)):

            add_special = self.define_special_starting_point(var)
            if bool(add_special):
                starting_point.update(add_special)
                continue

            if self.list_pams[var] == 'U':
                start_converted = self.starts[var]
            if self.list_pams[var] == 'LU':
                start_converted = np.log2(self.starts[var])
            starting_point[self.variable_sampler[var]] = start_converted

        return starting_point

        AHI AHI check indexes


    def return_priors(self, theta):
        """ return prior is defined here because, differently from other functions that can be esecuted more than once
        on the same variable,  the prior for a given parameter should be computed and added to the log_chi2 only one """
        prior_out = 0.00
        variable_value = self.convert(theta)
        for var in list(set(self.prior_pams) & set(variable_value)):
            prior_out += giveback_priors(self.prior_kind[var], self.prior_pams[var], variable_value[var])
        return prior_out

    def index_recenter_bounds(self):
        ind_list = []
        for var in list(set(self.recenter_pams) & set(self.variable_sampler)):
                ind_list.append(self.variable_sampler[var])
        return ind_list


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

    recenter_pams ={'f', 'o', 'lN'}

    def define_special_starting_point(self, var):
        '''eccentricity and argument of pericenter require a special treatment
         since they can be provided as fixed individual values or may need to be combined
         in ecosw and esinw if are both free variables'''

        starting_point = {}
        if not (var == "e" or var == "o"):
            return starting_point

        if 'ecoso' in self.variable_sampler and \
                        'esino' in self.variable_sampler:

            if 'e' in self.starts and 'o' in self.starts:
                starting_point[self.variable_sampler['ecoso']] = \
                    np.sqrt(self.starts['e']) * np.cos(self.starts['o'])
                starting_point[self.variable_sampler['esino']] = \
                    np.sqrt(self.starts['e']) * np.sin(self.starts['o'])

            elif 'ecoso' in self.starts and 'esino' in self.starts:
                starting_point[self.variable_sampler['ecoso']] = self.starts['ecoso']
                starting_point[self.variable_sampler['ecoso']] = self.starts['esino']

        return starting_point
