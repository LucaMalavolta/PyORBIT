from common import *


class AbstractCommon(object):
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

    def define_special_variables_bounds(self, ndim, var):
        return ndim, []

    def define_variables_bounds(self, ndim, variable_list):
        """ Bounds are defined in this class, where all the Planet-related variables are stored
            Bounds and parameter index CANNOT be defined in the Common class: we don't know a priori which parameters
            will be actually used in the complete model.
        """
        bounds_list = []
        for var in variable_list:
            '''We check for each parameter (except eccentricity and omega) if the variable is a
                fixed value or a free variable, and move the parameter into the requested space
                Notice that 'e' and 'w' are not yet included in list_pams[pl_name] at this stage
            '''

            ndim, bounds_special = self.define_special_variables_bounds(ndim, var)
            if len(bounds_special) > 0:
                bounds_list.extend(bounds_special)
                continue

            if var in self.fix_list:
                if var not in self.transformation:
                    self.transformation[var] = get_fix_val
                    self.fixed.append(self.fix_list[var][0])
                    self.variable_index[var] = self.nfix
                    self.variable_sampler[var] = self.nfix
                    self.nfix += 1
            elif var not in self.transformation:
                '''If no bounds have been specified in the input file, we use the default ones
                    Bounds must be provided in any case to avoid a failure of PyDE '''
                if var in self.bounds:
                    bounds_tmp = self.bounds[var]
                else:
                    bounds_tmp = self.default_bounds[var]

                if self.list_pams[var] == 'U':
                    self.transformation[var] = get_var_val
                    bounds_list.append(bounds_tmp)
                elif self.list_pams[var] == 'LU':
                    self.transformation[var] = get_var_exp
                    bounds_list.append(np.log2(bounds_tmp))

                self.variable_index[var] = ndim
                self.variable_sampler[var] = ndim
                ndim += 1

        return ndim, bounds_list

    def convert(self, theta):
        variable_value = {}
        for var in self.variable_index:
            variable_value[var] = self.transformation[var](theta, self.fixed, self.variable_index[var])
        return variable_value

    def define_special_starting_point(self, starting_point, var):
        return False

    def define_starting_point(self, starting_point):

        for var in list(set(self.starts) and set(self.variable_sampler)):

            if self.define_special_starting_point(starting_point, var): continue

            if self.list_pams[var] == 'U':
                start_converted = self.starts[var]
            if self.list_pams[var] == 'LU':
                start_converted = np.log2(self.starts[var])
            starting_point[self.variable_sampler[var]] = start_converted

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

    def special_index_recenter_bounds(self):
        return []

    def special_fix_population(self, pop_mean, population):
        return population


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

    recenter_pams = {'f', 'o', 'lN'}

    """ Variable used by trades"""
    period_average = None

    def define_special_variables_bounds(self, ndim, var):

        bounds_list = []

        if var == 'P':
            if var in self.fix_list:
                self.period_average = self.fix_list['P'][0]
            else:
                if var in self.bounds:
                    bounds_tmp = self.bounds[var]
                else:
                    bounds_tmp = self.default_bounds[var]
                self.period_average = np.average(bounds_tmp)
            return ndim, bounds_list

        if not(var == "e" or var == "o"):
            return ndim, bounds_list

        if 'e' in self.fix_list or \
           'o' in self.fix_list:
            return bounds_list

        if 'coso' in self.variable_sampler or \
            'esino' in self.variable_sampler:
            return ndim, bounds_list

        self.transformation['e'] = get_2var_e
        self.variable_index['e'] = [ndim, ndim + 1]
        self.transformation['o'] = get_2var_o
        self.variable_index['o'] = [ndim, ndim + 1]

        self.variable_sampler['ecoso'] = ndim
        self.variable_sampler['esino'] = ndim + 1
        bounds_list.append(self.default_bounds['ecoso'])
        bounds_list.append(self.default_bounds['esino'])
        ndim += 2

        return ndim, bounds_list

    def define_special_starting_point(self, starting_point, var):
        '''eccentricity and argument of pericenter require a special treatment
         since they can be provided as fixed individual values or may need to be combined
         in ecosw and esinw if are both free variables'''

        if not (var == "e" or var == "o"):
            return False

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

        return True

    def special_fix_population(self, population):

        n_pop = np.size(population, axis=0)
        if 'esino' in self.variable_sampler and \
           'ecoso' in self.variable_sampler:
                esino_list = self.variable_sampler['esino']
                ecoso_list = self.variable_sampler['ecoso']
                e_pops = population[:, esino_list] ** 2 + population[:, ecoso_list] ** 2
                o_pops = np.arctan2(population[:, esino_list], population[:, ecoso_list], dtype=np.double)
                # e_mean = (self[planet_name].bounds['e'][0] +
                # self[planet_name].bounds['e'][1]) / 2.
                for ii in xrange(0, n_pop):
                    if not self.bounds['e'][0] + 0.02 <= e_pops[ii] < \
                                    self.bounds['e'][1] - 0.02:
                        e_random = np.random.uniform(self.bounds['e'][0],
                                                     self.bounds['e'][1])
                        population[ii, esino_list] = np.sqrt(e_random) * np.sin(o_pops[ii])
                        population[ii, ecoso_list] = np.sqrt(e_random) * np.cos(o_pops[ii])

