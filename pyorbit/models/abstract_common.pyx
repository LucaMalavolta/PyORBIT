from ..classes.common import *


class AbstractCommon(object):
    """

        Comments to be updated

    """

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

            if var not in self.bounds:
                self.bounds[var] = self.default_bounds[var]

            if var in self.fix_list:
                if var not in self.transformation:
                    self.transformation[var] = get_fix_val
                    self.fixed.append(self.fix_list[var][0])
                    self.variable_index[var] = self.nfix
                    #self.variable_sampler[var] = self.nfix
                    self.nfix += 1
            elif var not in self.transformation:
                '''If no bounds have been specified in the input file, we use the default ones
                    Bounds must be provided in any case to avoid a failure of PyDE '''
                #if var in self.bounds:
                #    bounds_tmp = self.bounds[var]
                #else:
                #    bounds_tmp = self.default_bounds[var]
                #
                #if self.list_pams[var] == 'U':
                #    self.transformation[var] = get_var_val
                #    bounds_list.append(bounds_tmp)
                #elif self.list_pams[var] == 'LU':
                #    self.transformation[var] = get_var_exp
                #    bounds_list.append(np.log2(bounds_tmp))

                if self.list_pams[var] == 'U':
                    self.transformation[var] = get_var_val
                    bounds_list.append(self.bounds[var])
                elif self.list_pams[var] == 'LU':
                    self.transformation[var] = get_var_exp
                    bounds_list.append(np.log2(self.bounds[var]))

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
        """Compute the prior probability for a given set of input parameters

        return_priors is defined in the common models because, differently from other functions that can be executed
        more than once on the same variable, the prior for a given parameter should be computed and added to the
        log_chi2 only once

        Args:
            theta: the set of parameters created by the solver
        Returns:
            prior_out: prior probability, to be added to the posterior probability
        """

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


# Temporary fix to avoid running the GAPS activity analysis again
from .planets import CommonPlanets
from .activity import CommonActivity
