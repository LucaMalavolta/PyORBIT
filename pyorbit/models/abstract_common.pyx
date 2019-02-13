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
        self.variables = {}

        self.bounds = {}
        self.starts = {}
        self.spaces = {}

        self.fix_list = {}
        self.fixed = []
        self.nfix = 0

        self.prior_kind = {}
        self.prior_pams = {}

    def common_initialization_with_dataset(self, dataset):
        """ Sometimes the common variables still need to be initialized with values coming from a datasets"""
        pass

    def define_special_variable_properties(self, ndim, output_lists, var):
        return ndim, output_lists, False

    def define_variable_properties(self, ndim, output_lists, variable_list):
        """ Bounds are defined in this class, where all the Planet-related variables are stored
            Bounds and parameter index CANNOT be defined in the Common class: we don't know a priori which parameters
            will be actually used in the complete model.
        """

        for var in list(set(variable_list) & set(self.list_pams)):
            '''We check for each parameter (except eccentricity and omega) if the variable is a
                fixed value or a free variable, and move the parameter into the requested spaces
                Notice that 'e' and 'w' are not yet included in list_pams[pl_name] at this stage
            '''

            ndim, output_lists, applied = self.define_special_variable_properties(ndim, output_lists, var)
            if applied:
                continue

            if var not in self.bounds:
                self.bounds[var] = self.default_bounds[var]

            if var not in self.spaces:
                self.spaces[var] = self.default_spaces[var]

            if var in self.fix_list:
                if var not in self.transformation:
                    self.transformation[var] = get_fix_val
                    if self.fix_list[var] is 'default' and var in self.default_fixed:
                            self.fixed.append(get_2darray_from_val(self.default_fixed[var])[0])
                    else:
                        self.fixed.append(self.fix_list[var][0])
                    self.variable_index[var] = self.nfix
                    self.nfix += 1
                    self.prior_kind[var] = 'None'
                    self.prior_pams[var] = []

            elif var not in self.transformation:
                '''If no bounds have been specified in the input file, we use the default ones
                    Bounds must be provided in any case to avoid a failure of PyDE '''

                if self.spaces[var] == 'Linear':
                    self.transformation[var] = get_var_val
                    output_lists['bounds'].append(self.bounds[var])

                elif self.spaces[var] == 'Logarithmic':
                    self.transformation[var] = get_var_exp
                    output_lists['bounds'].append(np.log2(self.bounds[var]))

                if var not in self.prior_pams:
                    self.prior_kind[var] = self.default_priors[var][0]
                    self.prior_pams[var] = self.default_priors[var][1]

                nested_coeff = nested_sampling_prior_prepare(self.prior_kind[var],
                                                              output_lists['bounds'][-1],
                                                              self.prior_pams[var],
                                                              self.spaces[var])

                output_lists['spaces'].append(self.spaces[var])
                output_lists['priors'].append([self.prior_kind[var], self.prior_pams[var], nested_coeff])

                self.variable_index[var] = ndim
                self.variable_sampler[var] = ndim
                ndim += 1

        return ndim, output_lists

    def convert(self, theta):
        variable_value = {}
        for var in self.variable_index:
            variable_value[var] = self.transformation[var](theta, self.fixed, self.variable_index[var])
        return variable_value

    def define_special_starting_point(self, starting_point, var):
        return False

    def define_starting_point(self, starting_point):

        for var_sampler in list(set(self.starts) and set(self.variable_sampler)):
            if self.define_special_starting_point(starting_point, var_sampler): continue

            if self.spaces[var_sampler] == 'Linear':
                start_converted = self.starts[var_sampler]
            if self.spaces[var_sampler] == 'Logarithmic':
                start_converted = np.log2(self.starts[var_sampler])
            starting_point[self.variable_sampler[var_sampler]] = start_converted

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

        for var in variable_value:
            prior_out += giveback_priors(self.prior_kind[var],
                                         self.bounds[var],
                                         self.prior_pams[var],
                                         variable_value[var])

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
