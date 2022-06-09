from pyorbit.subroutines.common import np
from pyorbit.subroutines.common import get_var_exp, get_var_log
from pyorbit.subroutines.common import get_var_exp_base2, get_var_log_base2
from pyorbit.subroutines.common import get_var_exp_base10, get_var_log_base10
from pyorbit.subroutines.common import get_var_exp_natural, get_var_log_natural
from pyorbit.subroutines.common import get_var_val
from pyorbit.subroutines.common import get_fix_val
from pyorbit.subroutines.common import get_2darray_from_val
from pyorbit.subroutines.common import giveback_priors
from pyorbit.subroutines.common import nested_sampling_prior_prepare


class AbstractCommon(object):
    """

        Comments to be updated

    """

    def __init__(self, common_ref):
        self.common_ref = common_ref
        self.planet_ref = common_ref
        self.stellar_ref = 'star_parameters'

        self.list_pams_common = set()
        self.list_pams_dataset = set()

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

        self.multivariate_priors = False
        self.multivariate_vars = []
        self.multivariate_func = None
        self.multivariate_med = None
        self.multivariate_cov = None

        """ reverting the new parameter definition style to the old set
            of variables
        """
        if hasattr(self, 'parameters_dictionary'):
            self.list_pams = set()
            self.default_bounds = {}
            self.default_priors = {}
            self.default_spaces = {}
            self.default_fixed = {}

            for var_name, var_dict in self.parameters_dictionary.items():
                self.list_pams.update([var_name])
                self.default_bounds[var_name] = var_dict['bounds']
                self.default_priors[var_name] = var_dict['priors']
                self.default_spaces[var_name] = var_dict['spaces']
                self.default_fixed[var_name] = var_dict['fixed']

    def define_special_variable_properties(self, ndim, output_lists, var):
        return ndim, output_lists, False

    def define_variable_properties(self, ndim, output_lists, variable_list):
        """ Bounds are defined in this class, where all the Planet-related
            variables are stored
            Bounds and parameter index CANNOT be defined in the Common class:
            we don't know a priori which parameters
            will be actually used in the complete model.
        """

        for var in list(set(variable_list) & set(self.list_pams)):
            """We check for each parameter (except eccentricity and omega) if
                the variable is a fixed value or a free variable, and move the
                parameter into the requested spaces. Notice that 'e' and 'w'
                are not yet included in list_pams[pl_name] at this stage
            """

            ndim, output_lists, applied = \
                self.define_special_variable_properties(
                    ndim, output_lists, var)
            if applied:
                continue

            if var not in self.bounds:
                self.bounds[var] = self.default_bounds[var]

            if var not in self.spaces:
                self.spaces[var] = self.default_spaces[var]

            if var in self.fix_list:
                if var not in self.transformation:
                    self.transformation[var] = get_fix_val

                    # Workaround to preserve compatibility with Python 2.x
                    if isinstance(self.fix_list[var], type('string')) \
                            and var in self.default_fixed:
                        self.fixed.append(get_2darray_from_val(
                            self.default_fixed[var])[0])
                    else:
                        self.fixed.append(self.fix_list[var][0])

                    self.variable_index[var] = self.nfix
                    self.nfix += 1
                    self.prior_kind[var] = 'None'
                    self.prior_pams[var] = []

            elif var not in self.transformation:
                """ If no bounds have been specified in the input file, we use
                    the default ones. Bounds must be provided in any case to
                    avoid a failure of PyDE
                """

                if self.spaces[var] == 'Linear':
                    self.transformation[var] = get_var_val
                    output_lists['bounds'].append(self.bounds[var])
                elif self.spaces[var] == 'Log_Natural':
                    self.transformation[var] = get_var_log_natural
                    output_lists['bounds'].append(
                        np.log(self.bounds[var]))
                elif self.spaces[var] == 'Log_Base2':
                    self.transformation[var] = get_var_exp_base2
                    output_lists['bounds'].append(
                        np.log2(self.bounds[var]))
                elif self.spaces[var] == 'Log_Base10':
                    self.transformation[var] = get_var_exp_base10
                    output_lists['bounds'].append(
                        np.log10(self.bounds[var]))
                elif self.spaces[var] == 'Logarithmic':
                    self.transformation[var] = get_var_exp
                    output_lists['bounds'].append(
                        np.log2(self.bounds[var]))

                if var not in self.prior_pams:
                    self.prior_kind[var] = self.default_priors[var][0]
                    self.prior_pams[var] = self.default_priors[var][1]

                nested_coeff = \
                    nested_sampling_prior_prepare(
                        self.prior_kind[var],
                        output_lists['bounds'][-1],
                        self.prior_pams[var],
                        self.spaces[var])

                output_lists['spaces'].append(self.spaces[var])
                output_lists['priors'].append(
                    [self.prior_kind[var], self.prior_pams[var], nested_coeff])

                self.variable_index[var] = ndim
                self.variable_sampler[var] = ndim
                ndim += 1

        return ndim, output_lists

    def convert(self, theta):
        variable_value = {}
        for var in self.variable_index:
            variable_value[var] = self.transformation[var](
                theta, self.fixed, self.variable_index[var])
        return variable_value

    def define_special_starting_point(self, starting_point, var):
        return False

    def define_starting_point(self, starting_point):

        for var_sampler in list(set(self.starts)
                                and set(self.variable_sampler)):

            if self.define_special_starting_point(starting_point, var_sampler):
                continue

            if not self.starts.get(var_sampler, False): continue

            if self.spaces[var_sampler] == 'Linear':
                start_converted = self.starts[var_sampler]
            elif self.spaces[var_sampler] == 'Log_Natural':
                start_converted = np.log(self.starts[var_sampler])
            elif self.spaces[var_sampler] == 'Log_Base2':
                start_converted = np.log2(self.starts[var_sampler])
            elif self.spaces[var_sampler] == 'Log_Base10':
                start_converted = np.log10(self.starts[var_sampler])
            elif self.spaces[var_sampler] == 'Logarithmic':
                start_converted = np.log2(self.starts[var_sampler])
            starting_point[self.variable_sampler[var_sampler]
                           ] = start_converted

    def return_priors(self, theta):
        """Compute the prior probability for a given set of input parameters

        return_priors is defined in the common models because, differently
        from other functions that can be executed more than once on the same
        variable, the prior for a given parameter should be computed and added
        to the log_chi2 only once

        Args:
            theta: the set of parameters used by the solver
        Returns:
            prior_out: prior probability, to be added to the posterior prob.
        """

        prior_out = 0.00

        variable_value = self.convert(theta)

        """ Preserving backcompatibility with version 8
        #TODO: to be simplified in the next version
        """
        if getattr(self, 'multivariate_priors', False):
            multi_var = [variable_value[ii] for ii in self.multivariate_vars]
            pdf = self.multivariate_func.pdf(multi_var)
            if pdf > 0:
                prior_out += np.log(self.multivariate_func.pdf(multi_var))
            else:
                return -np.inf
        else:
            self.multivariate_vars = []

        for var in variable_value:

            if var in self.multivariate_vars: continue

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

    def _transfer_priors(self, mc, var_original, var_addition):


        if var_original not in self.list_pams:
            return

        self.list_pams.update([var_addition])

        if var_original in self.bounds:
            var_update = self.bounds[var_original]
        else:
            var_update = self.default_bounds[var_original]

        self.bounds.update({var_addition: var_update})

        if var_original in self.spaces:
            var_update = self.spaces[var_original]
        else:
            var_update = self.default_spaces[var_original]

        self.spaces.update({var_addition: var_update})

        if var_original in self.prior_pams:
            var_update0 = self.prior_kind[var_original]
            var_update1 = self.prior_pams[var_original]
        else:
            var_update0 = self.default_priors[var_original][0]
            var_update1 = self.default_priors[var_original][1]

        self.prior_kind.update({var_addition: var_update0})
        self.prior_pams.update({var_addition: var_update1})

