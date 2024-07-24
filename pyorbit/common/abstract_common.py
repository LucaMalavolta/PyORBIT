from pyorbit.subroutines.common import np, OrderedSet
from pyorbit.subroutines.common import get_var_exp, get_var_log
from pyorbit.subroutines.common import get_var_exp_base2, get_var_log_base2
from pyorbit.subroutines.common import get_var_exp_base10, get_var_log_base10
from pyorbit.subroutines.common import get_var_exp_natural, get_var_log_natural
from pyorbit.subroutines.common import get_var_val
from pyorbit.subroutines.common import get_fix_val
from pyorbit.subroutines.common import get_2darray_from_val
from pyorbit.subroutines.common import giveback_priors
from pyorbit.subroutines.common import nested_sampling_prior_prepare
from pyorbit.subroutines.common import get_var_arccosine, get_var_cosine, get_var_arcsine, get_var_sine
from pyorbit.subroutines import constants

class AbstractCommon(object):
    """

        Comments to be updated

    """

    def __init__(self, common_ref):
        self.common_ref = common_ref
        self.planet_ref = common_ref
        self.stellar_ref = 'star_parameters'
        self.multiple_planets = []

        self.list_pams_common = OrderedSet()
        self.list_pams_dataset = OrderedSet()

        self.sampler_parameters = {}

        self.transformation = {}
        self.parameter_index = {}
        self.parameters = {}

        self.bounds = {}
        self.starts = {}
        self.spaces = {}

        self.fix_list = {}
        self.fixed = []
        self.nfix = 0

        self.prior_kind = {}
        self.prior_pams = {}

        self.multivariate_priors = False
        self.multivariate_pams = []
        self.multivariate_func = None
        self.multivariate_med = None
        self.multivariate_cov = None

        """ reverting the new parameter definition style to the old set
            of parameters
        """
        if hasattr(self, 'parameters_dictionary'):
            self.list_pams = OrderedSet()
            self.default_bounds = {}
            self.default_priors = {}
            self.default_spaces = {}
            self.default_fixed = {}

            for par_name, par_dict in self.parameters_dictionary.items():
                self.list_pams.update([par_name])
                self.default_bounds[par_name] = par_dict['bounds']
                self.default_priors[par_name] = par_dict['priors']
                self.default_spaces[par_name] = par_dict['spaces']
                self.default_fixed[par_name] = par_dict['fixed']

    def initialize_model(self, mc, **kwargs):
        pass

    #def define_special_parameter_properties(self, ndim, output_lists, var):
    #    return ndim, output_lists, False

    def define_derived_parameters(self):
        pass

    def define_parameter_properties(self, ndim, output_lists, parameters_list):
        """ Bounds are defined in this class, where all the Planet-related
            parameters are stored
            Bounds and parameter index CANNOT be defined in the Common class:
            we don't know a priori which parameters
            will be actually used in the complete model.
        """

        for pam in list(OrderedSet(parameters_list) & OrderedSet(self.list_pams)):
            """We check for each parameter (except eccentricity and omega) if
                the parameter is a fixed value or a free parameter, and move the
                parameter into the requested spaces. Notice that 'e' and 'w'
                are not yet included in list_pams[pl_name] at this stage
            """


            # ndim, output_lists, applied = \
            #     self.define_special_parameter_properties(
            #        ndim, output_lists, pam)
            #if applied:
            #    continue

            if pam not in self.bounds:
                self.bounds[pam] = self.default_bounds[pam]

            if pam not in self.spaces:
                self.spaces[pam] = self.default_spaces[pam]

            if pam in self.fix_list:
                if pam not in self.transformation:
                    self.transformation[pam] = get_fix_val

                    # Workaround to preserve compatibility with Python 2.x
                    if isinstance(self.fix_list[pam], type('string')) \
                            and pam in self.default_fixed:
                        self.fixed.append(get_2darray_from_val(
                            self.default_fixed[pam])[0])
                    else:
                        self.fixed.append(self.fix_list[pam][0])

                    self.parameter_index[pam] = self.nfix
                    self.nfix += 1
                    self.prior_kind[pam] = 'None'
                    self.prior_pams[pam] = []

            elif pam not in self.transformation:
                """ If no bounds have been specified in the input file, we use
                    the default ones. Bounds must be provided in any case to
                    avoid a failure of PyDE
                """

                if self.spaces[pam] == 'Linear':
                    self.transformation[pam] = get_var_val
                    output_lists['bounds'].append(self.bounds[pam])
                elif self.spaces[pam] == 'Log_Natural':
                    self.transformation[pam] = get_var_exp_natural
                    output_lists['bounds'].append(
                        np.log(self.bounds[pam]))
                elif self.spaces[pam] == 'Log_Base2':
                    self.transformation[pam] = get_var_exp_base2
                    output_lists['bounds'].append(
                        np.log2(self.bounds[pam]))
                elif self.spaces[pam] == 'Log_Base10':
                    self.transformation[pam] = get_var_exp_base10
                    output_lists['bounds'].append(
                        np.log10(self.bounds[pam]))
                elif self.spaces[pam] == 'Logarithmic':
                    self.transformation[pam] = get_var_exp
                    output_lists['bounds'].append(
                        np.log2(self.bounds[pam]))
                elif self.spaces[pam] == 'Sine_Angle':
                    self.transformation[pam] = get_var_arcsine
                    output_lists['bounds'].append(
                        np.sin(self.bounds[pam]*constants.deg2rad))
                elif self.spaces[pam] == 'Cosine_Angle':
                    self.transformation[pam] = get_var_arccosine
                    output_lists['bounds'].append(
                        np.cos(self.bounds[pam]*constants.deg2rad))


                if pam not in self.prior_pams:
                    self.prior_kind[pam] = self.default_priors[pam][0]
                    self.prior_pams[pam] = self.default_priors[pam][1]

                nested_coeff = \
                    nested_sampling_prior_prepare(
                        self.prior_kind[pam],
                        output_lists['bounds'][-1],
                        self.prior_pams[pam],
                        self.spaces[pam])

                output_lists['spaces'].append(self.spaces[pam])
                output_lists['priors'].append(
                    [self.prior_kind[pam], self.prior_pams[pam], nested_coeff])

                self.parameter_index[pam] = ndim
                self.sampler_parameters[pam] = ndim
                ndim += 1

        self.define_derived_parameters()

        return ndim, output_lists

    def convert(self, theta):
        parameter_values = {}
        for pam in self.parameter_index:
            parameter_values[pam] = self.transformation[pam](
                theta, self.fixed, self.parameter_index[pam])
        return parameter_values

    def convert_with_name(self, theta, name):
        parameter_values = {}
        for pam in self.parameter_index:
            parameter_values[name + '_' + pam] = self.transformation[pam](
                theta, self.fixed, self.parameter_index[pam])
        return parameter_values

    def define_starting_point_from_derived(self, starting_point, var):
        return False

    def define_starting_point(self, starting_point):

        for sampler_pam in list(OrderedSet(self.starts)
                                and OrderedSet(self.sampler_parameters)):

            if self.define_starting_point_from_derived(starting_point, sampler_pam):
                continue

            if not self.starts.get(sampler_pam, False): continue

            if self.spaces[sampler_pam] == 'Linear':
                start_converted = self.starts[sampler_pam]
            elif self.spaces[sampler_pam] == 'Log_Natural':
                start_converted = np.log(self.starts[sampler_pam])
            elif self.spaces[sampler_pam] == 'Log_Base2':
                start_converted = np.log2(self.starts[sampler_pam])
            elif self.spaces[sampler_pam] == 'Log_Base10':
                start_converted = np.log10(self.starts[sampler_pam])
            elif self.spaces[sampler_pam] == 'Logarithmic':
                start_converted = np.log2(self.starts[sampler_pam])
            elif self.spaces[sampler_pam] == 'Sine_Angle':
                start_converted = np.arcsin(self.starts[sampler_pam])*constants.rad2deg
            elif self.spaces[sampler_pam] == 'Cosine_Angle':
                start_converted = np.arccos(self.starts[sampler_pam])*constants.rad2deg
            starting_point[self.sampler_parameters[sampler_pam]
                           ] = start_converted


    def return_priors(self, theta):
        """Compute the prior probability for a given set of input parameters

        return_priors is defined in the common models because, differently
        from other functions that can be executed more than once on the same
        parameter, the prior for a given parameter should be computed and added
        to the log_chi2 only once

        Args:
            theta: the set of parameters used by the solver
        Returns:
            prior_out: prior probability, to be added to the posterior prob.
        """

        prior_out = 0.00

        parameter_value = self.convert(theta)

        """ Preserving backcompatibility with version 8
        #TODO: to be simplified in the next version
        """
        if getattr(self, 'multivariate_priors', False):
            multi_var = [parameter_value[ii] for ii in self.multivariate_pams]
            pdf = self.multivariate_func.pdf(multi_var)
            if pdf > 0:
                prior_out += np.log(self.multivariate_func.pdf(multi_var))
            else:
                return -np.inf
        else:
            self.multivariate_pams = []

        for pam in parameter_value:

            if pam in self.multivariate_pams: continue

            prior_out += giveback_priors(self.prior_kind[pam],
                                         self.bounds[pam],
                                         self.prior_pams[pam],
                                         parameter_value[pam])
        return prior_out

    def index_recenter_bounds(self):
        ind_list = []
        for pam in list(OrderedSet(self.recenter_pams) & OrderedSet(self.sampler_parameters)):
            ind_list.append(self.sampler_parameters[pam])
        return ind_list

    def _transfer_priors(self, mc, pam_original, pam_addition):
        #! deprecated module


        if pam_original not in self.list_pams:
            return

        self.list_pams.update([pam_addition])

        if pam_original in self.bounds:
            pam_update = self.bounds[pam_original]
        else:
            pam_update = self.default_bounds[pam_original]

        self.bounds.update({pam_addition: pam_update})

        if pam_original in self.spaces:
            pam_update = self.spaces[pam_original]
        else:
            pam_update = self.default_spaces[pam_original]

        self.spaces.update({pam_addition: pam_update})

        if pam_original in self.prior_pams:
            pam_update0 = self.prior_kind[pam_original]
            pam_update1 = self.prior_pams[pam_original]
        else:
            pam_update0 = self.default_priors[pam_original][0]
            pam_update1 = self.default_priors[pam_original][1]

        self.prior_kind.update({pam_addition: pam_update0})
        self.prior_pams.update({pam_addition: pam_update1})

