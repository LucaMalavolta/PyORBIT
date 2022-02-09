from pyorbit.subroutines.common import \
    np,\
    get_var_exp,\
    get_var_val,\
    get_fix_val,\
    get_2darray_from_val,\
    giveback_priors,\
    nested_sampling_prior_prepare

from pyorbit.subroutines.common import get_var_exp, get_var_log
from pyorbit.subroutines.common import get_var_exp_base2, get_var_log_base2
from pyorbit.subroutines.common import get_var_exp_base10, get_var_log_base10
from pyorbit.subroutines.common import get_var_exp_natural, get_var_log_natural


class AbstractModel(object):

    def __init__(self, model_name, common_ref):
        self.model_name = model_name

        try:
            if len(common_ref) > 0:
                self.common_ref = np.atleast_1d(common_ref).tolist()
            else:
                self.common_ref = []
        except (NameError, TypeError):
            self.common_ref = []

        """ Set of variables shared among all the datasets relying on the same
        physical (i.e., common) objects. Variables belonging to different common
        objects can be listed together.
        When making new objects, always pay attention to avoid duplicate names
        for variables
        """
        self.list_pams_common = set()

        """ Dataset-specific variables must be listed here. A given model will
        be re-computed for each dataset using the corresponding values of the
        variables listed here, while all the other varables will be taken by the
        common objects of reference. As an example, all Gaussian Process object
        will have the same period (listed in list_pams_common) but each dataset
        will be characterized by its own covariance amplitude (listed in
        list_pams_dataset).
        """
        self.list_pams_dataset = set()

        """ For some circular variables it is convenient to recenter the
        boundaries so that tha mode is near the center of the interval defined
        by the boundaries, in such a way the computation of the median and the
        confidence interval are less prone to errors
        """
        self.recenter_pams_dataset = set()

        self.unitary_model = False
        self.normalization_model = False

        self.planet_ref = common_ref
        self.stellar_ref = 'star_parameters'

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
        self.spaces = {}

        self.default_bounds = {}
        self.default_spaces = {}
        self.default_priors = {}

        self.model_conf = None

        self.multivariate_priors = {}
        self.multivariate_vars = {}
        self.multivariate_func = {}
        self.multivariate_med = {}
        self.multivariate_cov = {}

    def initialize_model(self, mc, **kwargs):
        pass

    def change_variable_status(self, mc, **kwargs):

        dataset_vars = kwargs.get('dataset_variables', [])
        for var in dataset_vars:
            self.list_pams_common.discard(var)
            self.list_pams_dataset.update([var])
        common_vars = kwargs.get('common_variables', [])
        for var in common_vars:
            self.list_pams_dataset.discard(var)
            self.list_pams_common.update([var])

    def initialize_model_dataset(self, mc, dataset, **kwargs):
        pass

    def define_special_variable_properties(self,
                                           ndim,
                                           output_lists,
                                           dataset_name,
                                           var):
        return ndim, output_lists, False

    def define_variable_properties(self, ndim, output_lists, dataset_name):
        """[summary]
            Boundaries are defined in this class, where all the dataset-related
            variables are stored. Bounds and parameter index CANNOT be defined
            in the Common class: we don't know a priori which parameters
            will be actually used in the complete model.

        Args:
            ndim ([type]): [description]
            output_lists ([type]): [description]
            dataset_name ([type]): [description]

        Returns:
            [type]: [description]
        """

        self.transformation[dataset_name] = {}
        self.variable_index[dataset_name] = {}
        self.variable_sampler[dataset_name] = {}

        if dataset_name not in self.bounds.keys():
            self.bounds[dataset_name] = {}

        if dataset_name not in self.spaces.keys():
            self.spaces[dataset_name] = {}

        for var in self.list_pams_dataset:
            ndim, output_lists, applied = \
                self.define_special_variable_properties(
                    ndim, output_lists, dataset_name, var)
            if applied:
                continue

            if var not in self.bounds[dataset_name]:
                self.bounds[dataset_name][var] = self.default_bounds[var]

            if var not in self.spaces[dataset_name]:
                self.spaces[dataset_name][var] = self.default_spaces[var]

            if var in self.fix_list[dataset_name]:

                self.transformation[dataset_name][var] = get_fix_val
                self.variable_index[dataset_name][var] = self.nfix
                self.prior_kind[dataset_name][var] = 'None'
                self.prior_pams[dataset_name][var] = []

                # Workaround to preserve compatibility with Python 2.x
                if isinstance(self.fix_list[dataset_name][var],
                              type('string')) \
                        and var in self.default_fixed:
                    self.fixed.append(
                        get_2darray_from_val(self.default_fixed[var])[0])
                else:
                    self.fixed.append(self.fix_list[dataset_name][var][0])

                self.nfix += 1
            else:
                if self.spaces[dataset_name][var] == 'Linear':
                    self.transformation[dataset_name][var] = get_var_val
                    output_lists['bounds'].append(
                        self.bounds[dataset_name][var])
                elif self.spaces[dataset_name][var] == 'Log_Natural':
                    self.transformation[dataset_name][var] = get_var_log_natural
                    output_lists['bounds'].append(
                        np.log(self.bounds[dataset_name][var]))
                elif self.spaces[dataset_name][var] == 'Log_Base2':
                    self.transformation[dataset_name][var] = get_var_exp_base2
                    output_lists['bounds'].append(
                        np.log2(self.bounds[dataset_name][var]))
                elif self.spaces[dataset_name][var] == 'Log_Base10':
                    self.transformation[dataset_name][var] = get_var_exp_base10
                    output_lists['bounds'].append(
                        np.log10(self.bounds[dataset_name][var]))
                elif self.spaces[dataset_name][var] == 'Logarithmic':
                    self.transformation[dataset_name][var] = get_var_exp
                    output_lists['bounds'].append(
                        np.log2(self.bounds[dataset_name][var]))

                if var not in self.prior_pams[dataset_name]:
                    self.prior_kind[dataset_name][var] = \
                        self.default_priors[var][0]
                    self.prior_pams[dataset_name][var] = \
                        self.default_priors[var][1]

                nested_coeff = \
                    nested_sampling_prior_prepare(
                        self.prior_kind[dataset_name][var],
                        output_lists['bounds'][-1],
                        self.prior_pams[dataset_name][var],
                        self.spaces[dataset_name][var])

                output_lists['spaces'].append(self.spaces[dataset_name][var])
                output_lists['priors'].append(
                    [self.prior_kind[dataset_name][var],
                     self.prior_pams[dataset_name][var],
                     nested_coeff])

                self.variable_index[dataset_name][var] = ndim
                self.variable_sampler[dataset_name][var] = ndim
                ndim += 1

        return ndim, output_lists

    def define_special_starting_point(self, starting_point, dataset_name, var):
        return False

    def define_starting_point(self, starting_point, dataset_name):

        if not bool(self.starts):
            return

        for var in self.starts[dataset_name]:
            if self.define_special_starting_point(
                    starting_point, dataset_name, var):
                continue

            if self.spaces[dataset_name][var] == 'Linear':
                start_converted = self.starts[dataset_name][var]
            elif self.spaces[dataset_name][var] == 'Log_Natural':
                start_converted = np.log(self.starts[dataset_name][var])
            elif self.spaces[dataset_name][var] == 'Log_Base2':
                start_converted = np.log2(self.starts[dataset_name][var])
            elif self.spaces[dataset_name][var] == 'Log_Base10':
                start_converted = np.log10(self.starts[dataset_name][var])
            elif self.spaces[dataset_name][var] == 'Logarithmic':
                start_converted = np.log2(self.starts[dataset_name][var])
            starting_point[self.variable_sampler[dataset_name]
                           [var]] = start_converted

    def convert(self, theta, dataset_name):
        variable_value = {}
        # If we need the parameters for the prior, we are not providing any
        # name for the dataset

        for var in self.list_pams_dataset:
            variable_value[var] = self.transformation[dataset_name][var](
                theta, self.fixed, self.variable_index[dataset_name][var])
        return variable_value

    def return_priors(self, theta, dataset_name):
        prior_out = 0.00
        variable_value = self.convert(theta, dataset_name)

        """ Preserving backcompatibility with version 8
        #TODO: to be simplified in the next version
        """

        if getattr(self, 'multivariate_priors', False):
            if len(set(self.multivariate_vars[dataset_name]) & set(self.list_pams_dataset[dataset_name])) > 0:
                multi_var = [variable_value[ii]
                             for ii in self.multivariate_vars[dataset_name]]
                pdf = self.multivariate_func[dataset_name].pdf(multi_var)
                if pdf > 0:
                    prior_out += np.log(
                        self.multivariate_func[dataset_name].pdf(multi_var))
                else:
                    return -np.inf
            else:
                self.multivariate_vars[dataset_name] = []
        else:
            self.multivariate_vars = {dataset_name: []}

        for var in self.list_pams_dataset:

            if var in self.multivariate_vars[dataset_name]:
                continue

            prior_out += giveback_priors(
                self.prior_kind[dataset_name][var],
                self.bounds[dataset_name][var],
                self.prior_pams[dataset_name][var],
                variable_value[var])
        """
        for var in list(set(self.list_pams_dataset) and \
         set(self.prior_pams[dataset_name])):

            prior_out += giveback_priors(self.prior_kind[dataset_name][var],
                                         self.bounds[dataset_name][var],
                                         self.prior_pams[dataset_name][var],
                                         variable_value[var])
        """
        return prior_out

    def index_recenter_bounds(self, dataset_name):
        ind_list = []
        for var in list(set(self.recenter_pams_dataset)
                        & set(self.variable_sampler[dataset_name])):
            ind_list.append(self.variable_sampler[dataset_name][var])

        return ind_list

    def special_index_recenter_bounds(self, dataset_name):
        return []

    def special_fix_population(self, pop_mean, population, dataset_name):
        return population

    def compute(self, variable_value, dataset, x0_input=None):
        return 0.00

    def _subset_transfer_priors(self, mc, dataset, var_original, var_subset):

        self.list_pams_dataset.update([var_subset])

        for common_model in self.common_ref:
            if var_original not in mc.common_models[common_model].list_pams:
                continue

            if var_original in self.bounds[dataset.name_ref]:
                var_update = self.bounds[dataset.name_ref][var_original]
            else:
                var_update = mc.common_models[common_model].default_bounds[var_original]

            self.bounds[dataset.name_ref].update({var_subset: var_update})

            if var_original in self.spaces[dataset.name_ref]:
                var_update = self.spaces[dataset.name_ref][var_original]
            else:
                var_update = mc.common_models[common_model].default_spaces[var_original]

            self.spaces[dataset.name_ref].update({var_subset: var_update})

            if var_original in self.prior_pams[dataset.name_ref]:
                var_update0 = self.prior_kind[dataset.name_ref][var_original]
                var_update1 = self.prior_pams[dataset.name_ref][var_original]
            else:
                var_update0 = mc.common_models[common_model].default_priors[var_original][0]
                var_update1 = mc.common_models[common_model].default_priors[var_original][1]

            self.prior_kind[dataset.name_ref].update({var_subset: var_update0})
            self.prior_pams[dataset.name_ref].update({var_subset: var_update1})

