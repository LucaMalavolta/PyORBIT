from ..classes.common import *


class AbstractModel():
    """

        Comments to be updated

    """
    def __init__(self, model_name, common_ref):
        self.model_name = model_name
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
        self.spaces = {}

        self.model_conf = None

    def initialize_model(self, mc, **kwargs):

        pass

    def setup_dataset(self, dataset, **kwargs):
        pass

    def common_initialization_with_dataset(self, dataset):
        """ Initialization with dataset even if no dataset-specific parameter is present"""
        pass

    def define_special_variable_properties(self, ndim, output_lists, dataset_name, var):
        return ndim, output_lists, False

    def define_variable_properties(self, ndim, output_lists, dataset_name):
        """ Bounds are defined in this class, where all the Planet-related variables are stored
            Bounds and parameter index CANNOT be defined in the Common class: we don't know a priori which parameters
            will be actually used in the complete model.
        """

        self.transformation[dataset_name] = {}
        self.variable_index[dataset_name] = {}
        self.variable_sampler[dataset_name] = {}

        if dataset_name not in self.bounds.keys():
            self.bounds[dataset_name] = {}

        if dataset_name not in self.spaces.keys():
            self.spaces[dataset_name] = {}

        for var in self.list_pams_dataset:

            ndim, output_lists, applied = self.define_special_variable_properties(ndim, output_lists, dataset_name, var)
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
                self.fixed.append(self.fix_list[dataset_name][var][0])
                self.nfix += 1
            else:
                if self.spaces[dataset_name][var] == 'Linear':
                    self.transformation[dataset_name][var] = get_var_val
                    output_lists['bounds'].append(self.bounds[dataset_name][var])

                if self.spaces[dataset_name][var] == 'Logarithmic':
                    self.transformation[dataset_name][var] = get_var_exp
                    output_lists['bounds'].append(np.log2(self.bounds[dataset_name][var]))

                if var not in self.prior_pams[dataset_name]:
                    self.prior_kind[dataset_name][var] = self.default_priors[var][0]
                    self.prior_pams[dataset_name][var] = self.default_priors[var][1]

                output_lists['spaces'].append(self.spaces[dataset_name][var])
                output_lists['priors'].append([self.prior_kind[dataset_name][var], self.prior_pams[dataset_name][var]])

                self.variable_index[dataset_name][var] = ndim
                self.variable_sampler[dataset_name][var] = ndim
                ndim += 1

        return ndim, output_lists

    def define_special_starting_point(self, starting_point, dataset_name, var):
        return False

    def define_starting_point(self, starting_point, dataset_name):
        if not bool(self.starts): return

        print ' ---------> ', bool(self.starts), self.model_name, dataset_name,  self.starts
        for var in self.starts[dataset_name]:

            if self.define_special_starting_point(starting_point, dataset_name, var): continue

            if self.spaces[dataset_name][var] == 'Linear':
                start_converted = self.starts[dataset_name][var]
            if self.spaces[dataset_name][var] == 'Logarithmic':
                start_converted = np.log2(self.starts[dataset_name][var])
            starting_point[self.variable_sampler[dataset_name][var]] = start_converted

    def convert(self, theta, dataset_name):
        variable_value = {}
        # If we need the parameters for the prior, we are not providing any name for the dataset

        for var in self.list_pams_dataset:
            variable_value[var] = self.transformation[dataset_name][var](
                theta, self.fixed, self.variable_index[dataset_name][var])
        return variable_value

    def return_priors(self, theta, dataset_name):
        prior_out = 0.00
        variable_value = self.convert(theta, dataset_name)

        for var in self.list_pams_dataset:
            prior_out += giveback_priors(self.prior_kind[dataset_name][var],
                                         self.bounds[dataset_name][var],
                                         self.prior_pams[dataset_name][var],
                                         variable_value[var])
        """
        for var in list(set(self.list_pams_dataset) and set(self.prior_pams[dataset_name])):

            prior_out += giveback_priors(self.prior_kind[dataset_name][var],
                                         self.bounds[dataset_name][var],
                                         self.prior_pams[dataset_name][var],
                                         variable_value[var])
        """
        return prior_out

    def index_recenter_bounds(self, dataset_name):
        ind_list = []
        for var in list(set(self.recenter_pams_dataset) & set(self.variable_sampler[dataset_name])):
                ind_list.append(self.variable_sampler[dataset_name][var])

        return ind_list

    def special_index_recenter_bounds(self, dataset_name):
        return []

    def special_fix_population(self, pop_mean, population, dataset_name):
        return population

    def compute(self, variable_value, dataset, x0_input=None):
        return np.zeros(dataset.n, dtype=np.double)

    #def initialize(self):
    #    for var in self.list_pams_dataset:
    #        pam_names[self.variable_index[dataset_name][var]] =

