from ..classes.common import *


class AbstractModel():

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

        self.model_conf = None

    def setup_dataset(self, dataset, **kwargs):
        pass

    def define_special_variables_bounds(self, ndim, dataset_name, var):
        return ndim, []

    def define_variables_bounds(self, ndim, dataset_name):
        """ Bounds are defined in this class, where all the Planet-related variables are stored
            Bounds and parameter index CANNOT be defined in the Common class: we don't know a priori which parameters
            will be actually used in the complete model.
        """
        bounds_list = []

        self.transformation[dataset_name] = {}
        self.variable_index[dataset_name] = {}
        self.variable_sampler[dataset_name] = {}

        if dataset_name not in self.bounds.keys():
            self.bounds[dataset_name] = {}

        for var in self.list_pams_dataset:

            ndim, bounds_special = self.define_special_variables_bounds(ndim, dataset_name, var)
            if len(bounds_special) > 0:
                bounds_list.extend(bounds_special)
                continue

            if var not in self.bounds[dataset_name]:
                self.bounds[dataset_name][var] = self.default_bounds[var]

            if var in self.fix_list[dataset_name]:
                self.transformation[dataset_name][var] = get_fix_val
                self.variable_index[dataset_name][var] = self.nfix
                # self.variable_sampler[dataset_name][var] = self.nfix
                self.fixed.append(self.fix_list[dataset_name][var][0])
                self.nfix += 1
            else:
                if self.list_pams_dataset[var] == 'U':
                    self.transformation[dataset_name][var] = get_var_val
                    bounds_list.append(self.bounds[dataset_name][var])
                if self.list_pams_dataset[var] == 'LU':
                    self.transformation[dataset_name][var] = get_var_exp
                    bounds_list.append(np.log2(self.bounds[dataset_name][var]))
                self.variable_index[dataset_name][var] = ndim
                self.variable_sampler[dataset_name][var] = ndim
                ndim += 1
        return ndim, bounds_list

    def define_special_starting_point(self, starting_point, dataset_name, var):
        return False

    def define_starting_point(self, starting_point, dataset_name):
        if not bool(self.starts): return

        print ' ---------> ', bool(self.starts), self.model_name, dataset_name,  self.starts
        for var in self.starts[dataset_name]:

            if self.define_special_starting_point(starting_point, dataset_name, var): continue

            if self.list_pams_dataset[var] == 'U':
                start_converted = self.starts[dataset_name][var]
            if self.list_pams_dataset[var] == 'LU':
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

        for var in list(set(self.list_pams_dataset) and set(self.prior_pams[dataset_name])):
            prior_out += giveback_priors(self.prior_kind[dataset_name][var],
                                         self.prior_pams[dataset_name][var],
                                         variable_value[var])
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

