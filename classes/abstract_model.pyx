from common import *

class AbstractModel():

    def define_common_special_bounds(self, mc, var):
        return False

    def define_dataset_special_bounds(self, mc, var):
        return False

    def define_bounds(self, mc):
        """ Bounds are defined in this class, where all the Planet-related variables are stored
            Bounds and parameter index CANNOT be defined in the Common class: we don't know a priori which parameters
            will be actually used in the complete model.
        """

        for var in self.list_pams_common:
            '''We check for each parameter (except eccentricity and omega) if the variable is a
                fixed value or a free variable, and move the parameter into the requested space
                Notice that 'e' and 'w' are not yet included in list_pams[pl_name] at this stage
            '''

            if self.define_common_special_bounds(mc, var):
                continue

            if var in self.common_model.fix_list:
                if var not in self.common_model.transformation[var]:
                    self.common_model.transformation[var] = get_fix_val
                    self.common_model.fixed.append(self.common_model.fix_list[var])
                    self.common_model.variable_index[var] = mc.nfix
                    self.common_model.variable_sampler[var] = mc.nfix
                    self.common_model.nfix += 1
            elif var not in self.common_model.transformation:
                '''If no bounds have been specified in the input file, we use the default ones
                    Bounds must be provided in any case to avoid a failure of PyDE '''
                if var in self.common_model.bounds:
                    bounds_tmp = self.common_model.bounds[var]
                else:
                    bounds_tmp = self.common_model.default_bounds[var]

                if self.common_model.list_pams[var] == 'U':
                    self.common_model.transformation[var] = get_var_val
                    mc.bounds_list.append(bounds_tmp)
                elif self.common_model.list_pams[var] == 'LU':
                    self.common_model.transformation[var] = get_var_exp
                    mc.bounds_list.append(np.log2(bounds_tmp))

                self.common_model.variable_index[var] = mc.ndim
                self.common_model.variable_sampler[var] = mc.ndim
                mc.ndim += 1

        ''' Repeating the same procedure for dataset-related variables'''

        for dataset_name, dataset in mc.dataset_dict.items():
            if self.model_name in dataset.models:
                for var in self.list_pams_dataset:

                    if self.define_dataset_special_bounds(mc, var):
                        continue

                    if var in self.fix_list[dataset_name]:
                        self.transformation[dataset_name][var] = get_fix_val
                        self.variable_index[dataset_name][var] = self.nfix
                        self.variable_sampler[dataset_name][var] = self.nfix
                        self.fixed.append(self.fix_list[dataset_name][var])
                        self.nfix += 1
                    else:
                        if self.list_pams_dataset[var] == 'U':
                            self.transformation[dataset_name][var] = get_var_val
                            mc.bounds_list.append(self.bounds[dataset_name][var])
                        if self.list_pams_dataset[var] == 'LU':
                            self.transformation[dataset_name][var] = get_var_exp
                            mc.bounds_list.append(np.log2(self.bounds[dataset_name][var]))
                        self.variable_index[dataset_name][var] = mc.ndim
                        self.variable_sampler[dataset_name][var] = mc.ndim
                        mc.ndim += 1
    """ 
    def initialize(self, mc):

        for name in self.list_pams_common:
            if name in mc.variable_sampler[self.common_ref]:
                mc.pam_names[mc.variable_sampler[self.common_ref][name]] = name

        for dataset_name, dataset in mc.dataset_dict.items():
            for name in self.list_pams_dataset:
                if name in mc.variable_sampler[dataset_name]:
                    mc.pam_names[mc.variable_sampler[dataset_name][name]] = name

            if self.model_name in dataset.models:
                self.define_kernel(dataset)

    """
    def define_common_special_starting_point(self, mc, var):
        return False

    def define_dataset_special_starting_point(self, mc, var):
        return False

    def starting_point(self, mc):

        for var in self.common_model.starts:

            if self.define_common_special_starting_point(mc, var):
                continue

            if self.common_model.list_pams[var] == 'U':
                start_converted = self.common_model.starts[var]
            if self.common_model.list_pams[var] == 'LU':
                start_converted = np.log2(self.common_model.starts[var])
            mc.starting_point[self.common_model.variable_sampler[var]] = start_converted

        for dataset_name, dataset in mc.dataset_dict.items():
            if self.model_name in dataset.models and dataset_name in self.starts:
                for var in self.common_model.starts:

                    if self.define_dataset_special_starting_point(mc, var):
                        continue

                    if self.list_pams_dataset[var] == 'U':
                        start_converted = self.starts[dataset_name][var]
                    if self.list_pams_dataset[var] == 'LU':
                        start_converted = np.log2(self.starts[dataset_name][var])
                    mc.starting_point[self.variable_sampler[dataset_name][var]] = start_converted

    def convert(self, theta, dataset_name=None):
        variable_value = {}
        for var in self.list_pams_common:
            variable_value[var] = self.common_model.transformation[var](
                theta, self.fixed, self.common_model.variable_sampler[var])
        # If we need the parameters for the prior, we are not providing any name for the dataset
        if dataset_name is not None:
            for var in self.list_pams_dataset:
                variable_value[var] = self.transformation[dataset_name][var](
                    theta, self.fixed, self.variable_sampler[dataset_name][var])
        return variable_value

    def return_priors(self, theta, dataset_name):
        prior_out = 0.00
        variable_value = self.convert(theta, dataset_name)

        for var in self.list_pams_dataset:
            if var in self.prior_pams[dataset_name]:
                prior_out += giveback_priors(self.prior_kind[dataset_name][var],
                                             self.prior_pams[dataset_name][var],
                                             variable_value[var])
        return prior_out

