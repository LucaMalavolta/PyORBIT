from common import *

class CorrelationsCommonVariables:
    """Class to model the correlations between RVs and activity indexes"""
    def __init__(self):

        self.bounds = {}
        self.starts = {}

        self.order = 0
        self.order_ind = {}
        """order =0 means no correlation """

        self.default_bounds = np.asarray([-1000.0, 1000.0])

    def add_dataset(self, name_ref):
        # planet_name : name of the dataset, otherwise 'common'

        self.bounds[name_ref] = {}
        self.variables[name_ref] = {}
        self.fix_list[name_ref] = {}
        self.var_list[name_ref] = {}

        self.prior_kind[name_ref] = {}
        self.prior_pams[name_ref] = {}

        self.order[name_ref] = {}

    def define_bounds(self, mc):
        """ Bounds are defined in this class, where all the Planet-related variables are stored"""

        for dataset in mc.dataset_list:
            if 'correlation' in dataset.models:
                n_ord = self.order_ind[dataset.name_ref]
                for n_ord in xrange(1, self.order[dataset.name_ref] + 1):
                    var = 'correlation_' + repr(n_ord)
                    if var in self.fix_list[dataset.name_ref]:
                        self.variables[dataset.name_ref][var] = get_fix_val
                        self.var_list[dataset.name_ref][var] = self.nfix
                        self.fixed.append(self.fix_list[dataset.name_ref][var])
                        self.nfix += 1
                    else:
                        if var in self.bounds[dataset.name_ref][var]:
                            bounds_tmp = self.bounds[dataset.name_ref][var]
                        else:
                            bounds_tmp = self.default_bounds

                        if var in self.list_pams[dataset.name_ref]:
                            if self.list_pams[dataset.name_ref][var] == 'U':
                                self.variables[dataset.name_ref][var] = get_var_val
                                mc.bounds_list.append(bounds_tmp)
                            if self.list_pams[dataset.name_ref][var] == 'LU':
                                self.variables[dataset.name_ref][var] = get_var_exp
                                mc.bounds_list.append(np.log2(bounds_tmp))
                        else:
                            self.list_pams[dataset.name_ref][var] == 'U'
                            self.variables[dataset.name_ref][var] = get_var_val
                            mc.bounds_list.append(bounds_tmp)

                        self.var_list[dataset.name_ref][var] = mc.ndim
                        mc.variable_list[dataset.name_ref][var] = mc.ndim
                        mc.ndim += 1

    def starting_point(self, mc):

        for dataset in mc.dataset_list:
            if 'correlation' in dataset.models and dataset.name_ref in self.starts:
                    for var in self.starts[dataset.name_ref]:
                        if self.list_pams[dataset.name_ref][var] == 'U':
                            start_converted = self.starts[dataset.name_ref][var]
                        if self.list_pams[dataset.name_ref][var] == 'LU':
                            start_converted = np.log2(self.starts[dataset.name_ref][var])
                        mc.starting_point[mc.variable_list[dataset.name_ref][var]] = start_converted

    def convert(self, theta, d_name=None):
        dict_out = {}
        # If we need the parameters for the prior, we are not providing any name for the dataset
        if d_name is not None:
            for key in self.list_pams[d_name]:
                dict_out[key] = self.variables[d_name][key](theta, self.fixed, self.var_list[d_name][key])
        return dict_out

    def return_priors(self, theta, d_name=None):
        prior_out = 0.00
        key_pams = self.convert(theta, d_name)
        if d_name is None:
            for key in self.prior_pams[d_name]:
                prior_out += giveback_priors(self.prior_kind[d_name][key], self.prior_pams[d_name][key], key_pams[key])
        return prior_out


    def compute(self, theta, dataset):
        dict_pams = self.convert(theta)
        coeff = np.zeros(self.order+1)
        for var in self.list_pams:
            coeff[self.order_ind[var]] = dict_pams[var]
        return np.polynomial.polynomial.polyval(dataset.x0, coeff)

    def initialize(self, mc):

        for dataset in mc.dataset_list:
            for name in self.list_pams[dataset.name_ref]:
                if name in mc.variable_list[dataset.name_ref]:
                    mc.pam_names[mc.variable_list[dataset.name_ref][name]] = name

    def print_vars(self, mc, theta):

        for dataset in mc.dataset_list:
            for name in self.list_pams[dataset.name_ref]:
                if name in mc.variable_list[dataset.name_ref]:
                    var = self.variables[dataset.name_ref][name](theta, self.fixed,
                                                               self.var_list[dataset.name_ref][name])
                    print 'GaussianProcess ', dataset.name_ref, name, var, '(', theta[
                            self.var_list[dataset.name_ref][name]], ')'
        print

 CHECK VARIABLE LIST!!!!