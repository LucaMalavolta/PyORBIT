from common import *


class CorrelationsCommonVariables:
    """Class to model the correlations between RVs and activity indexes"""
    def __init__(self):

        self.bounds = {}
        self.starts = {}

        self.variables = {}
        self.var_list = {}
        self.fix_list = {}

        self.prior_kind = {}
        self.prior_pams = {}

        self.fixed = []
        self.nfix = 0
        self.list_pams = {}

        self.order = {}
        self.order_ind = {}
        """order =0 means no correlation """

        self.default_bounds = np.asarray([-1000000.0, 1000000.0])

        self.x_vals = {}
        self.x_mask = {}

    def add_dataset(self, name_ref):

        self.bounds[name_ref] = {}
        self.starts[name_ref] = {}

        self.variables[name_ref] = {}
        self.fix_list[name_ref] = {}
        self.var_list[name_ref] = {}

        self.prior_kind[name_ref] = {}
        self.prior_pams[name_ref] = {}

        self.list_pams[name_ref] = {}

        self.order[name_ref] = {}
        self.order_ind[name_ref] = {}

        self.x_vals[name_ref] = {}
        self.x_mask[name_ref] = {}

    def add_associated_dataset(self, mc, name_ref, name_asc, threshold=0.001):
        """To each dataset the correlated dataset is associated"""
        self.bounds[name_ref][name_asc] = {}
        self.variables[name_ref][name_asc] = {}
        self.fix_list[name_ref][name_asc] = {}
        self.var_list[name_ref][name_asc] = {}

        self.prior_kind[name_ref][name_asc] = {}
        self.prior_pams[name_ref][name_asc] = {}

        self.list_pams[name_ref][name_asc] = {}

        self.order[name_ref][name_asc] = 1
        self.order_ind[name_ref][name_asc] = {}

        self.x_vals[name_ref][name_asc] = np.zeros(mc.dataset_dict[name_ref].n, dtype=np.double)
        self.x_mask[name_ref][name_asc] = np.zeros(mc.dataset_dict[name_ref].n, dtype=bool)

        """ HERE: we must associated the data from name_asc dataset to the one from name_ref
            remove that part from dataset.pyx
            Add a None option for the dataset
            Fix input_parser to accomodate the new changes
            Jitter must not be included in the analysis, but how about the offset? 
            Or maybe I should just leave the zero point of the polynomial fit free?
        """

        for i_date, v_date in enumerate(mc.dataset_dict[name_asc].x):
            match = np.where(np.abs(v_date-mc.dataset_dict[name_ref].x) < threshold)[0]
            self.x_vals[name_ref][name_asc][match] = mc.dataset_dict[name_asc].x[i_date]
            self.x_mask[name_ref][name_asc][match] = True

    def define_bounds(self, mc):
        """ Bounds are defined in this class, where all the Planet-related variables are stored"""

        for name_ref in self.list_pams:
            for name_asc in self.list_pams[name_ref]:
                for n_ord in xrange(0, self.order[name_ref][name_asc] + 1):
                    var = 'correlation_' + repr(n_ord)
                    self.order_ind[name_ref][name_asc][var] = n_ord
                    if var in self.fix_list[name_ref][name_asc]:
                        self.variables[name_ref][name_asc][var] = get_fix_val
                        self.var_list[name_ref][name_asc][var] = self.nfix
                        self.fixed.append(self.fix_list[name_ref][name_asc][var])
                        self.nfix += 1
                    else:
                        if var in self.bounds[name_ref][name_asc]:
                            bounds_tmp = self.bounds[name_ref][name_asc][var]
                        else:
                            bounds_tmp = self.default_bounds

                        if var in self.list_pams[name_ref][name_asc]:
                            if self.list_pams[name_ref][name_asc][var] == 'U':
                                self.variables[name_ref][name_asc][var] = get_var_val
                                mc.bounds_list.append(bounds_tmp)
                            if self.list_pams[name_ref][name_asc][var] == 'LU':
                                self.variables[name_ref][name_asc][var] = get_var_exp
                                mc.bounds_list.append(np.log2(bounds_tmp))
                        else:
                            self.list_pams[name_ref][name_asc][var] = 'U'
                            self.variables[name_ref][name_asc][var] = get_var_val
                            mc.bounds_list.append(bounds_tmp)

                        self.var_list[name_ref][name_asc][var] = mc.ndim
                        mc.variable_list[name_ref][name_asc][var] = mc.ndim
                        mc.ndim += 1

    def starting_point(self, mc):

        for name_ref in self.list_pams:
            for name_asc in self.starts[name_ref]:
                if name_ref not in self.starts[name_ref]: continue

                for name_var in self.starts[name_ref][name_asc]:
                    if self.list_pams[name_ref][name_asc][name_var] == 'U':
                        start_converted = self.starts[name_ref][name_asc][name_var]
                    if self.list_pams[name_ref][name_asc][name_var] == 'LU':
                        start_converted = np.log2(self.starts[name_ref][name_asc][name_var])
                    mc.starting_point[mc.variable_list[name_ref][name_asc][name_var]] = start_converted

    def convert(self, theta, name_ref, name_asc):
        dict_out = {}
        # If we need the parameters for the prior, we are not providing any name for the dataset
        for key in self.list_pams[name_ref][name_asc]:
                dict_out[key] = self.variables[name_ref][name_asc][key](
                    theta, self.fixed, self.var_list[name_ref][name_asc][key])
        return dict_out

    def return_priors(self, theta):
        prior_out = 0.00
        for name_ref in self.prior_pams:
            for name_asc in self.prior_pams[name_ref]:
                key_pams = self.convert(theta, name_ref, name_asc)
                for key in self.prior_pams[name_ref][name_asc]:
                    prior_out += giveback_priors(
                        self.prior_kind[name_ref][name_asc][key],
                        self.prior_pams[name_ref][name_asc][key],
                        key_pams[key])
        return prior_out

    def compute(self, theta, dataset):

        name_ref = dataset.name_ref
        output = np.zeros(dataset.n)
        for name_asc in self.order[name_ref]:
            dict_pams = self.convert(theta, name_ref, name_asc)
            coeff = np.zeros(self.order[name_ref][name_asc]+1)
            for var in self.list_pams[name_ref][name_asc]:
                coeff[self.order_ind[name_ref][name_asc][var]] = dict_pams[var]

            output += np.where(self.x_mask[name_ref][name_asc],
                               np.polynomial.polynomial.polyval(self.x_vals[name_ref][name_asc], coeff),
                               0.0)
        return output

    def initialize(self, mc):
        for name_ref in self.list_pams:
            for name_asc in self.list_pams[name_ref]:
                for name_var in self.list_pams[name_ref][name_asc]:
                    if name_var in mc.variable_list[name_ref][name_asc]:
                        mc.pam_names[mc.variable_list[name_ref][name_asc][name_var]] = name_var

    def print_vars(self, mc, theta):
        for name_ref in self.list_pams:
            for name_asc in self.list_pams[name_ref]:
                for name_var in self.list_pams[name_ref][name_asc]:
                    if name_var in mc.variable_list[name_ref][name_asc]:
                        var = self.variables[name_ref][name_asc][name_var](theta, self.fixed,
                                                                 self.var_list[name_ref][name_asc][name_var])
                    print 'Correlation ', name_ref, name_asc, name_var, var, '(', theta[
                            self.var_list[name_ref][name_asc][name_var]], ')'
        print

