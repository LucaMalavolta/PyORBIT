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

                        if var in self.list_pams:
                            if self.list_pams[var] == 'U':
                                self.variables[dataset.name_ref][var] = get_var_val
                                mc.bounds_list.append(bounds_tmp)
                            if self.list_pams[var] == 'LU':
                                self.variables[dataset.name_ref][var] = get_var_exp
                                mc.bounds_list.append(np.log2(bounds_tmp))
                        else:
                            self.list_pams[var] == 'U'
                            self.variables[dataset.name_ref][var] = get_var_val
                            mc.bounds_list.append(bounds_tmp)

                        self.var_list[dataset.name_ref][var] = mc.ndim
                        mc.variable_list[dataset.name_ref][var] = mc.ndim
                        mc.ndim += 1



