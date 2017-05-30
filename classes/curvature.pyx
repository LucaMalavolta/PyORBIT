from common import *

class CurvatureCommonVariables:
    """Class to define the curvature of the dataset
    The zero-point of the curvature is always zero, since this parameter is already defined
    as the zero-point of each dataset """
    def __init__(self):

        self.bounds = {}
        self.starts = {}

        self.variables = {}
        self.var_list = {}
        self.fix_list = {}

        self.prior_kind = {}
        self.prior_pams = {}

        self.fix_list = {}
        self.fixed = []
        self.nfix = 0
        self.list_pams = {}

        self.order = 0
        self.order_ind = {}
        """order =0 means no trend
        RVtrend = a1*(t-Tref) + a2*(t-Tref)^2 + ... an*(t-Tref)^n """

        """These default bounds are used when the user does not define them in the yaml file"""
        self.default_bounds = {}

    def define_bounds(self, mc):

        mc.variable_list['Curvature'] = {}
        for n_ord in xrange(1, self.order+1):
            var = 'curvature_'+repr(n_ord)
            self.list_pams[var] = 'U'
            self.order_ind[var] = n_ord
            self.default_bounds[var] = np.asarray([-1000.0, 1000.0])

        for var in self.list_pams:
            if var in self.fix_list:
                self.variables[var] = get_fix_val
                self.var_list[var] = self.nfix
                self.fixed.append(self.fix_list[var])
                self.nfix += 1
            else:
                '''If no bounds have been specified in the input file, we use the default ones
                     Bounds must be provided in any case to avoid a failure of PyDE '''
                if var in self.bounds:
                    bounds_tmp = self.bounds[var]
                else:
                    bounds_tmp = self.default_bounds[var]

                if self.list_pams[var] == 'U':
                    self.variables[var] = get_var_val
                    mc.bounds_list.append(bounds_tmp)
                elif self.list_pams[var] == 'LU':
                    self.variables[var] = get_var_exp
                    mc.bounds_list.append(np.log2(bounds_tmp))

                self.var_list[var] = mc.ndim
                mc.variable_list['Curvature'][var] = mc.ndim
                mc.ndim += 1

    def starting_point(self, mc):

        """Default values are already set in the array"""

        for var in self.list_pams:
            if var in self.starts:
                if self.list_pams[var] == 'U':
                    start_converted = self.starts[var]
                elif self.list_pams[var] == 'LU':
                    start_converted = np.log2(self.starts[var])

                mc.starting_point[mc.variable_list['Curvature'][var]] = start_converted

    def return_priors(self, theta):
        prior_out = 0.00
        kep_pams = self.convert(theta)
        for key in self.prior_pams:
            prior_out += giveback_priors(self.prior_kind[key], self.prior_pams[key], kep_pams[key])
        return prior_out

    def convert(self, theta):
        dict_out = {}
        for key in self.list_pams:
            dict_out[key] = (self.variables[key](theta, self.fixed, self.var_list[key]))
        return dict_out

        # return self.variables[pl_name]['P'](theta, fixed, ), self.variables[pl_name]['K'](theta, fixed, 2), \
        #       self.variables[pl_name]['f'](theta, fixed, 2), self.variables[pl_name]['e'](theta, fixed, 2), \
        #       self.variables[pl_name]['o'](theta, fixed, 2)

    def compute(self, theta, dataset):
        dict_pams = self.convert(theta)
        coeff = np.zeros(self.order+1)
        for var in self.list_pams:
            coeff[self.order_ind[var]] = dict_pams[var]
        return np.polynomial.polynomial.polyval(dataset.x0, coeff)

    def model_curvature(self, dict_pams, x0):
        coeff = np.zeros(self.order+1)
        for var in self.list_pams:
            coeff[self.order_ind[var]] = dict_pams[var]
        return np.polynomial.polynomial.polyval(x0, coeff)

    def initialize(self, mc):
        for var in self.list_pams:
            mc.pam_names[mc.variable_list['Curvature'][var]] = var

    def print_vars(self, mc, theta):
        for var in self.list_pams:
            val = self.variables[var](theta, self.fixed, self.var_list[var])
            print 'Curvature ', var, val, self.var_list[var], '(', theta[self.var_list[var]], ')'
        print
