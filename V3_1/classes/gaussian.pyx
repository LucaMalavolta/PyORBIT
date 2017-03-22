from common import *

class GaussianProcessCommonVariables:
    def __init__(self):
        self.fix_list = {}
        self.bounds = {}
        self.starts = {}

        self.var_list = {}
        self.variables = {}

        self.prior_kind = {}
        self.prior_pams = {}

        self.fixed = []
        self.nfix = 0

        ''' Three parameters out of four are the same for all the datasets, since they are related to
        the properties of the physical process rather than the observed effects on a dataset
         From Grunblatt+2015, Affer+2016
         - theta: is usually related to the rotation period of the star( or one of its harmonics);
         - lambda: is the correlation decay timescale, and it can be related to the lifetime of the active regions.
         - omega: is the length scale of the periodic component, and can be linked to the size evolution of the active regions;
         - h: represents the amplitude of the correlations '''

        self.list_pams_common = {
            'Prot': 'U',
            'Pdec': 'U',
            'Oamp': 'LU'}

        self.list_pams_dataset = {'Hamp': 'U'}

        return

    def convert_val2gp(self, input_pam):
        """
        :param input_pam: dictonary with the 'physically meaningful' parameters of the GP kernel
        :return: dictonary with the parameters to be fed to 'george'
        WARNING: this subroutine is HIGHLY specific of your choice of the kernel! I reccomend to
        create a new Class with different transformations if you are planning of using a different
        kernel combination
        """
        return {
            'ln_P': np.log(input_pam['Prot']),
            'metric': input_pam['Pdec'] ** 2,
            'gamma': 2. / (input_pam['Oamp'] ** 2),
            'amp2': np.power(input_pam['Hamp'], 2)}

    def convert_gp2val(self, input_pam):
        """
        :param input_pam: dictonary with the parameters to be fed to 'george'
        :return: dictonary with the 'physically meaningful' parameters of the GP kernel
        WARNING: this subroutine is HIGHLY specific of your choice of the kernel! I reccomend to
        create a new Class with different transformations if you are planning of using a different
        kernel combination
        """
        return {
            'Prot': np.exp(input_pam['ln_P']),
            'Pdec': np.sqrt(input_pam['metric']),
            'Oamp': np.sqrt(2. / input_pam['gamma']),
            'Hamp': np.sqrt(input_pam['amp2'])}

    def add_dataset(self, name_ref='Common'):
        # planet_name : name of the dataset, otherwise 'common'

        self.bounds[name_ref] = {}
        self.variables[name_ref] = {}
        self.fix_list[name_ref] = {}
        self.var_list[name_ref] = {}

        self.prior_kind[name_ref] = {}
        self.prior_pams[name_ref] = {}

        if name_ref == 'Common':
            for name in self.list_pams_common:
                self.bounds['Common'][name] = [0.01, 100]
        else:
            for name in self.list_pams_dataset:
                self.bounds[name_ref][name] = [0.00001, 100]

    def define_bounds(self, mc):
        for var in self.list_pams_common:
            if var in self.fix_list['Common']:
                self.variables['Common'][var] = get_fix_val
                self.var_list['Common'][var] = self.nfix
                self.fixed.append(self.fix_list['Common'][var])
                self.nfix += 1
            else:
                if self.list_pams_common[var] == 'U':
                    self.variables['Common'][var] = get_var_val
                    mc.bounds_list.append(self.bounds['Common'][var])
                if self.list_pams_common[var] == 'LU':
                    self.variables['Common'][var] = get_var_exp
                    mc.bounds_list.append(np.log2(self.bounds['Common'][var]))

                self.var_list['Common'][var] = mc.ndim
                mc.variable_list['Common'][var] = mc.ndim
                mc.ndim += 1

        for dataset in mc.dataset_list:
            if 'gaussian' in dataset.models:
                for var in self.list_pams_dataset:
                    if var in self.fix_list[dataset.name_ref]:
                        self.variables[dataset.name_ref][var] = get_fix_val
                        self.var_list[dataset.name_ref][var] = self.nfix
                        self.fixed.append(self.fix_list[dataset.name_ref][var])
                        self.nfix += 1
                    else:
                        if self.list_pams_dataset[var] == 'U':
                            self.variables[dataset.name_ref][var] = get_var_val
                            mc.bounds_list.append(self.bounds[dataset.name_ref][var])
                        if self.list_pams_dataset[var] == 'LU':
                            self.variables[dataset.name_ref][var] = get_var_exp
                            mc.bounds_list.append(np.log2(self.bounds[dataset.name_ref][var]))
                        self.var_list[dataset.name_ref][var] = mc.ndim
                        mc.variable_list[dataset.name_ref][var] = mc.ndim
                        mc.ndim += 1

    def starting_point(self, mc):
        if 'Common' in self.starts:
            for var in self.starts['Common']:
                if self.list_pams_common[var] == 'U':
                    start_converted = self.starts['Common'][var]
                if self.list_pams_common[var] == 'LU':
                    start_converted = np.log2(self.starts['Common'][var])
                mc.starting_point[mc.variable_list['Common'][var]] = start_converted

        for dataset in mc.dataset_list:
            if 'gaussian' in dataset.models and dataset.name_ref in self.starts:
                    for var in self.starts[dataset.name_ref]:
                        if self.list_pams_common[var] == 'U':
                            start_converted = self.starts[dataset.name_ref][var]
                        if self.list_pams_common[var] == 'LU':
                            start_converted = np.log2(self.starts[dataset.name_ref][var])
                        mc.starting_point[mc.variable_list[dataset.name_ref][var]] = start_converted

    def convert(self, theta, d_name=None):
        dict_out = {}
        for key in self.list_pams_common:
            dict_out[key] = self.variables['Common'][key](theta, self.fixed, self.var_list['Common'][key])
        # If we need the parameters for the prior, we are not providing any name for the dataset
        if d_name is not None:
            for key in self.list_pams_dataset:
                dict_out[key] = self.variables[d_name][key](theta, self.fixed, self.var_list[d_name][key])
        return dict_out

    def return_priors(self, theta, d_name=None):
        prior_out = 0.00
        kep_pams = self.convert(theta, d_name)
        if d_name is None:
            for key in self.prior_pams['Common']:
                prior_out += giveback_priors(self.prior_kind['Common'][key], self.prior_pams['Common'][key],
                                             kep_pams[key])
        else:
            for key in self.prior_pams[d_name]:
                prior_out += giveback_priors(self.prior_kind[d_name][key], self.prior_pams[d_name][key], kep_pams[key])
        return prior_out

    def lnlk_compute(self, theta, dataset):
        # 2 steps:
        #   1) theta parameters must be converted in physical units (e.g. from logarithmic to linear space)
        #   2) physical values must be converted to {\tt george} input parameters
        gp_pams = self.convert_val2gp(self.convert(theta, dataset.name_ref))
        # gp_pams['ln_P] = ln_theta = ln_Period -> ExpSine2Kernel(gamma, ln_period)
        # gp_pams['metric'] = metric = r^2 = lambda**2  -> ExpSquaredKernel(metric=r^2)
        # gp_pams['gamma'] = Gamma =  1/ (2 omega**2) -> ExpSine2Kernel(gamma, ln_period)
        # gp_pams['amp2] = h^2 -> h^2 * ExpSquaredKernel * ExpSine2Kernel
        kernel = gp_pams['amp2'] * george.kernels.ExpSquaredKernel(metric=gp_pams['metric']) * \
                 george.kernels.ExpSine2Kernel(gamma=gp_pams['gamma'], ln_period=gp_pams['ln_P'])

        gp = george.GP(kernel)
        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        gp.compute(dataset.x0, env)

        return gp.lnlikelihood(dataset.y - dataset.model, quiet=True)

    def sample_compute(self, theta, dataset):

        gp_pams = self.convert_val2gp(self.convert(theta, dataset.name_ref))

        kernel = gp_pams['amp2'] * george.kernels.ExpSquaredKernel(metric=gp_pams['metric']) * \
                 george.kernels.ExpSine2Kernel(gamma=gp_pams['gamma'], ln_period=gp_pams['ln_P'])

        gp = george.GP(kernel)
        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        gp.compute(dataset.x0, env)

        return gp.sample_conditional(dataset.y - dataset.model, dataset.x)

    def print_vars(self, mc, theta):

        for name in self.list_pams_common:
            if name in mc.variable_list['Common']:
                mc.pam_names[mc.variable_list['Common'][name]] = name
                var = self.variables['Common'][name](theta, self.fixed, self.var_list['Common'][name])
                print 'GaussianProcess ', name, var, self.var_list['Common'][name], '(', theta[
                    self.var_list['Common'][name]], ')'

        for dataset in mc.dataset_list:
            for name in self.list_pams_dataset:
                if name in mc.variable_list[dataset.name_ref]:
                    mc.pam_names[mc.variable_list[dataset.name_ref][name]] = name
                    var = self.variables[dataset.name_ref][name](theta, self.fixed,
                                                                 self.var_list[dataset.name_ref][name])
                    print 'GaussianProcess ', dataset.name_ref, name, var, '(', theta[
                        self.var_list[dataset.name_ref][name]], ')'
