from common import *


class GaussianProcessAbstract:
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

        self.gp = {}

        self.model_type = 'gaussian'

        #self.list_pams_common = {
        #    'Prot': 'U',
        #    'Pdec': 'U',
        #    'Oamp': 'LU'}

        #self.list_pams_dataset = {'Hamp': 'U'}

        return

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

        for dataset_name, dataset in mc.dataset_dict.items():
            if self.model_name in dataset.models:
                for var in self.list_pams_dataset:
                    if var in self.fix_list[dataset_name]:
                        self.variables[dataset_name][var] = get_fix_val
                        self.var_list[dataset_name][var] = self.nfix
                        self.fixed.append(self.fix_list[dataset_name][var])
                        self.nfix += 1
                    else:
                        if self.list_pams_dataset[var] == 'U':
                            self.variables[dataset_name][var] = get_var_val
                            mc.bounds_list.append(self.bounds[dataset_name][var])
                        if self.list_pams_dataset[var] == 'LU':
                            self.variables[dataset_name][var] = get_var_exp
                            mc.bounds_list.append(np.log2(self.bounds[dataset_name][var]))
                        self.var_list[dataset_name][var] = mc.ndim
                        mc.variable_list[dataset_name][var] = mc.ndim
                        mc.ndim += 1

    def starting_point(self, mc):
        if 'Common' in self.starts:
            for var in self.starts['Common']:
                if self.list_pams_common[var] == 'U':
                    start_converted = self.starts['Common'][var]
                if self.list_pams_common[var] == 'LU':
                    start_converted = np.log2(self.starts['Common'][var])
                mc.starting_point[mc.variable_list['Common'][var]] = start_converted

        for dataset_name, dataset in mc.dataset_dict.items():
            if self.model_name in dataset.models and dataset_name in self.starts:
                    for var in self.starts[dataset_name]:
                        if self.list_pams_dataset[var] == 'U':
                            start_converted = self.starts[dataset_name][var]
                        if self.list_pams_dataset[var] == 'LU':
                            start_converted = np.log2(self.starts[dataset_name][var])
                        mc.starting_point[mc.variable_list[dataset_name][var]] = start_converted

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
        key_pams = self.convert(theta, d_name)
        if d_name is None:
            for key in self.prior_pams['Common']:
                prior_out += giveback_priors(self.prior_kind['Common'][key], self.prior_pams['Common'][key],
                                             key_pams[key])
        else:
            for key in self.prior_pams[d_name]:
                prior_out += giveback_priors(self.prior_kind[d_name][key], self.prior_pams[d_name][key], key_pams[key])
        return prior_out

    def lnlk_compute(self, theta, dataset):
        """ 2 steps:
           1) theta parameters must be converted in physical units (e.g. from logarithmic to linear space)
           2) physical values must be converted to {\tt george} input parameters
        """
        gp_pams = self.convert_val2gp(self.convert(theta, dataset.name_ref))

        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(dataset.x0, env)

        return self.gp[dataset.name_ref].log_likelihood(dataset.y - dataset.model, quiet=True)

    def sample_compute(self, theta, dataset):

        gp_pams = self.convert_val2gp(self.convert(theta, dataset.name_ref))

        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(dataset.x0, env)

        return self.gp[dataset.name_ref].sample_conditional(dataset.y - dataset.model, dataset.x0)

    def initialize(self, mc):

        for name in self.list_pams_common:
            if name in mc.variable_list['Common']:
                mc.pam_names[mc.variable_list['Common'][name]] = name

        for dataset_name, dataset in mc.dataset_dict.items():
            for name in self.list_pams_dataset:
                if name in mc.variable_list[dataset_name]:
                    mc.pam_names[mc.variable_list[dataset_name][name]] = name

            if self.model_name in dataset.models:
                self.define_kernel(dataset)


    def print_vars(self, mc, theta):

        for name in self.list_pams_common:
            if name in mc.variable_list['Common']:
                var = self.variables['Common'][name](theta, self.fixed, self.var_list['Common'][name])
                print 'GaussianProcess ', name, var, self.var_list['Common'][name], '(', theta[
                        self.var_list['Common'][name]], ')'

        #for dataset in mc.dataset_dict.itervalues():
        for dataset_name in mc.dataset_dict.iterkeys():
            for name in self.list_pams_dataset:
                if name in mc.variable_list[dataset_name]:
                    var = self.variables[dataset_name][name](theta, self.fixed,
                                                               self.var_list[dataset_name][name])
                    print 'GaussianProcess ', dataset_name, name, var, '(', theta[
                            self.var_list[dataset_name][name]], ')'
        print


class GaussianProcess_QuasiPeriodicActivity(GaussianProcessAbstract):
    ''' Three parameters out of four are the same for all the datasets, since they are related to
    the properties of the physical process rather than the observed effects on a dataset
     From Grunblatt+2015, Affer+2016
     - theta: is usually related to the rotation period of the star( or one of its harmonics);
     - lambda: is the correlation decay timescale, and it can be related to the lifetime of the active regions.
     - omega: is the length scale of the periodic component, and can be linked to the size evolution of the active regions;
     - h: represents the amplitude of the correlations '''

    model_name = 'gaussian'
    common_ref = 'Common'

    list_pams_common = {
        'Prot': 'U',
        'Pdec': 'U',
        'Oamp': 'LU'}
    list_pams_dataset = {'Hamp': 'U'}

    n_pams = 4

    """ Indexing is determined by the way the kernel is constructed"""
    gp_pams_index = {
        'Hamp': 0, # amp2
        'Pdec': 1, # metric
        'Oamp': 2, # gamma
        'Prot': 3 # ln_P
    }

    # AAAAHHHHHH!!!!!
    #
    # self.kernel = gp_pams['amp2'] * george.kernels.ExpSquaredKernel(metric=gp_pams['metric']) * \
    #         george.kernels.ExpSine2Kernel(gamma=gp_pams['gamma'], log_period=gp_pams['ln_P'])

    def convert_val2gp(self, input_pams):
        """
        :param input_pam: dictonary with the 'physically meaningful' parameters of the GP kernel
        :return: dictonary with the parameters to be fed to 'george'
        WARNING: this subroutine is HIGHLY specific of your choice of the kernel! I reccomend to
        create a new Class with different transformations if you are planning of using a different
        kernel combination
        """
        output_pams = np.zeros(self.n_pams, dtype=np.double)

        """ You must check _george_ documentation (and possibily do a lot of testing) to know how to convert physical 
        values to the parameter vector accepted by george.set_parameter_vector() function. Note: these values may be 
        different from ones accepted by the kernel
        """
        output_pams[self.gp_pams_index['Hamp']] = np.log(input_pams['Hamp'])*2
        output_pams[self.gp_pams_index['Pdec']] = np.log(input_pams['Pdec'])*2
        output_pams[self.gp_pams_index['Oamp']] = 1. / (2.*input_pams['Oamp'] ** 2)
        output_pams[self.gp_pams_index['Prot']] = np.log(input_pams['Prot'])

        return output_pams

    def convert_gp2val(self, input_pams):
        """
        :param input_pam: dictonary with the parameters to be fed to 'george'
        :return: dictonary with the 'physically meaningful' parameters of the GP kernel
        WARNING: this subroutine is HIGHLY specific of your choice of the kernel! I reccomend to
        create a new Class with different transformations if you are planning of using a different
        kernel combination
        """
        return {
            'Hamp': np.exp(input_pams[self.gp_pams_index['Hamp']]/2.0),
            'Pdec': np.exp(input_pams[self.gp_pams_index['Pdec']] / 2.0),
            'Oamp': np.sqrt(1. / (2.*input_pams[self.gp_pams_index['Oamp']])),
            'Prot': np.exp(input_pams[self.gp_pams_index['Prot']])
        }

    def define_kernel(self, dataset):

        gp_pams = np.ones(self.n_pams)
        """ Kernel initialized with fake values... don't worry, they'll be overwritten soon"""
        self.kernel = np.exp(gp_pams[0]) * \
                      george.kernels.ExpSquaredKernel(metric=np.exp(gp_pams[1])) * \
                      george.kernels.ExpSine2Kernel(gamma=gp_pams[1], log_period=gp_pams[2])



        """
         gp_pams[0] = h^2 -> h^2 * ExpSquaredKernel * ExpSine2Kernel
           -> set_parameter_vector() accepts the natural logarithm of this value
         gp_pams[1] = metric = r^2 = lambda**2  -> ExpSquaredKernel(metric=r^2)
           -> set_parameter_vector() accepts the natural logarithm of this value
         gp_pams[2] = Gamma =  1/ (2 omega**2) -> ExpSine2Kernel(gamma, ln_period)
         gp_pams[3] = ln_theta = ln_Period -> ExpSine2Kernel(gamma, ln_period)
         
        """

        self.gp[dataset.name_ref] = george.GP(self.kernel)

        """ I've decided to add the jitter in quadrature instead of using a constant kernel to allow the use of 
        different / selective jitter within the dataset
        """
        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)

        self.gp[dataset.name_ref].compute(dataset.x0, env)

        return