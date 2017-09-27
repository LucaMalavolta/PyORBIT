from common import *

class ModelContainer:
    def __init__(self):
        self.planet_dict = {}
        self.dynamical_dict = {}
        self.dynamical_t0_dict = {}
        self.dynamical_model = None

        self.dataset_dict = {}
        self.dataset_index = {}
        self.n_datasets = 0

        self.models = {}
        self.common_models = {}

        """ pyde/emcee variables
        wondering if I should move these somewhere else, for now they are staying here because they are 
        essentially harmless """
        self.pyde_dir_output = None
        self.emcee_dir_output = None
        self.polychord_dir_output = None

        self.emcee_parameters = {'nsave': 0, 'npop_mult': 2, 'thin': 1, 'nburn':0,
                                 'multirun': None, 'multirun_iter': 20}
        self.pyde_parameters = {'ngen': 1000, 'npop_mult': 2}

        # Default values, taken from the PyPolyChord wrapper in PolyChord official distribution, V1.9
        self.polychord_parameters = {'nlive_mult': 25,
                                     'num_repeats_mult': 5,
                                     'feedback': 1,
                                     'precision_criterion': 0.001,
                                     'max_ndead': -1,
                                     'boost_posterior': 0.0,
                                     'read_resume': True,
                                     'base_dir': 'polychord/',
                                     'shutdown_jitter': False}

        self.ndata = None
        self.ndof = None

        self.starting_point = None
        self.starting_point_flag = False
        self.recenter_bounds_flag = True

        self.bounds = None
        self.range = None
        self.ndim = 0
        self.pam_names = ''
        self.star_mass = [1.0000, 0.1000]
        self.star_radius = [1.0000, 0.1000]

        self.Tref = None

        """
        Values have been taken from TRADES
        These variables will be renamed in the next release, right now I'm keeping the original names
        to avoid breaking the code
        """
        self.G_grav = constants.Gsi  # Gravitational Constants in SI system [m^3/kg/s^2]
        self.G_ttvfast = constants.Giau  # G [AU^3/Msun/d^2]
        self.M_SJratio = constants.Msjup
        self.M_SEratio = constants.Msear
        self.M_JEratio = constants.Mjear

        self.R_SJratio = constants.Rsjup
        self.R_JEratio = constants.Rjear
        self.R_SEratio = constants.Rsjup * constants.Rjear

        self.Mu_sun = constants.Gsi * constants.Msun
        self.seconds_in_day = constants.d2s
        self.AU_km = constants.AU
        self.AUday2ms = self.AU_km / self.seconds_in_day * 1000.0

    def model_setup(self):
        self.n_datasets = len(self.dataset_dict)
        for dataset in self.dataset_dict.itervalues():

            if 'sinusoids' in dataset.models:
                self.scv.model_setup(dataset)

            dataset.model_reset()

    def initialize_model(self):

        # First step: setting up the models
        # Second step: define the number of variables and setting up the boundaries
        # Third step: initialize the variables names

        self.model_setup()
        self.create_variables_bounds()

        self.ndata = 0
        for dataset in self.dataset_dict.itervalues():
            if not dataset.models:
                continue
            self.ndata += dataset.n
        self.ndof = self.ndata - self.ndim

        #self.initialize()

    def create_variables_bounds(self):
        # This routine creates the boundary array and at the same time
        # creates a dictionary with the name of the arrays and their
        # positions in bounds/theta array so that they can be accessed
        # without using nested counters

        self.ndim = 0
        bounds_list = []
        for model in self.models.itervalues():
            self.ndim, bounds_ext = self.common_models[model.common_ref].define_variables_bounds(
                    self.ndim, model.list_pams_common)
            bounds_list.extend(bounds_ext)

        for dataset in self.dataset_dict.itervalues():
            self.ndim, bounds_ext = dataset.define_variables_bounds(self.ndim, dataset.list_pams)
            bounds_list.extend(bounds_ext)

            for model_name in dataset.models:
                self.ndim, bounds_ext = self.models[model_name].define_variables_bounds(self.ndim, dataset.name_ref)
                bounds_list.extend(bounds_ext)

        self.bounds = np.asarray(bounds_list)
        self.range = self.bounds[:, 1] - self.bounds[:, 0]

    """ 
    def initialize(self):
        # Function with two goals:
        # * Unfold and print out the output from theta
        # * give back a parameter name associated to each value in the result array

        self.pam_names = self.ndim * ['']

        for dataset in self.dataset_dict.itervalues():
            dataset.initialize(self)

        for model in self.models.itervalues():
            model.initialize(self)

    """
    def create_starting_point(self):

        self.starting_point = np.average(self.bounds, axis=1)

        for model in self.common_models.itervalues():
            model.define_starting_point(self.starting_point)

        for dataset in self.dataset_dict.itervalues():
            dataset.define_starting_point(self.starting_point)

            for model in dataset.models:
                self.models[model].define_starting_point(self.starting_point)

    def check_bounds(self, theta):
        for ii in xrange(0, self.ndim):
            if not (self.bounds[ii, 0] < theta[ii] < self.bounds[ii, 1]):
                return False

        period_storage = []
        for planet_name in self.planet_dict:

            """ Step 1: save the all planet periods into a list"""
            period_storage.extend(
                [self.common_models[planet_name].transformation['P'](theta, self.common_models[planet_name].fixed,
                                                               self.common_models[planet_name].variable_index['P'])])

            """ Step 2: check if the eccentricity is within the given range"""
            e = self.common_models[planet_name].transformation['e'](theta,
                                                              self.common_models[planet_name].fixed,
                                                              self.common_models[planet_name].variable_index['e'])
            if not self.common_models[planet_name].bounds['e'][0] <= e < self.common_models[planet_name].bounds['e'][1]:
                return False

        """ Step 4 check ofr overlapping periods (within 5% arbitrarily chosen)"""
        for i_n, i_v in enumerate(period_storage):
            if i_n == len(period_storage) - 1: break
            if np.amin(np.abs(period_storage[i_n + 1:] - i_v)) / i_v < 0.05:
                return False

        return True

    def __call__(self, theta):
        if not self.check_bounds(theta):
            return -np.inf

        if self.dynamical_model is not None:
            """ check if any keyword ahas get the output model from the dynamical tool
            we must do it here because all the planet are involved"""
            dynamical_output = self.dynamical_model.compute(self, theta)

        logchi2_out = 2. * self.ndof * np.log(2 * np.pi)

        for model in self.common_models.itervalues():
            logchi2_out += model.return_priors(theta)

        for dataset_name, dataset in self.dataset_dict.iteritems():
            dataset.model_reset()
            variable_values = dataset.convert(theta)
            dataset.compute(variable_values)

            if 'none' in dataset.models or 'None' in dataset.models:
                continue
            if not dataset.models:
                continue

            logchi2_gp_model = None
            for model_name in dataset.models:

                logchi2_out += self.models[model_name].return_priors(theta, dataset_name)

                if hasattr(self.models[model_name], 'internal_likelihood'):
                    logchi2_gp_model = model_name.copy()
                    continue
                #if self.models[model_name].model_class == 'gaussian_process':

                if dataset.dynamical:
                    dataset.model += dynamical_output[dataset_name]
                    continue

                common_ref = self.models[model_name].common_ref
                variable_values = self.common_models[common_ref].convert(theta)
                variable_values.update(self.models[model_name].convert(theta, dataset_name))
                dataset.model += self.models[model_name].compute(variable_values, dataset)

            """ Gaussian Process check MUST be the last one or the program will fail
             that's because for the GP to work we need to know the _deterministic_ part of the model 
             (i.e. the theoretical values you get when you feed your model with the parameter values) """
            if logchi2_gp_model:
                common_ref = self.models[logchi2_gp_model].common_ref
                variable_values = self.common_models[common_ref].convert(theta)
                variable_values.update(self.models[logchi2_gp_model].convert(theta, dataset_name))
                logchi2_out += self.models[logchi2_gp_model].lnlk_compute(variable_values, dataset)
            else:
                logchi2_out += dataset.model_logchi2()

        return logchi2_out

    def recenter_bounds(self, pop_mean, recenter=True):
        # This function recenters the bounds limits for circular variables
        # Also, it extends the range of a variable if the output of PyDE is a fixed number

        ind_list = []

        for model in self.common_models.itervalues():
            ind_list.extend(model.special_index_recenter_bounds())
            ind_list.extend(model.index_recenter_bounds())

        for dataset in self.dataset_dict.itervalues():
            for model in dataset.models:
                ind_list.extend(self.models[model].special_index_recenter_bounds(dataset.name_ref))
                ind_list.extend(self.models[model].index_recenter_bounds(dataset.name_ref))

        if not recenter:
            return ind_list

        if ind_list:
            tmp_range = (self.bounds[:, 1] - self.bounds[:, 0]) / 2
            replace_bounds = np.zeros([self.ndim, 2])
            replace_bounds[:, 0] = pop_mean - tmp_range
            replace_bounds[:, 1] = pop_mean + tmp_range
            self.bounds[ind_list, :] =  replace_bounds[ind_list, :]

    def fix_population(self, pop_mean, population):

        ind_list = self.recenter_bounds(pop_mean, recenter=False)
        n_pop = np.size(population, axis=0)

        if ind_list:
            for var_ind in ind_list:
                fix_sel = (population[:, var_ind] <= self.bounds[var_ind, 0]) | (
                    population[:, var_ind] >= self.bounds[var_ind, 1])
                population[fix_sel, var_ind] = pop_mean[var_ind]

        for ii in xrange(0, self.ndim):
            if np.amax(population[:, ii]) - np.amin(population[:, ii]) < 10e-7:
                range_restricted = (self.bounds[ii, 1] - self.bounds[ii, 0]) / 100.
                min_bound = np.maximum((pop_mean[ii] - range_restricted / 2.0), self.bounds[ii, 0])
                max_bound = np.minimum((pop_mean[ii] + range_restricted / 2.0), self.bounds[ii, 1])
                population[:, ii] = np.random.uniform(min_bound, max_bound, n_pop)

        return population

    def results_resumen(self, theta, skip_theta=False):
        # Function with two goals:
        # * Unfold and print out the output from theta
        # * give back a parameter name associated to each value in the result array

        print
        print '===================================================================================================='
        print '     ------------------------------------------------------------------------------------------     '
        print '===================================================================================================='
        print
        for dataset_name, dataset in self.dataset_dict.iteritems():
            print '----- dataset: ', dataset_name
            print_theta_bounds(dataset.variable_sampler, theta, self.bounds, skip_theta)

            for model_name in dataset.models:
                print '---------- ', dataset_name, '     ----- model: ', model_name
                print_theta_bounds(self.models[model_name].variable_sampler[dataset_name], theta, self.bounds, skip_theta)

        for model in self.common_models.itervalues():
            print '----- common model: ', model.common_ref
            print_theta_bounds(model.variable_sampler, theta, self.bounds, skip_theta)

        if skip_theta:
            return

        print '===================================================================================================='
        print '===================================================================================================='
        print

        for dataset_name, dataset in self.dataset_dict.iteritems():
            print '----- dataset: ', dataset_name
            variable_values = dataset.convert(theta)
            print_dictionary(variable_values)

            print
            for model_name in dataset.models:
                print '---------- ', dataset_name, '     ----- model: ', model_name
                variable_values = self.models[model_name].convert(theta, dataset_name)
                print_dictionary(variable_values)

        for model in self.common_models.itervalues():
            print '----- common model: ', model.common_ref
            variable_values = model.convert(theta)
            print_dictionary(variable_values)

        print '===================================================================================================='
        print '     ------------------------------------------------------------------------------------------     '
        print '===================================================================================================='
        print

    def get_theta_dictionary(self):
        # * give back a parameter name associated to each value in the result array

        theta_dictionary = {}
        for dataset_name, dataset in self.dataset_dict.iteritems():
            for var, i in dataset.variable_sampler.iteritems():
                try:
                    theta_dictionary[dataset_name + '_' + var] = i
                except:
                    theta_dictionary[repr(dataset_name) + '_' + var] = i

            for model_name in dataset.models:
                for var, i in self.models[model_name].variable_sampler[dataset_name].iteritems():
                    try:
                        theta_dictionary[dataset_name + '_' + model_name +  '_' + var] = i
                    except:
                        theta_dictionary[repr(dataset_name) + '_' + model_name + '_' + var] = i

        for model in self.common_models.itervalues():
            for var, i in model.variable_sampler.iteritems():
                theta_dictionary[model.common_ref + '_' + var] = i

        return theta_dictionary

    def get_model(self, theta):
        model_out = {}
        if self.dynamical_model is not None:
            """ check if any keyword ahas get the output model from the dynamical tool
            we must do it here because all the planet are involved"""
            dynamical_output = self.dynamical_model.compute(self, theta)

        logchi2_out = 2. * self.ndof * np.log(2 * np.pi)

        for model in self.common_models.itervalues():
            logchi2_out += model.return_priors(theta)

        for dataset_name, dataset in self.dataset_dict.iteritems():
            model_out[dataset_name] = {}
            dataset.model_reset()

            if hasattr(self, 'deepcopy_for_plot'):
                model_out[dataset_name]['systematics'] = np.zeros(dataset.n, dtype=np.double)
                model_out[dataset_name]['jitter'] = np.zeros(dataset.n, dtype=np.double)
                model_out[dataset_name]['complete'] = np.zeros(dataset.n, dtype=np.double)
            else:
                variable_values = dataset.convert(theta)
                dataset.compute(variable_values)
                model_out[dataset_name]['systematics'] = dataset.model.copy()
                model_out[dataset_name]['jitter'] = dataset.jitter.copy()
                model_out[dataset_name]['complete'] = dataset.model.copy()


            if 'none' in dataset.models or 'None' in dataset.models:
                continue
            if not dataset.models:
                continue

            logchi2_gp_model = None

            for model_name in dataset.models:

                logchi2_out += self.models[model_name].return_priors(theta, dataset_name)

                if hasattr(self.models[model_name], 'internal_likelihood'):
                    logchi2_gp_model = model_name.copy()
                    continue

                if dataset.dynamical:
                    dataset.model += dynamical_output[dataset_name]

                    model_out[dataset_name]['complete'] += dataset.model
                    model_out[dataset_name][model_name] = dynamical_output[dataset_name].copy()
                    continue

                common_ref = self.models[model_name].common_ref
                variable_values = self.common_models[common_ref].convert(theta)
                variable_values.update(self.models[model_name].convert(theta, dataset_name))

                model_out[dataset_name][model_name] = self.models[model_name].compute(variable_values, dataset)

                dataset.model += model_out[dataset_name][model_name]
                model_out[dataset_name]['complete'] += model_out[dataset_name][model_name]

            """ Gaussian Process check MUST be the last one or the program will fail
             that's because for the GP to work we need to know the _deterministic_ part of the model 
             (i.e. the theoretical values you get when you feed your model with the parameter values) """
            if logchi2_gp_model:
                common_ref = self.models[logchi2_gp_model].common_ref
                variable_values = self.common_models[common_ref].convert(theta)
                variable_values.update(self.models[logchi2_gp_model].convert(theta, dataset))
                logchi2_out += self.models[logchi2_gp_model].lnlk_compute(variable_values, dataset)

                model_out[dataset_name][model_name] = \
                    self.models[logchi2_gp_model].sample_conditional(variable_values, dataset)
                model_out[dataset_name]['complete'] += model_out[dataset_name][model_name]
            else:
                logchi2_out += dataset.model_logchi2()

        # workaround to avoid memory leaks from GP module
        #gc.collect()


        return model_out, logchi2_out


def print_theta_bounds(i_dict, theta, bounds, skip_theta=False):
    format_string = "%10s  %4d  %12f ([%10f, %10f])"
    format_string_notheta = "%10s  %4d  ([%10f, %10f])"
    format_string_long = "%10s  %4d  %12f   %12f %12f (15-84 p) ([%9f, %9f])"

    for var, i in i_dict.iteritems():
        if skip_theta:
            print format_string_notheta % (var, i, bounds[i, 0], bounds[i, 1])
        elif len(np.shape(theta)) == 2:
            perc0, perc1, perc2 = np.percentile(theta[:, i], [15.865, 50, 84.135], axis=0)
            print format_string_long %(var, i, perc1, perc0-perc1, perc2-perc1, bounds[i, 0], bounds[i, 1])
        else:
            print format_string % (var, i, theta[i], bounds[i, 0], bounds[i, 1])
    print


def print_dictionary(variable_values):
    format_string_long = "%10s   %15f   %15f %15f (15-84 p)"
    format_string = "%10s   %15f "
    for var_names, var_vals in variable_values.iteritems():
        if np.size(var_vals) > 1:
            perc0, perc1, perc2 = np.percentile(var_vals, [15.865, 50, 84.135], axis=0)
            print format_string_long %(var_names,  perc1, perc0-perc1, perc2-perc1)
        else:
            print format_string % (var_names, var_vals)
    print


class ModelContainerPolyChord(ModelContainer):
    def polychord_priors(self, cube):
        theta = (self.bounds[:, 1] - self.bounds[:, 0]) * cube + self.bounds[:, 0]
        return theta.tolist()

    def polychord_call(self, theta1):
        theta = np.empty(self.ndim)
        for i in xrange(0, self.ndim):
            theta[i] = theta1[i]
        phi = [0.0] * 0
        chi_out = self(theta)
        if chi_out < -0.5e10:
            return -0.5e10, phi
        return chi_out, phi


class ModelContainerMultiNest(ModelContainer):
    def multinest_priors(self, cube, ndim, nparams):
        # cube[:] = (self.bounds[:, 1] - self.bounds[:, 0]) * cube[:] + self.bounds[:, 0]
        for i in xrange(0, ndim):
            cube[i] = (self.bounds[i, 1] - self.bounds[i, 0]) * cube[i] + self.bounds[i, 0]

    def multinest_call(self, theta1, ndim, nparams):
        # Workaround for variable selection: if a variable as null index
        # (i.e. it has not been included in the model)
        # the numpy array will give back an empty list, the ctype will give back an error
        theta = np.empty(ndim)
        for i in xrange(0, ndim):
            theta[i] = theta1[i]
        chi_out = self(theta)
        if chi_out < -0.5e10:
            return -0.5e10
        return chi_out


class ModelContainerDnest4(ModelContainer):
    def from_prior(self):
        return self.starting_point

    def perturb(self, params):
        """
        Unlike in C++, this takes a numpy array of parameters as input,
        and modifies it in-place. The return value is still logH.
        """
        logH = 0.0
        which = np.random.randint(np.size(params))

        logH -= -0.5 * (params[which] / self.range[which]) ** 2
        params_mod = params[which] + self.range[which] * dnest4.randh()
        params[which] = dnest4.wrap(params_mod, self.bounds[which, 0], self.bounds[which, 1])

        logH += -0.5 * (params[which] / self.range[which]) ** 2

        # if which == 0:
        #    print which,  params[0], np.exp2(params[0])

        return logH

    def log_likelihood(self, theta):
        return self(theta)
