from common import *

__all__ = ["ModelContainer"]


class ModelContainer:
    def __init__(self):
        self.planet_dict = {}
        self.dynamical_dict = {}
        self.dynamical_t0_dict = {}
        self.dynamical_model = None

        self.dataset_dict = {}

        self.models = {}
        self.common_models = {}

        """ pyde/emcee variables
        wondering if I should move these somewhere else, for now they are staying here because they are 
        essentially harmless """
        self.pyde_dir_output = None
        self.emcee_dir_output = None

        self.emcee_parameters = {'nsave': 0, 'npop_mult': 2, 'thin': 1, 'nburn':0,
                                 'multirun': None, 'multirun_iter': 20, 'version': '2.2.1'}

        self.pyde_parameters = {'ngen': 1000, 'npop_mult': 2}

        self.polychord_parameters = {'nlive_mult': 25,
                                'num_repeats_mult': 5,
                                'feedback': 1,
                                'precision_criterion': 0.001,
                                'sampling_efficiency': 0.8,
                                'max_ndead': -1,
                                'boost_posterior': 0.0,
                                'read_resume': True,
                                'base_dir': './',
                                'shutdown_jitter': False}
        self.polychord_dir_output = None

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
        # First step: setting up the correct associations between models and dataset

        for model_name, model in self.models.iteritems():
            if not model.model_conf:
                continue

            model.initialize_model(self, **model.model_conf)

            for dataset_name in list(set(model.model_conf) & set(self.dataset_dict)):
                model.setup_dataset(self.dataset_dict[dataset_name], **model.model_conf)

        for dataset_name, dataset in self.dataset_dict.iteritems():
            for model_name in dataset.models:
                self.models[model_name].common_initialization_with_dataset(dataset)

                common_model = self.models[model_name].common_ref
                if common_model:
                    self.common_models[common_model].common_initialization_with_dataset(dataset)

                #print self.models[model_name].model_conf['common']
                #for common_model in self.models[model_name].model_conf['common']:
                #    if hasattr(self.common_models[common_model], 'common_initialization_with_dataset'):
                #        self.common_models[common_model].setup_dataset(self.dataset_dict[dataset_name])


        if self.dynamical_model:
            self.dynamical_model.prepare(self)

    def create_variables_bounds(self):
        # This routine creates the boundary array and at the same time
        # creates a dictionary with the name of the arrays and their
        # positions in bounds/theta array so that they can be accessed
        # without using nested counters

        self.ndim = 0
        bounds_list = []
        for model in self.models.itervalues():
            if model.common_ref:
                model.default_bounds = self.common_models[model.common_ref].default_bounds
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

    def initialize_logchi2(self):

        # Second step: define the number of variables and setting up the boundaries
        # To be done only the first time ever
        self.ndata = 0
        for dataset in self.dataset_dict.itervalues():
            if not dataset.models:
                continue
            self.ndata += dataset.n
        self.ndof = self.ndata - self.ndim

    def create_starting_point(self):

        self.starting_point = np.average(self.bounds, axis=1)

        for model in self.common_models.itervalues():
            model.define_starting_point(self.starting_point)

        for dataset_name, dataset in self.dataset_dict.iteritems():
            dataset.define_starting_point(self.starting_point)

            for model in dataset.models:
                self.models[model].define_starting_point(self.starting_point, dataset_name)

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

        logchi2_out = (-0.5) * self.ndof * np.log(2 * np.pi)
        logchi2_out_copy = logchi2_out*1.

        for model in self.common_models.itervalues():
            logchi2_out += model.return_priors(theta)

        delayed_lnlk_computation = []

        for dataset_name, dataset in self.dataset_dict.iteritems():

            logchi2_gp_model = None

            dataset.model_reset()
            variable_values = dataset.convert(theta)
            dataset.compute(variable_values)

            if 'none' in dataset.models or 'None' in dataset.models:
                continue
            if not dataset.models:
                continue

            for model_name in dataset.models:
                if hasattr(self.models[model_name], 'common_jitter'):
                    common_ref = self.models[model_name].common_ref
                    variable_values = self.common_models[common_ref].convert(theta)
                    self.models[model_name].compute(variable_values, dataset)

            for model_name in dataset.models:

                logchi2_out += self.models[model_name].return_priors(theta, dataset_name)

                if hasattr(self.models[model_name], 'internal_likelihood'):
                    logchi2_gp_model = model_name
                    continue

                if dataset.dynamical:
                    dataset.model += dynamical_output[dataset_name]
                    continue

                if self.models[model_name].common_ref:
                    """ Taking the parameter values from the common model"""
                    common_ref = self.models[model_name].common_ref
                    variable_values = self.common_models[common_ref].convert(theta)
                else:
                    """ This model has no common model reference, i.e., it is strictly connected to the dataset"""
                    variable_values = {}

                variable_values.update(self.models[model_name].convert(theta, dataset_name))
                dataset.model += self.models[model_name].compute(variable_values, dataset)

            """ Gaussian Process check MUST be the last one or the program will fail
             that's because for the GP to work we need to know the _deterministic_ part of the model 
             (i.e. the theoretical values you get when you feed your model with the parameter values) """
            if logchi2_gp_model:
                common_ref = self.models[logchi2_gp_model].common_ref
                variable_values = self.common_models[common_ref].convert(theta)
                variable_values.update(self.models[logchi2_gp_model].convert(theta, dataset_name))

                """ GP Log-likelihood is not computed now because a single matrix must be created with 
                the joined dataset"""
                if hasattr(self.models[logchi2_gp_model], 'delayed_lnlk_computation'):
                    self.models[logchi2_gp_model].add_internal_dataset(variable_values, dataset,
                                                                   reset_status=delayed_lnlk_computation)
                    delayed_lnlk_computation.append(logchi2_gp_model)
                else:
                    logchi2_out += self.models[logchi2_gp_model].lnlk_compute(variable_values, dataset)
            else:
                logchi2_out += dataset.model_logchi2()

        """ In case there is more than one GP model"""
        for logchi2_gp_model in delayed_lnlk_computation:
            logchi2_out += self.models[logchi2_gp_model].lnlk_compute()

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

    def results_resumen(self, theta, skip_theta=False, compute_lnprob=False):
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

        if compute_lnprob:
            print
            print '===================================================================================================='
            print '===================================================================================================='
            print

            if len(np.shape(theta)) == 2:
                n_samples, n_values = np.shape(theta)
                logchi2_collection = np.zeros(n_samples)
                for i in xrange(0,n_samples):
                    logchi2_collection[i] = self(theta[i, :])
                perc0, perc1, perc2 = np.percentile(logchi2_collection, [15.865, 50, 84.135], axis=0)
                print ' LN probability: %12f   %12f %12f (15-84 p) ' % (perc1, perc0-perc1, perc2-perc1)
            else:
                print ' LN probability: %12f ' % (self(theta))

        print
        print '===================================================================================================='
        print '     ------------------------------------------------------------------------------------------     '
        print '===================================================================================================='
        print
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

    def get_model(self, theta, bjd_dict):
        model_out = {}
        model_x0 = {}

        delayed_lnlk_computation = {}

        if self.dynamical_model is not None:
            """ check if any keyword ahas get the output model from the dynamical tool
            we must do it here because all the planet are involved"""
            dynamical_output_x0 = self.dynamical_model.compute(self, theta, bjd_dict['full']['x0_plot'])
            dynamical_output = self.dynamical_model.compute(self, theta)

        for dataset_name, dataset in self.dataset_dict.iteritems():
            x0_plot = bjd_dict[dataset_name]['x0_plot']
            n_input = np.size(x0_plot)
            model_out[dataset_name] = {}
            model_x0[dataset_name] = {}
            dataset.model_reset()

            variable_values = dataset.convert(theta)
            dataset.compute(variable_values)

            for model_name in dataset.models:
                common_ref = self.models[model_name].common_ref
                variable_values = self.common_models[common_ref].convert(theta)
                if hasattr(self.models[model_name], 'common_jitter'):
                    self.models[model_name].compute(variable_values, dataset)
                if hasattr(self.models[model_name], 'common_offset'):
                    dataset.model += self.models[model_name].compute(variable_values, dataset)

            model_out[dataset_name]['systematics'] = dataset.model.copy()
            model_out[dataset_name]['jitter'] = dataset.jitter.copy()
            model_out[dataset_name]['complete'] = dataset.model.copy()

            model_x0[dataset_name]['complete'] = np.zeros(n_input, dtype=np.double)

            if 'none' in dataset.models or 'None' in dataset.models:
                continue
            if not dataset.models:
                continue

            logchi2_gp_model = None

            for model_name in dataset.models:

                if hasattr(self.models[model_name], 'internal_likelihood'):
                    logchi2_gp_model = model_name
                    continue

                if dataset.dynamical:
                    dataset.model += dynamical_output[dataset_name]
                    model_out[dataset_name][model_name] = dynamical_output[dataset_name].copy()
                    model_out[dataset_name]['complete'] += dynamical_output[dataset_name]

                    model_x0[dataset_name][model_name] = dynamical_output_x0[dataset_name].copy()
                    model_x0[dataset_name]['complete'] += dynamical_output_x0[dataset_name]
                    continue

                if self.models[model_name].common_ref:
                    """ Taking the parameter values from the common model"""
                    common_ref = self.models[model_name].common_ref
                    variable_values = self.common_models[common_ref].convert(theta)
                else:
                    """ This model has no common model reference, i.e., it is strictly connected to the dataset"""
                    variable_values = {}
                variable_values.update(self.models[model_name].convert(theta, dataset_name))

                dataset.model += self.models[model_name].compute(variable_values, dataset)

                model_out[dataset_name][model_name] = self.models[model_name].compute(variable_values, dataset)
                model_out[dataset_name]['complete'] += model_out[dataset_name][model_name]

                if hasattr(self.models[model_name], 'not_time_dependant'):
                    model_x0[dataset_name][model_name] = np.zeros(np.size(x0_plot), dtype=np.double)
                else:
                    model_x0[dataset_name][model_name] = \
                        self.models[model_name].compute(variable_values, dataset, x0_plot)
                    model_x0[dataset_name]['complete'] += model_x0[dataset_name][model_name]

            """ Gaussian Process check MUST be the last one or the program will fail
             that's because for the GP to work we need to know the _deterministic_ part of the model 
             (i.e. the theoretical values you get when you feed your model with the parameter values) """
            if logchi2_gp_model:
                common_ref = self.models[logchi2_gp_model].common_ref
                variable_values = self.common_models[common_ref].convert(theta)
                variable_values.update(self.models[logchi2_gp_model].convert(theta, dataset.name_ref))

                if hasattr(self.models[logchi2_gp_model], 'delayed_lnlk_computation'):
                    self.models[logchi2_gp_model].add_internal_dataset(variable_values, dataset,
                                                                   reset_status=delayed_lnlk_computation)
                    delayed_lnlk_computation[dataset.name_ref] = logchi2_gp_model

                else:
                    model_out[dataset_name][logchi2_gp_model] = \
                        self.models[logchi2_gp_model].sample_conditional(variable_values, dataset)
                    model_out[dataset_name]['complete'] += model_out[dataset_name][logchi2_gp_model]

                    model_x0[dataset_name][logchi2_gp_model], var  = \
                        self.models[logchi2_gp_model].sample_predict(variable_values, dataset, x0_plot)

                    model_x0[dataset_name][logchi2_gp_model + '_std'] = np.sqrt(var)
                    model_x0[dataset_name]['complete'] += model_x0[dataset_name][logchi2_gp_model]

        for dataset_name, logchi2_gp_model in delayed_lnlk_computation.iteritems():
            model_out[dataset_name][logchi2_gp_model] = \
                self.models[logchi2_gp_model].sample_conditional(dataset)
            model_out[dataset_name]['complete'] += model_out[dataset_name][logchi2_gp_model]

            model_x0[dataset_name][logchi2_gp_model], var = \
                self.models[logchi2_gp_model].sample_predict(dataset, x0_plot)

            model_x0[dataset_name][logchi2_gp_model + '_std'] = np.sqrt(var)
            model_x0[dataset_name]['complete'] += model_x0[dataset_name][logchi2_gp_model]

        # workaround to avoid memory leaks from GP module
        #gc.collect()

        return model_out, model_x0


def print_theta_bounds(i_dict, theta, bounds, skip_theta=False):
    format_string = "%10s  %4d  %12f ([%10f, %10f])"
    format_string_notheta = "%10s  %4d  ([%10f, %10f])"
    format_string_long = "%10s  %4d  %12f   %12f %12f (15-84 p) ([%9f, %9f])"

    for var, i in i_dict.iteritems():
        i_dict
        if skip_theta:
            print format_string_notheta % (var, i, bounds[i, 0], bounds[i, 1])
        elif len(np.shape(theta)) == 2:

            theta_med = compute_value_sigma(theta[:, i])

            perc0, perc1, perc2 = np.percentile(theta[:, i], [15.865, 50, 84.135], axis=0)
            print format_string_long %(var, i, theta_med[0], theta_med[2], theta_med[1], bounds[i, 0], bounds[i, 1])
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
    # Default values, taken from the PyPolyChord wrapper in PolyChord official distribution, V1.9
    polychord_parameters = {'nlive_mult': 25,
                                 'num_repeats_mult': 5,
                                 'feedback': 1,
                                 'precision_criterion': 0.001,
                                 'max_ndead': -1,
                                 'boost_posterior': 0.0,
                                 'read_resume': True,
                                 'base_dir': 'polychord/',
                                 'shutdown_jitter': False}
    polychord_dir_output = None

    def polychord_priors(self, cube):
        theta = (self.bounds[:, 1] - self.bounds[:, 0]) * cube + self.bounds[:, 0]
        return theta.tolist()

    def polychord_call(self, theta1):
        #theta = np.empty(self.ndim)
        #for i in xrange(0, self.ndim):
        #    theta[i] = theta1[i]
        theta = [theta1[i] for i in xrange(0, self.ndim)]
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


"""




class ModelContainerDnest4(ModelContainer):
    def from_prior(self):
        return self.starting_point

    def perturb(self, params):
        
        #Unlike in C++, this takes a numpy array of parameters as input,
        #and modifies it in-place. The return value is still logH.
        
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

"""