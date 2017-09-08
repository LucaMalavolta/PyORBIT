from common import *
from compute_RV import ComputeKeplerian, ComputeDynamical


class ModelContainer:
    def __init__(self):
        self.dataset_dict = {}
        self.dataset_index = {}
        self.n_datasets = 0

        self.t0_list = {}

        self.keplerian_model = ComputeKeplerian()
        self.dynamical_model = ComputeDynamical()

        # pyde/emcee variables
        self.emcee_parameters = {'nsave': 0, 'npop_mult': 2, 'thin': 1,
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

        self.models = {}

        self.model_list = []
        self.bounds_list = []

        self.starting_point = None
        self.starting_point_flag = False
        self.recenter_bounds_flag = True

        self.planet_name = ''

        self.variable_list = {'common': {}}
        self.bound_list = []
        self.bounds = None
        self.range = None
        self.ndim = 0
        self.pam_names = ''
        self.star_mass = [1.0000, 0.1000]
        self.star_radius = [1.0000, 0.1000]

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

    def create_model_list(self):

        for dataset in self.dataset_dict.itervalues():
            for data_model in dataset.models:
                if not (data_model in self.model_list):
                    self.model_list.append(data_model)


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
        self.create_bounds()

        self.ndata = 0
        for dataset in self.dataset_dict.itervalues():
            if 'none' in dataset.models or 'None' in dataset.models:
                continue
            self.ndata += dataset.n
        self.ndof = self.ndata - self.ndim

        self.initialize()

    def create_bounds(self):
        # This routine creates the boundary array and at the same time
        # creates a dictionary with the name of the arrays and their
        # positions in bounds/theta array so that they can be accessed
        # without using nested counters

        self.ndim = 0

        for dataset in self.dataset_dict.itervalues():
            dataset.define_bounds(self)

        for model in self.models.itervalues():
            model.define_bounds(self)

        self.bounds = np.asarray(self.bounds_list)
        self.range = self.bounds[:, 1] - self.bounds[:, 0]

    def initialize(self):
        # Function with two goals:
        # * Unfold and print out the output from theta
        # * give back a parameter name associated to each value in the result array

        self.pam_names = self.ndim * ['']

        for dataset in self.dataset_dict.itervalues():
            dataset.initialize(self)

        for model in self.models.itervalues():
            model.initialize(self)

    def create_starting_point(self):

        self.starting_point = np.average(self.bounds, axis=1)

        for dataset in self.dataset_dict.itervalues():
            dataset.starting_point(self)

        for model in self.models.itervalues():
            model.starting_point(self)

    def check_bounds(self, theta):
        for ii in xrange(0, self.ndim):
            if not (self.bounds[ii, 0] < theta[ii] < self.bounds[ii, 1]):
                return False

        period_storage = []
        for pl_name in self.models['planets'].planet_name:
            """ Step 1: save the all planet periods into a list"""
            period_storage.extend(
                [self.models['planets'].variables[pl_name]['P'](theta, self.models['planets'].fixed,
                                                               self.models['planets'].var_list[pl_name]['P'])])

            """ Step 2: check if the eccentricity is within the given range"""
            e = self.models['planets'].variables[pl_name]['e'](theta,
                                                              self.models['planets'].fixed,
                                                              self.models['planets'].var_list[pl_name]['e'])
            if not self.models['planets'].bounds[pl_name]['e'][0] <= e < self.models['planets'].bounds[pl_name]['e'][1]:
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
        logchi2_out = 2. * self.ndof * np.log(2 * np.pi)

        for model in self.models.itervalues():
            logchi2_out += model.return_priors(theta)

            "dynamical computations are performed and stored inside the planets classes"
            if model.model_class == 'planets':
                model.store_dynamical_model(self, theta)

        #if 'sinusoid' in self.model_list:
        #    logchi2_out += giveback_priors(
        #        self.scv.prior_kind['Prot'], self.scv.prior_pams['Prot'], theta[self.variable_list['Common']['Prot']])

        for dataset_name, dataset in self.dataset_dict.items():
            dataset.model_reset()
            dataset.model_offset(theta[self.variable_list[dataset_name]['offset']])
            dataset.model_jitter(theta[self.variable_list[dataset_name]['jitter']])
            dataset.model_linear(theta[self.variable_list[dataset_name]['linear']])

            if 'none' in dataset.models or 'None' in dataset.models:
                continue

            logchi2_gp_model = None
            for model in dataset.models:

                if self.models[model].model_class == 'gaussian_process':
                    logchi2_gp_model = model
                    continue
                else:
                    dataset.model += self.models[model].compute(theta, dataset)

                #if 'sinusoids' in dataset.models:
                #    dataset.model += self.scv.compute(self, theta, dataset)

            if 'Tcent' in dataset.kind:
                if dataset.planet_name in self.models['planets'].dynamical:
                    """ we have dynamical computations, so we include them in the model"""
                    dataset.model += self.models['planets'].dynamical_output[dataset_name]
                else:
                    dataset.model += dataset.compute(self, theta)

            # Gaussian Process check MUST be the last one or the program will fail
            # if 'gaussian_process' in dataset.models:
            if logchi2_gp_model is not None:
                logchi2_out += self.models[logchi2_gp_model].return_priors(theta, dataset_name)
                logchi2_out += self.models[logchi2_gp_model].lnlk_compute(theta, dataset)
            else:
                logchi2_out += dataset.model_logchi2()

             # workaround to avoid memory leaks from GP module
        gc.collect()

        return logchi2_out

    def results_resumen(self, theta, verbose=True):
        # Function with two goals:
        # * Unfold and print out the output from theta
        # * give back a parameter name associated to each value in the result array

        for dataset in self.dataset_dict.itervalues():
            dataset.print_vars(self, theta)

        for model in self.models.itervalues():
            model.print_vars(self, theta)

    def rv_make_model(self, theta, x_range, x_phase):
        # it return the RV model for a single planet, after removing the activity from the RV curve and removing
        # the offsets between the datasets

        model_actv = {}
        model_plan = {}
        model_dsys = {}
        model_curv = {}

        ''' NEW!
        create fake dataset with dataset.x0 and all the required parameters
        '''
        model_orbs = {'BJD': x_range * 0.0, 'pha': x_phase * 0.0}
        model_curv['BJD'] = x_range * 0.0
        model_curv['pha'] = x_phase * 0.0
        model_plan['BJD'] = {}
        model_plan['pha'] = {}

        if bool(self.models['planets'].dynamical):
            """Check if the dynamical option has been activated: full RV curve will be computed using
            the dynamical integrator"""
            dyn_output = self.dynamical_model.compute(self, self.models['planets'], theta)
            dyn_output_fullorbit = self.dynamical_model.compute(self, self.models['planets'],
                                                                theta,
                                                                full_orbit=(x_range - self.Tref))
            model_orbs['BJD'] += dyn_output_fullorbit['full_orbit']
            print 'WARNING: phase plot generated using non-interacting model!!!'
            print model_orbs['BJD'][:10]
        # computing the orbit for the full dataset

        for pl_name in self.models['planets'].planet_name:

            dyn_flag = (pl_name in self.models['planets'].dynamical)
            if dyn_flag:
                dict_pams = self.models['planets'].kepler_from_dynamical(self, theta, pl_name)
            else:
                dict_pams = self.models['planets'].convert(pl_name, theta)

            model_plan['BJD'][pl_name] = self.keplerian_model.model(dict_pams, x_range - self.Tref)
            model_plan['pha'][pl_name] = self.keplerian_model.model(dict_pams, x_phase * dict_pams['P'])

            model_orbs['pha'] += model_plan['pha'][pl_name]
            if not dyn_flag:
                model_orbs['BJD'] += model_plan['BJD'][pl_name]

        for model in self.models.itervalues():
            if model.model_class == 'curvature':
                curv_pams = model.convert(theta)
                model_curv['BJD'] += model.model_curvature(curv_pams, x_range - self.Tref)

        for dataset_name, dataset in self.dataset_dict.items():

            model_actv[dataset_name] = np.zeros(dataset.n)
            model_orbs[dataset_name] = np.zeros(dataset.n)
            model_curv[dataset_name] = np.zeros(dataset.n)
            model_plan[dataset_name] = {}

            dataset.model_reset()
            dataset.model_jitter(theta[self.variable_list[dataset_name]['jitter']])
            dataset.model_offset(theta[self.variable_list[dataset_name]['offset']])
            dataset.model_linear(theta[self.variable_list[dataset_name]['linear']])

            model_dsys[dataset_name] = dataset.model.copy()

            if 'none' in dataset.models or 'None' in dataset.models:
                continue

            '''
            for model in dataset.models:

                    dataset_out[model] = self.models[model].compute(theta, dataset)
            '''

            for model in dataset.models:
                if model == 'planets':
                    if bool(self.models['planets'].dynamical):
                        """Dynamical models (with all the interacting planets included in the resulting RVs)"""
                        model_orbs[dataset_name] += dyn_output[dataset_name]
                        dataset.model += model_orbs[dataset_name]

                    for pl_name in self.models['planets'].planet_name:
                        dyn_flag = (pl_name in self.models['planets'].dynamical)
                        if dyn_flag:
                            dict_pams = self.models['planets'].kepler_from_dynamical(self, theta, pl_name)
                            model_plan[dataset_name][pl_name] = self.keplerian_model.model(dict_pams, dataset.x0)
                        else:
                            model_plan[dataset_name][pl_name] = \
                                self.keplerian_model.compute(self.models['planets'], theta, dataset, pl_name)
                            model_orbs[dataset_name] += model_plan[dataset_name][pl_name]
                            dataset.model += model_plan[dataset_name][pl_name]

                elif self.models[model].model_class == 'curvature':
                    model_curv[dataset_name] += self.models[model].model_curvature(curv_pams, dataset.x0)
                    """ Saving to dataset model in case it's required by GP computation"""
                    dataset.model += model_curv[dataset_name]
                elif self.models[model].model_class == 'correlation':
                    model_actv[dataset_name] += self.models[model].compute(theta, dataset)
                    dataset.model += model_actv[dataset_name]
                elif self.models[model].model_class == 'gaussian_process':
                    model_actv[dataset_name] += self.models[model].sample_compute(theta, dataset)

            #if 'sinusoids' in dataset.models:
            #    model_actv[dataset_name] += self.scv.compute(self, theta, dataset)
            #    dataset.model += model_actv[dataset_name]

        return model_dsys, model_plan, model_orbs, model_actv, model_curv

    def recenter_bounds(self, pop_mean, population):
        # This function recenters the bounds limits for circular variables
        # Also, it extends the range of a variable if the output of PyDE is a fixed number
        ind_list = []
        n_pop = np.size(population, axis=0)
        if 'planets' in self.model_list:
            for pl_name in self.models['planets'].planet_name:

                if 'esino' in self.variable_list[pl_name]:
                    esino_list = self.variable_list[pl_name]['esino']
                    ecoso_list = self.variable_list[pl_name]['ecoso']
                    e_pops = population[:, esino_list] ** 2 + population[:, ecoso_list] ** 2
                    o_pops = np.arctan2(population[:, esino_list], population[:, ecoso_list], dtype=np.double)
                    # e_mean = (self.models['planets'].bounds[pl_name]['e'][0] +
                    # self.models['planets'].bounds[pl_name]['e'][1]) / 2.
                    for ii in xrange(0, n_pop):
                        if not self.models['planets'].bounds[pl_name]['e'][0] + 0.02 <= e_pops[ii] < \
                                        self.models['planets'].bounds[pl_name]['e'][1] - 0.02:
                            e_random = np.random.uniform(self.models['planets'].bounds[pl_name]['e'][0],
                                                         self.models['planets'].bounds[pl_name]['e'][1])
                            population[ii, esino_list] = np.sqrt(e_random) * np.sin(o_pops[ii])
                            population[ii, ecoso_list] = np.sqrt(e_random) * np.cos(o_pops[ii])

                if 'f' in self.variable_list[pl_name]:
                    ind_list.append(self.variable_list[pl_name]['f'])

                if 'o' in self.variable_list[pl_name]:
                    ind_list.append(self.variable_list[pl_name]['o'])

                if 'lN' in self.variable_list[pl_name]:
                    ind_list.append(self.variable_list[pl_name]['lN'])

        if 'sinusoids' in self.model_list:
            for jj in range(0, self.scv.n_seasons):
                ind_list.extend(self.variable_list[self.scv.season_name[jj] + '_pha'])
            for dataset in self.dataset_dict.itervalues():
                for jj in range(0, self.scv.n_seasons):
                    if dataset.season_flag[jj]:
                        # ind_list.extend(self.variable_list[dataset.planet_name][self.scv.season_name[jj] + '_amp'])
                        if self.scv.use_offset[dataset.kind]:
                            ind_list.append(self.variable_list[dataset.kind][self.scv.season_name[jj] + '_off'])

        if np.size(ind_list) > 0:
            tmp_range = (self.bounds[:, 1] - self.bounds[:, 0]) / 2
            for var_ind in ind_list:
                self.bounds[var_ind, :] = pop_mean[var_ind] + [-tmp_range[var_ind], tmp_range[var_ind]]
                fix_sel = (population[:, var_ind] <= self.bounds[var_ind, 0]) | (
                    population[:, var_ind] >= self.bounds[var_ind, 1])
                population[fix_sel, var_ind] = pop_mean[var_ind]

        for ii in xrange(0, self.ndim):
            if np.amax(population[:, ii]) - np.amin(population[:, ii]) < 10e-7:
                range_restricted = (self.bounds[ii, 1] - self.bounds[ii, 0]) / 100.
                min_bound = np.maximum((pop_mean[ii] - range_restricted / 2.0), self.bounds[ii, 0])
                max_bound = np.minimum((pop_mean[ii] + range_restricted / 2.0), self.bounds[ii, 1])
                population[:, ii] = np.random.uniform(min_bound, max_bound, n_pop)


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
