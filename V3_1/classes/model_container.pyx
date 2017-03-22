from common import *
from planets import PlanetsCommonVariables
from gaussian import GaussianProcessCommonVariables
from curvature import CurvatureCommonVariables
from sinusoids import SinusoidsCommonVariables
from dataset import Dataset, TransitCentralTimes
sys.path.append('/home/malavolta/CODE/trades/pytrades')
#sys.path.append('/Users/malavolta/Astro/CODE/trades/pytrades')
from pytrades_lib import pytrades
import constants

class ModelContainer:
    def __init__(self):
        self.dataset_list = []
        self.n_datasets = 0

        self.t0_list = {}

        self.scv = SinusoidsCommonVariables()
        self.pcv = PlanetsCommonVariables()
        self.gcv = GaussianProcessCommonVariables()
        self.ccv = CurvatureCommonVariables()

        # pyde/emcee variables
        self.ngen = 0
        self.nsave = 0
        self.nsteps = 0
        self.nburn = 0
        self.npop_multi = 0
        self.nwalkers = 0
        self.thin = 1

        self.model_list = []
        self.bounds_list = []

        self.starting_point = None
        self.starting_point_flag = False
        self.recenter_bounds_flag = True

        self.planet_name = ''

        self.variable_list = {'Common': {}}
        self.bound_list = []
        self.bounds = 0
        self.ndim = 0
        self.pam_names = ''
        self.star_mass = [1.0000,  0.1000]
        self.star_radius = [1.0000,  0.1000]

        """
        Values have been taken from TRADES
        These variables will be renamed in the next release, right now I'm keeping the original names
        to avoid breaking the code
        """
        self.G_grav = constants.Gsi # Gravitational Constants in SI system [m^3/kg/s^2]
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
        self.n_datasets = np.size(self.dataset_list)
        for dataset in self.dataset_list:

            if 'sinusoids' in dataset.models:
                self.scv.setup_model_sinusoids(dataset)

            dataset.model_reset()
            for data_model in dataset.models:
                if not (data_model in self.model_list):
                    self.model_list.append(data_model)

    def create_bounds(self):
        # This routine creates the boundary array and at the same time
        # creates a dictionary with the name of the arrays and their
        # positions in bounds/theta array so that they can be accessed
        # without using nested counters

        self.ndim = 0

        for dataset in self.dataset_list:
            dataset.define_bounds(self)

        if 'kepler' in self.model_list:
            self.pcv.define_bounds(self)

        if 'curvature' in self.model_list:
            self.ccv.define_bounds(self)

        if 'sinusoids' in self.model_list:
            self.scv.define_bounds(self)

        if 'gaussian' in self.model_list:
            self.gcv.define_bounds(self)

        self.bounds = np.asarray(self.bounds_list)

    def create_starting_point(self):

        self.starting_point = np.average(self.bounds, axis=1)

        for dataset in self.dataset_list:
            dataset.starting_point(self)

        if 'kepler' in self.model_list:
            self.pcv.starting_point(self)

        if 'curvature' in self.model_list:
            self.ccv.starting_point(self)

        if 'sinusoids' in self.model_list:
            self.scv.starting_point(self)

        if 'gaussian' in self.model_list:
            self.gcv.starting_point(self)

    def check_bounds(self, theta):
        for ii in xrange(0, self.ndim):
            if not (self.bounds[ii, 0] < theta[ii] < self.bounds[ii, 1]):
                return False
        for pl_name in self.pcv.planet_name:
            # if not bool(self.pcv.circular[pl_name]):
            e = self.pcv.variables[pl_name]['e'](theta, self.pcv.fixed, self.pcv.var_list[pl_name]['e'])
            if not self.pcv.bounds[pl_name]['e'][0] <= e < self.pcv.bounds[pl_name]['e'][1]:
                return False

        return True

    def __call__(self, theta):
        if not self.check_bounds(theta):
            return -np.inf
        chi2_out = 0.0

        if 'kepler' in self.model_list:
            for pl_name in self.pcv.planet_name:
                chi2_out += self.pcv.return_priors(pl_name, theta)

        if bool(self.pcv.dynamical):
            """ check if any keyword ahas get the output model from the dynamical tool
            we must do it here because all the planet are involved"""
            dyn_output = self.pcv.compute_dynamical(self, theta)

        if 'curvature' in self.model_list:
            chi2_out += self.ccv.return_priors(theta)

        if 'sinusoid' in self.model_list:
            chi2_out += giveback_priors(
                self.scv.prior_kind['Prot'], self.scv.prior_pams['Prot'], theta[self.variable_list['Common']['Prot']])

        if 'gaussian' in self.model_list:
            chi2_out += self.gcv.return_priors(theta)

        for dataset in self.dataset_list:
            dataset.model_reset()
            dataset.model_offset(theta[self.variable_list[dataset.name_ref]['offset']])
            dataset.model_jitter(theta[self.variable_list[dataset.name_ref]['jitter']])
            dataset.model_linear(theta[self.variable_list[dataset.name_ref]['linear']])

            if 'kepler' in dataset.models:
                if bool(self.pcv.dynamical):
                    """ we have dynamical computations, so we include them in the model"""
                    dataset.model += dyn_output[dataset.name_ref]
                for pl_name in self.pcv.planet_name:
                    """ we check if there is any planet which model has been obtained by assuming non-intercating
                    keplerians, and then we compute the expected RVs"""
                    pass
                    if pl_name not in self.pcv.dynamical:
                        dataset.model += self.pcv.compute(theta, dataset, pl_name)

            if 'sinusoids' in dataset.models:
                dataset.model += self.scv.compute(self, theta, dataset)

            if 'curvature' in dataset.models:
                dataset.model += self.ccv.compute(theta, dataset)

            if 'Tcent' in dataset.models:
                if dataset.planet_name in self.pcv.dynamical:
                    """ we have dynamical computations, so we include them in the model"""
                    dataset.model += dyn_output[dataset.name_ref]
                else:
                    dataset.model += dataset.compute(self, theta)

            # Gaussian Process check MUST be the last one or the program will fail
            if 'gaussian' in dataset.models:
                chi2_out += self.gcv.return_priors(theta, dataset.name_ref)
                chi2_out += self.gcv.lnlk_compute(theta, dataset)
            else:
                chi2_out += dataset.model_logchi2()

        return chi2_out

    def results_resumen(self, theta):
        # Function with two goals:
        # * Unfold and print out the output from theta
        # * give back a parameter name associated to each value in the result array

        self.pam_names = self.ndim * ['']

        for dataset in self.dataset_list:
            dataset.print_vars(self, theta)
            print

        if 'kepler' in self.model_list:
            self.pcv.print_vars(self, theta)
            print

        if 'curvature' in self.model_list:
            self.ccv.print_vars(self, theta)
            print

        if 'sinusoids' in self.model_list:
            self.scv.print_vars(self, theta)
            print

        if 'gaussian' in self.model_list:
            self.gcv.print_vars(self, theta)
            print

    def rv_make_model(self, theta, x_range, x_phase):
        # it return the RV model for a single planet, after removing the activity from the RV curve and removing
        # the offsets between the datasets

        model_actv = {}
        model_plan = {}
        model_dsys = {}
        model_curv = {}

        model_orbs = {'BJD': x_range * 0.0, 'pha':x_phase * 0.0}
        model_curv['BJD'] = x_range * 0.0
        model_curv['pha'] = x_phase * 0.0
        model_plan['BJD'] = {}
        model_plan['pha'] = {}

        if bool(self.pcv.dynamical):
            """Check if the dynamical option has been activated: full RV curve will be computed using
            the dynamical integrator"""
            dyn_output = self.pcv.compute_dynamical(self, theta)
            dyn_output_fullorbit = self.pcv.compute_dynamical(self, theta, full_orbit=(x_range-self.Tref))
            model_orbs['BJD'] += dyn_output_fullorbit['full_orbit']
            print 'WARNING: phase plot generated using non-interacting model!!!'
            print model_orbs['BJD'][:10]
        # computing the orbit for the full dataset
        for pl_name in self.pcv.planet_name:

            dyn_flag = (pl_name in self.pcv.dynamical)
            if dyn_flag:
                dict_pams = self.pcv.kepler_from_dynamical(self, theta, pl_name)
            else:
                dict_pams = self.pcv.convert(pl_name, theta)

            model_plan['BJD'][pl_name] = self.pcv.model_kepler(dict_pams, x_range - self.Tref)
            model_plan['pha'][pl_name] = self.pcv.model_kepler(dict_pams, x_phase * dict_pams['P'])

            model_orbs['pha'] += model_plan['pha'][pl_name]
            if not dyn_flag:
                model_orbs['BJD'] += model_plan['BJD'][pl_name]

        if 'curvature' in self.model_list:
            curv_pams = self.ccv.convert(theta)
            model_curv['BJD'] = self.ccv.model_curvature(curv_pams, x_range - self.Tref)

        for dataset in self.dataset_list:

            model_actv[dataset.name_ref] = np.zeros(dataset.n)
            model_orbs[dataset.name_ref] = np.zeros(dataset.n)
            model_curv[dataset.name_ref] = np.zeros(dataset.n)
            model_plan[dataset.name_ref] = {}

            dataset.model_reset()
            dataset.model_offset(theta[self.variable_list[dataset.name_ref]['offset']])
            dataset.model_linear(theta[self.variable_list[dataset.name_ref]['linear']])

            model_dsys[dataset.name_ref] = dataset.model

            if 'curvature' in dataset.models:
                model_curv[dataset.name_ref] = self.ccv.model_curvature(curv_pams, dataset.x0)

            if 'kepler' in dataset.models:
                if bool(self.pcv.dynamical):
                    """Dynamical models (with all the interacting planets included in the resulting RVs)"""
                    model_orbs[dataset.name_ref] += dyn_output[dataset.name_ref]

                for pl_name in self.pcv.planet_name:
                    dyn_flag = (pl_name in self.pcv.dynamical)
                    if dyn_flag:
                        dict_pams = self.pcv.kepler_from_dynamical(self, theta, pl_name)
                        model_plan[dataset.name_ref][pl_name] = self.pcv.model_kepler(dict_pams, dataset.x0)
                    else:
                        model_plan[dataset.name_ref][pl_name] = self.pcv.compute(theta, dataset, pl_name)
                        model_orbs[dataset.name_ref] += model_plan[dataset.name_ref][pl_name]

            if 'sinusoids' in dataset.models:
                model_actv[dataset.name_ref] += self.scv.compute(self, theta, dataset)

        return model_dsys, model_plan, model_orbs, model_actv, model_curv

    # This function recenters the bounds limits for circular variables
    # Also, it extends the range of a variable if the output of PyDE is a fixed number
    def recenter_bounds(self, pop_mean, population):
        ind_list = []

        if 'kepler' in self.model_list:
            for pl_name in self.pcv.planet_name:

                if 'esino' in self.variable_list[pl_name]:
                    esino_list = self.variable_list[pl_name]['esino']
                    ecoso_list = self.variable_list[pl_name]['ecoso']
                    e_pops = population[:, esino_list] ** 2 + population[:, ecoso_list] ** 2
                    o_pops = np.arctan2(population[:, esino_list], population[:, ecoso_list], dtype=np.double)
                    # e_mean = (self.pcv.bounds[pl_name]['e'][0] + self.pcv.bounds[pl_name]['e'][1]) / 2.
                    for ii in xrange(0, self.nwalkers):
                        if not self.pcv.bounds[pl_name]['e'][0] + 0.02 <= e_pops[ii] < \
                                        self.pcv.bounds[pl_name]['e'][1] - 0.02:
                            e_random = np.random.uniform(self.pcv.bounds[pl_name]['e'][0],
                                                         self.pcv.bounds[pl_name]['e'][1])
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
            for dataset in self.dataset_list:
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
                population[:, ii] = np.random.uniform(min_bound, max_bound, self.nwalkers)


class ModelContainerPolyChord(ModelContainer):
    def polychord_priors(self, cube):
        theta = (self.bounds[:, 1] - self.bounds[:, 0]) * cube + self.bounds[:, 0]
        #print theta.tolist()
        return theta.tolist()

    def polychord_call(self, theta1):
        theta = np.empty(self.ndim)
        for i in xrange(0, self.ndim):
            theta[i] = theta1[i]
        #print theta
        phi = [0.0] * self.nDerived
        chi_out = self(theta)
        if chi_out < -0.5e10:
            return -0.5e10, phi
        return chi_out, phi


class ModelContainerMultiNest(ModelContainer):
    def multinest_priors(self, cube, ndim, nparams):
        #cube[:] = (self.bounds[:, 1] - self.bounds[:, 0]) * cube[:] + self.bounds[:, 0]
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


class ModelContainerSPOTPY(ModelContainer):
    def spotpy_priors(self, theta_med):
        pams_priors = []
        for ii in xrange(0, self.ndim):
            pams_priors.append(
                spotpy.parameter.Uniform(
                    'theta' + repr(ii), self.bounds[ii, 0], self.bounds[ii, 1], optguess=theta_med[ii]))

        if 'kepler' in self.model_list:
            for pl_name in self.pcv.planet_name:
                for key in self.pcv.prior_pams[pl_name]:
                    ind = self.variable_list[pl_name][key]
                    pams_priors[ind] = self.assign_priors(self.pcv.prior_kind[pl_name][key], 'theta' + repr(ind),
                                                          self.pcv.prior_pams[pl_name][key], theta_med[ind])

        if 'sinusoids' in self.model_list:
            if 'Prot' in self.scv.prior_pams:
                ind = self.variable_list['Common']['Prot']
                pams_priors[ind] = self.assign_priors(self.scv.prior_kind[key], 'theta' + repr(ind),
                                                      self.scv.prior_pams[key], theta_med[ind])

        if 'gaussian' in self.model_list:
            for key in self.gcv.prior_pams['Common']:
                ind = self.variable_list[key]
                pams_priors[ind] = self.assign_priors(self.gcv.prior_kind[key], 'theta' + repr(ind),
                                                      self.gcv.prior_pams['Common'][key], theta_med[ind])
            for dataset in self.dataset_list:
                for key in self.gcv.prior_pams[dataset.name_ref]:
                    ind = self.variable_list[dataset.name_ref][key]
                    pams_priors[ind] = self.assign_priors(self.pcv.prior_kind[dataset.name_ref][key],
                                                          'theta' + repr(ind),
                                                          self.pcv.prior_pams[dataset.name_ref][key], theta_med[ind])
        return pams_priors

    def assign_priors(self, key_name, prior_name, prior_val, prior_med):
        if key_name == 'Normal' or key_name == 'Gaussian':
            return spotpy.parameter.Normal(prior_name, prior_val[0], prior_val[1], optguess=prior_med)
        if key_name == 'logNormal':
            return spotpy.parameter.logNormal(prior_name, prior_val[0], prior_val[1], optguess=prior_med)
        if key_name == 'Exponential':
            return spotpy.parameter.Exponential(prior_name, prior_val[0], optguess=prior_med)
        return spotpy.parameter.Uniform(prior_name, prior_val[0], prior_val[1], optguess=prior_med)

    def output_concatenated(self):
        output_array = []
        for dataset in self.dataset_list:
            output_array = np.concatenate((output_array, dataset.y), axis=0)
        return output_array

    def __call__(self, theta):
        output_array = []
        for dataset in self.dataset_list:
            dataset.model_reset()
            dataset.model_offset(theta[self.variable_list[dataset.name_ref]['offset']])
            dataset.model_jitter(theta[self.variable_list[dataset.name_ref]['jitter']])
            dataset.model_linear(theta[self.variable_list[dataset.name_ref]['linear']])

            if 'kepler' in dataset.models:
                for pl_name in self.pcv.planet_name:
                    dataset.model += self.pcv.compute(theta, dataset, pl_name)

            if 'sinusoids' in dataset.models:
                dataset.model += self.scv.compute(self, theta, dataset)

            if 'Tcent' in dataset.models:
                dataset.model += dataset.compute(self, theta)

            # Gaussian Process check MUST be the last one or the program will fail
            if 'gaussian' in dataset.models:
                dataset.model += self.gcv.sample_compute(theta, dataset)
            output_array = np.concatenate((output_array, dataset.model), axis=0)

        return output_array
