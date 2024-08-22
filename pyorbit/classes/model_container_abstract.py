from pyorbit.subroutines.common import *
from packaging.version import parse as parse_version
import pyorbit
import time

__all__ = ["ModelContainer"]

class ModelContainer(object):
    """

    """

    def __init__(self):

        self.pyorbit_version = pyorbit.__version__

        # FIXME: change the renaming and link directly the variables in constants module

        # These dictionaries are dedicated to RVs and T0s computed via dynamical integration
        # Since the individual contribution of each planet cannot be disentagled, the output of
        # the dynamical integration cannot be associated to a model.
        self.dynamical_dict = {}
        self.dynamical_t0_dict = {}
        self.dynamical_model = None

        self.dataset_dict = {}

        # Dictionaries including the common models and data-specific models
        self.common_models = {}
        self.models = {}

        self.include_priors = True

        self.ndata = None
        self.ndof = None

        self.starting_point = None
        self.starting_point_flag = False
        self.recenter_bounds_flag = True

        self.bounds = None
        self.spaces = None
        self.priors = None
        self.range = None
        self.ndim = 0
        self.pam_names = ''

        self.Tref = None

        self.ordered_planets = {}

    def model_setup(self):

        # Fixing back-compatibility issues within  v8 sub-versions
        if not hasattr(self, 'ordered_planets'):
            self.ordered_planets = {}

        # First step: setting up the correct associations between models and dataset

        for model_name, model in self.common_models.items():
            try:
                model_conf = model.model_conf
                if model_conf is None:
                    model_conf = {}
            except:
                model_conf = {}

            model.initialize_model(self, **model_conf)

        for model_name, model in self.models.items():

            try:
                model_conf = model.model_conf
                if model_conf is None:
                    model_conf = {}
            except:
                model_conf = {}

            model.initialize_model(self, **model_conf)
            model.change_parameter_status(self, **model_conf)

            for dataset_name in list(OrderedSet(model_conf) & OrderedSet(self.dataset_dict)):
                model.initialize_model_dataset(
                    self, self.dataset_dict[dataset_name], **model_conf)

        if self.dynamical_model:
            self.dynamical_model.to_be_initialized = True

            # self.dynamical_model.prepare(self)

    def boundaries_setup(self):
        # This routine setup the boundary array and at the same time
        # define a dictionary with the name of the arrays and their
        # positions in bounds/theta array so that they can be accessed
        # without using nested counters

        self.ndim = 0
        output_lists = {'bounds': [],
                        'spaces': [],
                        'priors': [],
                        }

        for model_name, model in self.models.items():
            if len(model.common_ref) > 0:
                for common_ref in model.common_ref:
                    model.default_bounds.update(
                        self.common_models[common_ref].default_bounds)
                    model.default_spaces.update(
                        self.common_models[common_ref].default_spaces)
                    model.default_priors.update(
                        self.common_models[common_ref].default_priors)
                    self.ndim, output_lists = self.common_models[common_ref].define_parameter_properties(
                        self.ndim, output_lists, model.list_pams_common)
            else:
                pass


        for dataset_name, dataset in self.dataset_dict.items():
            self.ndim, output_lists = dataset.define_parameter_properties(
                self.ndim, output_lists, dataset.list_pams)

            for model_name in dataset.models:
                self.ndim, output_lists = self.models[model_name].define_parameter_properties(
                    self.ndim, output_lists, dataset.name_ref)


        self.bounds = np.asarray(output_lists['bounds'])
        self.spaces = output_lists['spaces']
        self.priors = output_lists['priors']
        self.range = self.bounds[:, 1] - self.bounds[:, 0]


    def initialize_logchi2(self):

        # Second step: define the number of variables and setting up the boundaries
        # To be done only the first time ever
        self.ndata = 0
        for dataset_name, dataset in self.dataset_dict.items():
            if not dataset.models:
                continue
            self.ndata += dataset.n
        self.ndof = self.ndata - self.ndim

    def starting_points_setup(self):

        self.starting_point = np.average(self.bounds, axis=1)
        self.starting_point[:] = None

        for model_name, model in self.common_models.items():
            model.define_starting_point(self.starting_point)

        for dataset_name, dataset in self.dataset_dict.items():
            dataset.define_starting_point(self.starting_point)

            for model in dataset.models:
                self.models[model].define_starting_point(
                    self.starting_point, dataset_name)

    def check_bounds(self, theta):

        for ii in range(0, self.ndim):
            if not (self.bounds[ii, 0] < theta[ii] < self.bounds[ii, 1]):
                return False

        period_storage = []
        period_storage_ordered = np.zeros(len(self.ordered_planets))

        for model_name, model in self.common_models.items():

            if model.model_class == 'planet':

                """ Step 1: retrieve the planet period"""
                period = model.transformation['P'](theta, model.fixed,
                                                   model.parameter_index['P'])

                """ Step 2: save the all planet periods into a list"""
                period_storage.extend([period])

                """ Step 3: save the period of the planet in the ordered list"""
                if model_name in self.ordered_planets:
                    period_storage_ordered[self.ordered_planets[model_name]] = period
                    #print('    ', model_name, self.ordered_planets[model_name], period_storage_ordered)

                """ Step 4: check if the eccentricity is within the given range"""
                if 'e' in model.parameter_index:
                    e = model.transformation['e'](theta,
                                                  model.fixed,
                                                  model.parameter_index['e'])

                    if not model.bounds['e'][0] <= e < model.bounds['e'][1]:
                        # print('eccentricity>1')
                        # print()
                        return False

                """ Step 5: check if the impact parameter is below 1 + Rp/Rs """
                if 'b' in model.parameter_index and 'R_Rs' in model.parameter_index:
                    b = model.transformation['b'](theta,
                                                  model.fixed,
                                                  model.parameter_index['b'])
                    R = model.transformation['R_Rs'](theta,
                                                  model.fixed,
                                                  model.parameter_index['R_Rs'])
                    if not b <= 1 + R:
                        return False

                """ Step 6 eclipse depth must be greater than the amplitude of
                phase variation (during the eclipse we only have the stellar line) """
                if 'phase_amp' in model.parameter_index and 'delta_occ' in model.parameter_index:
                    phase_amp = model.transformation['delta_occ'](theta,
                                                  model.fixed,
                                                  model.parameter_index['phase_amp'])
                    delta_occ = model.transformation['delta_occ'](theta,
                                                  model.fixed,
                                                  model.parameter_index['delta_occ'])
                    if phase_amp > delta_occ: return False


        #""" Step 6 eclipse depth must be greater than the amplitude of
        #phase variation (during the eclipse we only have the stellar line) """
        #for dataset_name, dataset in self.dataset_dict.items():
        #    if 'phase_amp' in dataset.parameter_index and 'delta_occ' in dataset.parameter_index:
        #        phase_amp = dataset.transformation['delta_occ'](theta,
        #                                        dataset.fixed,
        #                                        dataset.parameter_index['phase_amp'])
        #        delta_occ = dataset.transformation['delta_occ'](theta,
        #                                        dataset.fixed,
        #                                        dataset.parameter_index['delta_occ'])
        #        if phase_amp > delta_occ: return False

        """ Step 7 check for overlapping periods (within 2.5% arbitrarily chosen)"""
        for i_n, i_v in enumerate(period_storage):
            if i_n == len(period_storage) - 1:
                break
            if np.amin(np.abs(period_storage[i_n + 1:] - i_v)) / i_v < 0.025:
                # print('overlapping periods  detected')
                # print()
                return False

        """ Step 8 check if the planet are ordered"""
        for i_n, i_v in enumerate(period_storage_ordered):

            if i_n == len(period_storage_ordered) - 1:
                break
            if np.amin(period_storage_ordered[i_n + 1:] - i_v) < 0.0:
                # print('inverted order detected')
                # print()
                return False

        return True

    def __call__(self, theta, include_priors=True):

        #start = time.time()

        log_priors, log_likelihood = self.log_priors_likelihood(theta)
        #end = time.time()
        #print("Computation took {0:.16f} seconds".format(end - start))

        if self.include_priors and include_priors:
            return log_priors + log_likelihood
        else:
            return log_likelihood

    def log_priors_likelihood(self, theta, return_priors=True):

        log_priors = 0.00
        log_likelihood = 0.00

        """
        Constant term added either by dataset.model_logchi2() or gp.log_likelihood()
        """
        if not self.check_bounds(theta):
            if return_priors is False:
                return -np.inf
            else:
                return -np.inf, -np.inf

        if self.dynamical_model is not None:
            """ check if any keyword ahas get the output model from the dynamical tool
            we must do it here because all the planet are involved"""
            dynamical_output = self.dynamical_model.compute(self, theta)

        for model_name, model in self.common_models.items():
            log_priors += model.return_priors(theta)

        delayed_lnlk_computation = []
        residuals_analysis = {}

        for dataset_name, dataset in self.dataset_dict.items():

            logchi2_gp_model = None
            compute_gp_residuals = False

            dataset.model_reset()
            parameter_values = dataset.convert(theta)
            dataset.compute(parameter_values)

            log_priors += dataset.return_priors(theta)

            if 'none' in dataset.models or 'None' in dataset.models:
                continue
            if not dataset.models:
                continue

            skip_loglikelihood = False

            for model_name in dataset.models:

                log_priors += self.models[model_name].return_priors(
                    theta, dataset_name)

                parameter_values = {}
                for common_ref in self.models[model_name].common_ref:
                    parameter_values.update(
                        self.common_models[common_ref].convert(theta))

                #TODO: remove try-except starting from version 11 !!
                try:
                    for planet_name in self.models[model_name].multiple_planets:
                        parameter_values.update(
                            self.common_models[planet_name].convert_with_name(theta, planet_name))
                except:
                    pass

                parameter_values.update(
                    self.models[model_name].convert(theta, dataset_name))

                """ residuals will be computed following the definition in Dataset class
                    This section has never been tested
                """
                if getattr(self.models[model_name], 'residuals_analysis', False):

                    skip_loglikelihood = True
                    if self.models[model_name].gp_before_correlation:
                        compute_gp_residuals = True


                    try:
                        residuals_analysis[model_name]['parameter_values'] = parameter_values
                    except:
                        residuals_analysis[model_name] = {'parameter_values': parameter_values}

                    if dataset_name == self.models[model_name].x_dataset:
                        residuals_dataset_label = 'x'
                        residuals_analysis[model_name]['x_gp_model'] = None
                        residuals_analysis[model_name]['x_gp_parameters'] = None
                    else:
                        # Subtract the priors previously added, to avoid including the priors twice
                        log_priors -= self.models[model_name].return_priors(theta, dataset_name)

                        residuals_dataset_label = 'y'
                        residuals_analysis[model_name]['y_gp_model'] = None
                        residuals_analysis[model_name]['y_gp_parameters'] = None
                    continue

                if getattr(self.models[model_name], 'internal_likelihood', False):
                    logchi2_gp_model = model_name
                    continue

                # if getattr(self.models[model_name], 'model_class', None) is 'common_jitter':
                if getattr(self.models[model_name], 'jitter_model', False):
                    dataset.jitter += self.models[model_name].compute(
                        parameter_values, dataset)
                    continue

                if getattr(dataset, 'dynamical', False):
                    dataset.external_model = dynamical_output[dataset_name]

                if dataset.normalization_model is None and (self.models[model_name].unitary_model or self.models[model_name].normalization_model):
                    dataset.normalization_model = np.ones(dataset.n, dtype=np.double)

                if self.models[model_name].unitary_model:
                    dataset.unitary_model += self.models[model_name].compute(
                    parameter_values, dataset)
                elif self.models[model_name].normalization_model:
                    dataset.normalization_model *= self.models[model_name].compute(
                        parameter_values, dataset)
                else:
                    dataset.additive_model += self.models[model_name].compute(
                        parameter_values, dataset)

            dataset.compute_model()
            dataset.compute_residuals()

            """ Gaussian Process check MUST be the last one or the program will fail
            that's because for the GP to work we need to know the _deterministic_ part of the model
            (i.e. the theoretical values you get when you feed your model with the parameter values) """

            if logchi2_gp_model:

                parameter_values = {}
                for common_ref in self.models[logchi2_gp_model].common_ref:
                    parameter_values.update(
                        self.common_models[common_ref].convert(theta))

                parameter_values.update(
                    self.models[logchi2_gp_model].convert(theta, dataset_name))

                """ GP Log-likelihood is not computed now because a single matrix must be
                    computed with the joint dataset"""
                if hasattr(self.models[logchi2_gp_model], 'delayed_lnlk_computation'):

                    self.models[logchi2_gp_model].add_internal_dataset(parameter_values, dataset)
                    if logchi2_gp_model not in delayed_lnlk_computation:
                        delayed_lnlk_computation.append(logchi2_gp_model)
                elif skip_loglikelihood:
                    residuals_analysis[model_name][residuals_dataset_label+'_gp_model'] = logchi2_gp_model
                    residuals_analysis[model_name][residuals_dataset_label+'_gp_parameters'] = parameter_values
                else:
                    log_likelihood += self.models[logchi2_gp_model].lnlk_compute(
                        parameter_values, dataset)

                if compute_gp_residuals:
                    dataset.residuals_for_regression -= self.models[logchi2_gp_model].sample_predict(parameter_values, dataset)
                    residuals_analysis[model_name][residuals_dataset_label+'_gp_model'] = None

            elif not skip_loglikelihood:
                log_likelihood += dataset.model_logchi2()

        """ Correlation between residuals of two datasets through orthogonal distance regression
            it must be coomputed after any other model has been removed from the
            independent variable, if not provided as ancillary dataset
        """
        for model_name in residuals_analysis:
            parameter_values =  residuals_analysis[model_name]['parameter_values']
            x_dataset = self.dataset_dict[self.models[model_name].x_dataset]
            y_dataset = self.dataset_dict[self.models[model_name].y_dataset]
            modelout_xx, modelout_yy = self.models[model_name].compute(parameter_values, x_dataset, y_dataset)

            x_dataset.residuals -= modelout_xx
            if residuals_analysis[model_name]['x_gp_model']:
                logchi2_gp_model = residuals_analysis[model_name]['x_gp_model']
                parameter_values = residuals_analysis[model_name]['x_gp_parameters']
                log_likelihood += self.models[logchi2_gp_model].lnlk_compute(parameter_values, x_dataset)
            else:
                log_likelihood += x_dataset.model_logchi2()

            y_dataset.residuals -= modelout_yy
            if residuals_analysis[model_name]['y_gp_model']:
                logchi2_gp_model = residuals_analysis[model_name]['y_gp_model']
                parameter_values = residuals_analysis[model_name]['y_gp_parameters']
                log_likelihood += self.models[logchi2_gp_model].lnlk_compute(parameter_values, y_dataset)
            else:
                log_likelihood += y_dataset.model_logchi2()


        """ In case there is more than one GP model """
        for logchi2_gp_model in delayed_lnlk_computation:
            log_likelihood += self.models[logchi2_gp_model].lnlk_compute()

        """ check for finite log_priors and log_likelihood"""
        if np.isnan(log_priors) or np.isnan(log_likelihood):
            log_likelihood = -np.inf
            log_priors = -np.inf

        if return_priors is False:
            return log_likelihood
        else:
            return log_priors, log_likelihood

    def recenter_bounds(self, pop_mean, recenter=True):
        # This function recenters the bounds limits for circular parameters
        # Also, it extends the range of a parameter if the output of PyDE is a fixed number

        ind_list = []

        for model_name, model in self.common_models.items():
            ind_list.extend(model.index_recenter_bounds())

        for dataset_name, dataset in self.dataset_dict.items():
            for model in dataset.models:
                ind_list.extend(
                    self.models[model].index_recenter_bounds(dataset.name_ref))

        if not recenter:
            return ind_list

        if ind_list:
            tmp_range = (self.bounds[:, 1] - self.bounds[:, 0]) / 2
            replace_bounds = np.zeros([self.ndim, 2])
            replace_bounds[:, 0] = pop_mean - tmp_range
            replace_bounds[:, 1] = pop_mean + tmp_range
            self.bounds[ind_list, :] = replace_bounds[ind_list, :]

    def fix_population(self, pop_mean, population):

        ind_list = self.recenter_bounds(pop_mean, recenter=False)
        n_pop = np.size(population, axis=0)

        if ind_list:
            for par_ind in ind_list:
                fix_sel = (population[:, par_ind] <= self.bounds[par_ind, 0]) | (
                    population[:, par_ind] >= self.bounds[par_ind, 1])
                population[fix_sel, par_ind] = pop_mean[par_ind]

        for ii in range(0, self.ndim):
            if np.amax(population[:, ii]) - np.amin(population[:, ii]) < 10e-14:
                average_pops = np.average(population[:, ii])
                population[:, ii] = np.random.normal(
                    average_pops, 10e-12, n_pop)

        return population


