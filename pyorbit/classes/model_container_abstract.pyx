from common import *

__all__ = ["ModelContainer"]


class ModelContainer(object):
    """

    """

    def __init__(self):

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
        self.use_threading_pool = True

        self.bounds = None
        self.spaces = None
        self.priors = None
        self.range = None
        self.ndim = 0
        self.pam_names = ''
        self.star_mass = [1.0000, 0.1000]
        self.star_radius = [1.0000, 0.1000]

        self.Tref = None

    def model_setup(self):
        # First step: setting up the correct associations between models and dataset

        for model_name, model in self.models.iteritems():
            try:
                model.initialize_model(self, **model.model_conf)

                for dataset_name in list(set(model.model_conf) & set(self.dataset_dict)):
                    model.setup_dataset(self.dataset_dict[dataset_name], **model.model_conf)
            except:
                model.initialize_model(self, **{})

        for dataset_name, dataset in self.dataset_dict.iteritems():
            for model_name in dataset.models:
                self.models[model_name].common_initialization_with_dataset(dataset)
                try:
                    for common_model in self.models[model_name].common_ref:
                        self.common_models[common_model].common_initialization_with_dataset(dataset)
                except:
                    pass

        if self.dynamical_model:
            self.dynamical_model.prepare(self)

    def create_variables_bounds(self):
        # This routine creates the boundary array and at the same time
        # creates a dictionary with the name of the arrays and their
        # positions in bounds/theta array so that they can be accessed
        # without using nested counters

        self.ndim = 0
        output_lists = {'bounds': [],
                        'spaces': [],
                        'priors': [],
                        }

        for model in self.models.itervalues():
            if len(model.common_ref) > 0:

                for common_ref in model.common_ref:
                    model.default_bounds.update(self.common_models[common_ref].default_bounds)
                    model.default_spaces.update(self.common_models[common_ref].default_spaces)
                    model.default_priors.update(self.common_models[common_ref].default_priors)
                    self.ndim, output_lists = self.common_models[common_ref].define_variable_properties(
                        self.ndim, output_lists, model.list_pams_common)

            else:
                pass

        for dataset in self.dataset_dict.itervalues():
            self.ndim, output_lists = dataset.define_variable_properties(self.ndim, output_lists, dataset.list_pams)

            for model_name in dataset.models:
                self.ndim, output_lists = self.models[model_name].define_variable_properties(self.ndim, output_lists, dataset.name_ref)

        self.bounds = np.asarray(output_lists['bounds'])
        self.spaces = output_lists['spaces']
        self.priors = output_lists['priors']
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
        for model_name, model in self.common_models.items():

            if model.model_class == 'planet':

                """ Step 1: save the all planet periods into a list"""
                period_storage.extend(
                    [model.transformation['P'](theta, model.fixed,
                                               model.variable_index['P'])])

                """ Step 2: check if the eccentricity is within the given range"""
                e = model.transformation['e'](theta,
                                              model.fixed,
                                              model.variable_index['e'])
                if not model.bounds['e'][0] <= e < model.bounds['e'][1]:
                    return False

        """ Step 4 check ofr overlapping periods (within 2.5% arbitrarily chosen)"""
        for i_n, i_v in enumerate(period_storage):
            if i_n == len(period_storage) - 1: break
            if np.amin(np.abs(period_storage[i_n + 1:] - i_v)) / i_v < 0.025:
                return False

        return True

    def __call__(self, theta, include_priors=True):
        log_priors, log_likelihood = self.log_priors_likelihood(theta)

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

        for model in self.common_models.itervalues():
            log_priors += model.return_priors(theta)

        delayed_lnlk_computation = []

        for dataset_name, dataset in self.dataset_dict.iteritems():

            logchi2_gp_model = None

            dataset.model_reset()
            variable_values = dataset.convert(theta)
            dataset.compute(variable_values)

            log_priors += dataset.return_priors(theta)

            if 'none' in dataset.models or 'None' in dataset.models:
                continue
            if not dataset.models:
                continue

            for model_name in dataset.models:

                log_priors += self.models[model_name].return_priors(theta, dataset_name)

                variable_values = {}
                for common_ref in self.models[model_name].common_ref:
                    variable_values.update(self.common_models[common_ref].convert(theta))

                #try:
                #    """ Taking the parameter values from the common models"""
                #    for common_ref in self.models[model_name].common_ref:
                #        variable_values.update(self.common_models[common_ref].convert(theta))
                #except:
                #    """ This model has no common model reference, i.e., it is strictly connected to the dataset"""
                #    pass

                variable_values.update(self.models[model_name].convert(theta, dataset_name))

                """ residuals will be computed following the definition in Dataset class
                """

                if getattr(self.models[model_name], 'internal_likelihood', False):
                    logchi2_gp_model = model_name
                    continue

                if getattr(self.models[model_name], 'model_class', None) is 'common_jitter':
                    dataset.jitter = self.models[model_name].compute(variable_values, dataset)
                    continue

                if getattr(dataset, 'dynamical', False):
                    dataset.external_model = dynamical_output[dataset_name]

                if getattr(self.models[model_name], 'unitary_model', False):
                    dataset.unitary_model += self.models[model_name].compute(variable_values, dataset)
                    if dataset.normalization_model is None:
                        dataset.normalization_model = np.ones(dataset.n, dtype=np.double)

                elif getattr(self.models[model_name], 'normalization_model', False):
                    if dataset.normalization_model is None:
                        dataset.normalization_model = np.ones(dataset.n, dtype=np.double)
                    dataset.normalization_model *= self.models[model_name].compute(variable_values, dataset)
                else:
                    dataset.additive_model += self.models[model_name].compute(variable_values, dataset)

            dataset.compute_model()
            dataset.compute_residuals()

            """ Gaussian Process check MUST be the last one or the program will fail
             that's because for the GP to work we need to know the _deterministic_ part of the model 
             (i.e. the theoretical values you get when you feed your model with the parameter values) """
            if logchi2_gp_model:

                variable_values = {}
                for common_ref in self.models[logchi2_gp_model].common_ref:
                        variable_values.update(self.common_models[common_ref].convert(theta))

                variable_values.update(self.models[logchi2_gp_model].convert(theta, dataset_name))

                """ GP Log-likelihood is not computed now because a single matrix must be created with 
                the joined dataset"""
                if hasattr(self.models[logchi2_gp_model], 'delayed_lnlk_computation'):

                    self.models[logchi2_gp_model].add_internal_dataset(variable_values, dataset,
                                                                   reset_status=delayed_lnlk_computation)
                    delayed_lnlk_computation.append(logchi2_gp_model)
                else:
                    log_likelihood += self.models[logchi2_gp_model].lnlk_compute(variable_values, dataset)
            else:
                log_likelihood += dataset.model_logchi2()

        """ In case there is more than one GP model"""
        for logchi2_gp_model in delayed_lnlk_computation:
            log_likelihood += self.models[logchi2_gp_model].lnlk_compute()

        if return_priors is False:
            return log_likelihood
        else:
            return log_priors, log_likelihood

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

