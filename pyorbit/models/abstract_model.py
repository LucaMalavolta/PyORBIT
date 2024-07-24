from pyorbit.subroutines.common import \
    np,\
    get_var_exp,\
    get_var_val,\
    get_fix_val,\
    get_2darray_from_val,\
    giveback_priors,\
    nested_sampling_prior_prepare, OrderedSet

from pyorbit.subroutines.common import get_var_exp, get_var_log
from pyorbit.subroutines.common import get_var_exp_base2, get_var_log_base2
from pyorbit.subroutines.common import get_var_exp_base10, get_var_log_base10
from pyorbit.subroutines.common import get_var_exp_natural, get_var_log_natural
from pyorbit.subroutines.common import get_var_arccosine, get_var_cosine, get_var_arcsine, get_var_sine
from pyorbit.subroutines import constants

class AbstractModel(object):

    def __init__(self, model_name, common_ref):
        self.model_name = model_name

        try:
            if len(common_ref) > 0:
                self.common_ref = np.atleast_1d(common_ref).tolist()
            else:
                self.common_ref = []
        except (NameError, TypeError):
            self.common_ref = []

        """ Set of parameters shared among all the datasets relying on the same
        physical (i.e., common) objects. Parameters belonging to different common
        objects can be listed together.
        When making new objects, always pay attention to avoid duplicate names
        for parameters
        """
        self.list_pams_common = OrderedSet()

        """ Dataset-specific parameters must be listed here. A given model will
        be re-computed for each dataset using the corresponding values of the
        parameters listed here, while all the other parameters will be taken by the
        common objects of reference. As an example, all Gaussian Process object
        will have the same period (listed in list_pams_common) but each dataset
        will be characterized by its own covariance amplitude (listed in
        list_pams_dataset).
        """
        self.list_pams_dataset = OrderedSet()

        """ For some circular parameters it is convenient to recenter the
        boundaries so that tha mode is near the center of the interval defined
        by the boundaries, in such a way the computation of the median and the
        confidence interval are less prone to errors
        """
        self.recenter_pams_dataset = OrderedSet()

        self.unitary_model = False
        self.normalization_model = False

        self.planet_ref = common_ref
        self.stellar_ref = 'star_parameters'

        self.planet_list = [] # used only by very complex models

        self.sampler_parameters = {}

        self.transformation = {}
        self.parameter_index = {}
        self.bounds = {}
        self.parameters = {}

        self.starts = {}

        self.fix_list = {}
        self.fixed = []
        self.nfix = 0

        self.prior_kind = {}
        self.prior_pams = {}
        self.spaces = {}

        self.default_bounds = {}
        self.default_spaces = {}
        self.default_priors = {}

        self.model_conf = None

        self.multivariate_priors = {}
        self.multivariate_pams = {}
        self.multivariate_func = {}
        self.multivariate_med = {}
        self.multivariate_cov = {}

        self.include_zero_point = False
        self.exclude_zero_point = False

        ## New attribute
        self.multiple_planets = False

    def initialize_model(self, mc, **kwargs):
        pass

    def change_parameter_status(self, mc, **kwargs):

        dataset_pams = kwargs.get('dataset_parameters', [])
        for par in dataset_pams:
            self.list_pams_common.discard(par)
            self.list_pams_dataset.update([par])
        common_pams = kwargs.get('common_parameters', [])
        for par in common_pams:
            self.list_pams_dataset.discard(par)
            self.list_pams_common.update([par])

    def initialize_model_dataset(self, mc, dataset, **kwargs):
        pass

    def define_parameter_properties(self, ndim, output_lists, dataset_name):
        """[summary]
            Boundaries are defined in this class, where all the dataset-related
            parameters are stored. Bounds and parameter index CANNOT be defined
            in the Common class: we don't know a priori which parameters
            will be actually used in the complete model.

        Args:
            ndim ([type]): [description]
            output_lists ([type]): [description]
            dataset_name ([type]): [description]

        Returns:
            [type]: [description]
        """

        self.transformation[dataset_name] = {}
        self.parameter_index[dataset_name] = {}
        self.sampler_parameters[dataset_name] = {}

        if dataset_name not in self.bounds.keys():
            self.bounds[dataset_name] = {}

        if dataset_name not in self.spaces.keys():
            self.spaces[dataset_name] = {}

        for par in self.list_pams_dataset:
            if par not in self.bounds[dataset_name]:
                self.bounds[dataset_name][par] = self.default_bounds[par]

            if par not in self.spaces[dataset_name]:
                self.spaces[dataset_name][par] = self.default_spaces[par]

            if par in self.fix_list[dataset_name]:

                self.transformation[dataset_name][par] = get_fix_val
                self.parameter_index[dataset_name][par] = self.nfix
                self.prior_kind[dataset_name][par] = 'None'
                self.prior_pams[dataset_name][par] = []

                # Workaround to preserve compatibility with Python 2.x
                if isinstance(self.fix_list[dataset_name][par],
                            type('string')) \
                        and par in self.default_fixed:
                    self.fixed.append(
                        get_2darray_from_val(self.default_fixed[par])[0])
                else:
                    self.fixed.append(self.fix_list[dataset_name][par][0])

                self.nfix += 1
            else:
                if self.spaces[dataset_name][par] == 'Linear':
                    self.transformation[dataset_name][par] = get_var_val
                    output_lists['bounds'].append(
                        self.bounds[dataset_name][par])
                elif self.spaces[dataset_name][par] == 'Log_Natural':
                    self.transformation[dataset_name][par] = get_var_exp_natural
                    output_lists['bounds'].append(
                        np.log(self.bounds[dataset_name][par]))
                elif self.spaces[dataset_name][par] == 'Log_Base2':
                    self.transformation[dataset_name][par] = get_var_exp_base2
                    output_lists['bounds'].append(
                        np.log2(self.bounds[dataset_name][par]))
                elif self.spaces[dataset_name][par] == 'Log_Base10':
                    self.transformation[dataset_name][par] = get_var_exp_base10
                    output_lists['bounds'].append(
                        np.log10(self.bounds[dataset_name][par]))
                elif self.spaces[dataset_name][par] == 'Logarithmic':
                    self.transformation[dataset_name][par] = get_var_exp
                    output_lists['bounds'].append(
                        np.log2(self.bounds[dataset_name][par]))
                elif self.spaces[dataset_name][par] == 'Sine_Angle':
                    self.transformation[dataset_name][par] = get_var_arcsine
                    output_lists['bounds'].append(
                        np.sin(self.bounds[dataset_name][par]*constants.deg2rad))
                elif self.spaces[dataset_name][par] == 'Cosine_Angle':
                    self.transformation[dataset_name][par] = get_var_arccosine
                    output_lists['bounds'].append(
                        np.cos(self.bounds[dataset_name][par]*constants.deg2rad))

                if par not in self.prior_pams[dataset_name]:
                    self.prior_kind[dataset_name][par] = \
                        self.default_priors[par][0]
                    self.prior_pams[dataset_name][par] = \
                        self.default_priors[par][1]

                nested_coeff = \
                    nested_sampling_prior_prepare(
                        self.prior_kind[dataset_name][par],
                        output_lists['bounds'][-1],
                        self.prior_pams[dataset_name][par],
                        self.spaces[dataset_name][par])

                output_lists['spaces'].append(self.spaces[dataset_name][par])
                output_lists['priors'].append(
                    [self.prior_kind[dataset_name][par],
                        self.prior_pams[dataset_name][par],
                        nested_coeff])

                self.parameter_index[dataset_name][par] = ndim
                self.sampler_parameters[dataset_name][par] = ndim
                ndim += 1

        return ndim, output_lists

    def define_starting_point_from_derived(self, starting_point, dataset_name, par):
        return False

    def define_starting_point(self, starting_point, dataset_name):

        if not bool(self.starts):
            return

        for par in self.starts[dataset_name]:
            if self.define_starting_point_from_derived(
                    starting_point, dataset_name, par):
                continue

            if self.spaces[dataset_name][par] == 'Linear':
                start_converted = self.starts[dataset_name][par]
            elif self.spaces[dataset_name][par] == 'Log_Natural':
                start_converted = np.log(self.starts[dataset_name][par])
            elif self.spaces[dataset_name][par] == 'Log_Base2':
                start_converted = np.log2(self.starts[dataset_name][par])
            elif self.spaces[dataset_name][par] == 'Log_Base10':
                start_converted = np.log10(self.starts[dataset_name][par])
            elif self.spaces[dataset_name][par] == 'Logarithmic':
                start_converted = np.log2(self.starts[dataset_name][par])
            elif self.spaces[dataset_name][par] == 'Sine_Angle':
                start_converted = np.arcsin(self.starts[dataset_name][par])*constants.rad2deg
            elif self.spaces[dataset_name][par] == 'Cosine_Angle':
                start_converted = np.arccos(self.starts[dataset_name][par])*constants.rad2deg

            starting_point[self.sampler_parameters[dataset_name]
                            [par]] = start_converted

    def convert(self, theta, dataset_name):
        parameter_values = {}
        # If we need the parameters for the prior, we are not providing any
        # name for the dataset

        for par in self.list_pams_dataset:
            parameter_values[par] = self.transformation[dataset_name][par](
                theta, self.fixed, self.parameter_index[dataset_name][par])
        return parameter_values

    def return_priors(self, theta, dataset_name):
        prior_out = 0.00
        parameter_values = self.convert(theta, dataset_name)

        """ Preserving back-compatibility with version 8
        #TODO: to be simplified in the next version
        """

        if getattr(self, 'multivariate_priors', False):
            if len(OrderedSet(self.multivariate_pams[dataset_name]) & OrderedSet(self.list_pams_dataset[dataset_name])) > 0:
                multi_par = [parameter_values[ii]
                             for ii in self.multivariate_pams[dataset_name]]
                pdf = self.multivariate_func[dataset_name].pdf(multi_par)
                if pdf > 0:
                    prior_out += np.log(
                        self.multivariate_func[dataset_name].pdf(multi_par))
                else:
                    return -np.inf
            else:
                self.multivariate_pams[dataset_name] = []
        else:
            self.multivariate_pams = {dataset_name: []}

        for par in self.list_pams_dataset:

            if par in self.multivariate_pams[dataset_name]:
                continue

            prior_out += giveback_priors(
                self.prior_kind[dataset_name][par],
                self.bounds[dataset_name][par],
                self.prior_pams[dataset_name][par],
                parameter_values[par])

        return prior_out


    def index_recenter_bounds(self, dataset_name):
        ind_list = []
        for par in list(OrderedSet(self.recenter_pams_dataset)
                        & OrderedSet(self.sampler_parameters[dataset_name])):
            ind_list.append(self.sampler_parameters[dataset_name][par])

        return ind_list


    def compute(self, parameter_values, dataset, x0_input=None):
        return 0.00


    def transfer_parameter_properties(self, mc, dataset, par_original, par_subset, keywords={},
    common_pam=False, dataset_pam=False):

        if common_pam:
            self.list_pams_common.update([par_subset])
        elif dataset_pam:
            self.list_pams_dataset.update([par_subset])
        else:
            print('ERROR in transfer_prior subroutine')
            print('The new variable must be either common- or dataset-defined')
            quit()

        for common_model in self.common_ref:
            if par_original not in mc.common_models[common_model].list_pams:
                continue

            if common_pam:
                mc.common_models[common_model].list_pams.update([par_subset])

            if par_original in keywords.get('boundaries', {}):
                par_update = keywords['boundaries'][par_original]
            elif par_subset in keywords.get('boundaries', {}):
                par_update = keywords['boundaries'][par_subset]
            elif par_original in mc.common_models[common_model].bounds:
                par_update = mc.common_models[common_model].bounds[par_original]
            elif par_original in self.bounds[dataset.name_ref]:
                par_update = self.bounds[dataset.name_ref][par_original]
            else:
                par_update = mc.common_models[common_model].default_bounds[par_original]

            if common_pam:
                mc.common_models[common_model].bounds.update({par_subset: par_update})
            else:
                self.bounds[dataset.name_ref].update({par_subset: par_update})

            if par_original in keywords.get('spaces', {}):
                par_update = keywords['spaces'][par_original]
            elif par_subset in keywords.get('spaces', {}):
                par_update = keywords['spaces'][par_subset]
            elif par_original in mc.common_models[common_model].spaces:
                par_update = mc.common_models[common_model].spaces[par_original]
            elif par_original in self.spaces[dataset.name_ref]:
                par_update = self.spaces[dataset.name_ref][par_original]
            else:
                par_update = mc.common_models[common_model].default_spaces[par_original]

            if common_pam:
                mc.common_models[common_model].spaces.update({par_subset: par_update})
            else:
                self.spaces[dataset.name_ref].update({par_subset: par_update})

            if par_original in keywords.get('priors', {}):
                prior_pams = np.atleast_1d(keywords['priors'][par_original])
                par_update0 = prior_pams[0]
                try:
                    par_update1 = np.asarray(prior_pams[1:], dtype=np.double)
                except:
                    par_update1 = np.asarray([0.00], dtype=np.double)
            elif par_subset in keywords.get('priors', {}):
                prior_pams = np.atleast_1d(keywords['priors'][par_subset])
                par_update0 = prior_pams[0]
                try:
                    par_update1 = np.asarray(prior_pams[1:], dtype=np.double)
                except:
                    par_update1 = np.asarray([0.00], dtype=np.double)
            elif par_original in mc.common_models[common_model].prior_pams:
                par_update0 = mc.common_models[common_model].prior_kind[par_original]
                par_update1 = mc.common_models[common_model].prior_pams[par_original]
            elif par_original in self.prior_pams[dataset.name_ref]:
                par_update0 = self.prior_kind[dataset.name_ref][par_original]
                par_update1 = self.prior_pams[dataset.name_ref][par_original]
            else:
                par_update0 = mc.common_models[common_model].default_priors[par_original][0]
                par_update1 = mc.common_models[common_model].default_priors[par_original][1]

            if common_pam:
                mc.common_models[common_model].prior_kind.update({par_subset: par_update0})
                mc.common_models[common_model].prior_pams.update({par_subset: par_update1})
            else:
                self.prior_kind[dataset.name_ref].update({par_subset: par_update0})
                self.prior_pams[dataset.name_ref].update({par_subset: par_update1})
