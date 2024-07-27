from pyorbit.subroutines.common import np, OrderedSet
from pyorbit.common.abstract_common import AbstractCommon
from pyorbit.model_definitions import datatype_definition

from numpy.lib.recfunctions import append_fields, drop_fields

class Dataset(AbstractCommon):

    def __init__(self, model_name, kind, models):
        super(self.__class__, self).__init__(None)

        for kind_name, kind_list in datatype_definition.items():
            if kind in kind_list:
                self.kind = kind_name
                break

        # model kind:  'RV', 'PHOT', 'ACT'...
        self.models = models
        self.name_ref = model_name

        self.dynamical = False
        self.planet_name = None

        self.generic_list_pams = OrderedSet(['jitter', 'offset', 'subset'])

        self.generic_default_priors = {
            'jitter': ['Uniform'  , []],
            'offset': ['Uniform', []]}

        self.generic_default_spaces = {
            'jitter': 'Linear',
            'offset': 'Linear'}

        if self.kind == 'Tcent':
            self.generic_default_spaces['jitter'] = 'Logarithmic'
            self.generic_default_priors['jitter'] = ['Uniform', []]

        self.variable_compressed = {}
        self.variable_expanded = {}

        self.model_class = 'dataset'
        self.list_pams = OrderedSet()
        self.default_bounds = {}
        self.default_spaces = {}
        self.default_priors = {}
        self.default_fixed = {}
        self.recenter_pams = {}

        self.Tref = None
        self.residuals = None
        self.residuals_for_regression = None
        self.model = None
        # this model is compute externally and passed to the compute subroutine of the model
        self.external_model = None
        self.additive_model = None
        self.unitary_model = None
        self.normalization_model = None
        self.jitter = None
        self.mask = {}

        self.ancillary = None

        self.compute_plot = True


    def append_ancillary(self, input_file, input_array=False, input_array_str=False):
        """ Function to either read an ancillary file or pass the value of another
        dataset as ancillary

        we define as ancillary every dataset required to
        model some effect (usually, instrumental) but not object of modelling
        itself (i.e., the dataset does not enter in the log-likelihhod)
        Added in PyORBIT 9.2

        Args:
            input_file (string): name of the file, absolute or relative path
            data_input: True if the input_file is actually a ndarray
        """

        if not input_array:
            # read ancillary data from txt file when it is not passed as a ndarray
            input_array = np.genfromtxt(input_file, names=True)
            #! reading as string if needed - the selection of columns to be
            #! considered as string must happen inside the model
            input_array_str =np.genfromtxt(input_file, dtype=str)

        """ NOTE: string are only supported in ancillary data files
        """
        self.ancillary_str_index = {}
        self.ancillary_str = input_array.copy()

        if self.ancillary is not None:
            # Data ancillary has been already defined when reading the main files
            # we take the keywords from the file and add them to the existing
            for iname, name in enumerate(input_array.dtype.names):
                try:
                    self.ancillary = drop_fields(self.ancillary, name)
                    self.ancillary = append_fields(self.ancillary, name, input_array[name])

                except ValueError:
                    print('The ancillary input array is not a structured array')
                    print('https://numpy.org/doc/stable/user/basics.rec.html')
                    quit()

                try:
                    self.ancillary_str = drop_fields(self.ancillary_str, name)
                    self.ancillary_str = append_fields(self.ancillary_str, name, input_array_str[:,iname])
                except:
                    print('Ancillary file containing string must be provided as argument of append_Ancillary ')

        else:
            self.ancillary = input_array.copy()
            self.ancillary_str = input_array.copy()

            for iname, name in enumerate(input_array.dtype.names):
                self.ancillary_str = drop_fields(self.ancillary_str, name)
                self.ancillary_str = append_fields(self.ancillary_str, name, input_array_str[:,iname])


    def convert_dataset_from_file(self, input_file):

        """ Robust data reading, now encompassing the case of ancillary data
        columns embedded in the main file"""

        data0 = np.atleast_2d(np.loadtxt(input_file))
        data1 = np.genfromtxt(input_file, names=True)
        data_input = np.zeros([np.size(data0, axis=0), 6], dtype=np.double) - 1.

        if np.shape(data0)[0] == np.shape(data1)[0]:
            for i_name, v_name in enumerate(data1.dtype.names):
                if i_name<3:
                    data_input[:, i_name] = data1[v_name]
                elif v_name == 'jit' or v_name=='jitter' or v_name=='jitter_flag':
                    data_input[:, 3] = data1[v_name]
                elif v_name == 'off' or v_name=='offset' or v_name=='offset_flag':
                    data_input[:, 4] = data1[v_name]
                elif v_name == 'sub' or v_name=='subset' or v_name=='subset_flag':
                    data_input[:, 5] = data1[v_name]
                else:
                    self.ancillary = data1.copy()
                    self.ancillary_str = data1.copy()
        else:
            """ Fall back to previous behaviour of the code"""
            data_input[:, :np.size(data0, axis=1)] = data0[:, :]

        return data_input

    def define_dataset_base(self, data_input, update=False,
                            flag_shutdown_jitter=False):
        # Add a flag to save internally the input dataset

        if flag_shutdown_jitter:
            data_input[:, 4] = -1

        if not self.models:
            data_input[:, 3:] = -1

        if self.kind == 'transit_time':
            """ Special input reading from T0 files """
            self.n_transit = np.asarray(data_input[:, 0], dtype=np.int16)
            self.x = np.asarray(data_input[:, 1], dtype=np.double)
            self.e = np.asarray(data_input[:, 2], dtype=np.double)
            """ copy of self.y added for consistency with the rest of the code
            """
            self.y = self.x
        else:
            self.x = np.asarray(data_input[:, 0], dtype=np.double)
            self.y = np.asarray(data_input[:, 1], dtype=np.double)
            self.e = np.asarray(data_input[:, 2], dtype=np.double)

        self.n = np.size(self.x)
        self.n_shape = np.shape(self.y)

        if not update:
            if self.Tref is None:
                self.Tref = np.mean(self.x, dtype=np.double)

            """ Default boundaries are defined according to the characteristic
                of the dataset. They must be large enough to allow formost of
                the anomalous situations
            """
            x_range = np.max(self.x) - np.min(self.x)
            y_range = np.max(self.y) - np.min(self.y)
            y_trend = np.abs(y_range/x_range)
            y_diff = np.abs(np.mean(self.x)-self.Tref) * y_trend + 1000.
            self.generic_default_bounds = {'offset': [np.min(self.y) - 10.*y_diff, np.max(self.y) + 10.*y_diff],
                                           'jitter': [np.min(self.e)/100., 100 * np.max(self.e)]}

            if self.kind == 'Phot':
                self.generic_default_bounds['offset'][0] = - \
                    1000.0 * np.max(self.e)

            self._setup_systematic_dictionaries('jitter', data_input[:, 3])
            self._setup_systematic_dictionaries('offset', data_input[:, 4])

        self.x0 = self.x - self.Tref
        self._setup_systematic_mask('jitter', data_input[:, 3])
        self._setup_systematic_mask('offset', data_input[:, 4])

        if np.amax(data_input[:, 5]) > 0:
            sel = data_input[:, 5] >= -0.5
            self.submodel_minflag = np.int64(np.amin(data_input[sel, 5]))
            self.submodel_maxflag = np.int64(np.amax(data_input[:, 5])) + 1
            self.submodel_flag = np.int64(np.amax(data_input[:, 5])) + 1
            self.submodel_id = data_input[:, 5]
        else:
            self.submodel_minflag = None
            self.submodel_maxflag = None
            self.submodel_flag = None

        if np.amax(data_input[:, 3]) > 0:
            self.input_jitter = np.asarray(data_input[:, 3], dtype=np.double)
        else:
            self.input_jitter = None

        if np.amax(data_input[:, 4]) > 0:
            self.input_offset = np.asarray(data_input[:, 4], dtype=np.double)
        else:
            self.input_offset = None

        if np.amax(data_input[:, 5]) > 0:
            self.input_subset = np.asarray(data_input[:, 5], dtype=np.double)
        else:
            self.input_subset = None

        self.model_reset()

    def _setup_systematic_dictionaries(self, var_generic, dataset_vals):
        n_sys = np.max(dataset_vals.astype(np.int64)) + 1
        self.variable_compressed[var_generic] = {}
        for ii in range(0, n_sys):
            var = var_generic + '_' + repr(ii)
            self.list_pams.update([var])
            self.default_bounds[var] = self.generic_default_bounds[var_generic]
            self.default_spaces[var] = self.generic_default_spaces[var_generic]
            self.default_priors[var] = self.generic_default_priors[var_generic]
            self.variable_compressed[var_generic][var] = None
            self.variable_expanded[var] = var_generic

    def _setup_systematic_mask(self, var_generic, dataset_vals):
        n_sys = np.max(dataset_vals.astype(np.int64)) + 1
        if np.size(dataset_vals) == 1:
            dataset_vals = np.zeros(self.n_shape, dtype=np.int64)

        for ii in range(0, n_sys):
            var = var_generic + '_' + repr(ii)
            self.mask[var] = np.zeros(self.n_shape, dtype=bool)
            self.mask[var][(abs(dataset_vals - ii) < 0.1)] = True

    def _delete_systematic_dictionaries_mask(self, var_generic):
        if var_generic not in self.variable_compressed:
            return
        for var in self.variable_compressed[var_generic]:
            self.list_pams.discard(var)
            self.default_bounds.pop(var, None)
            self.default_spaces.pop(var, None)
            self.default_priors.pop(var, None)
            self.variable_expanded.pop(var, None)
            self.mask.pop(var)
        self.variable_compressed[var_generic] = {}

    def shutdown_jitter(self):
        self._delete_systematic_dictionaries_mask('jitter')
        self.input_jitter = None

    def shutdown_offset(self):
        self._delete_systematic_dictionaries_mask('offset')
        self.input_offset = None

    def common_Tref(self, Tref_in):
        self.Tref = Tref_in
        self.x0 = self.x - self.Tref
        return

    def model_reset(self):
        self.residuals = None
        self.model = None
        self.additive_model = np.zeros(self.n_shape, dtype=np.double)
        self.unitary_model = np.zeros(self.n_shape, dtype=np.double)
        self.normalization_model = None
        self.external_model = np.zeros(self.n_shape, dtype=np.double)
        self.jitter = np.zeros(self.n_shape, dtype=np.double)
        return

    def compute(self, variable_value):
        for var in self.list_pams:
            if self.variable_expanded[var] == 'jitter':
                self.jitter[self.mask[var]] += variable_value[var]
            elif self.variable_expanded[var] == 'offset':
                self.additive_model[self.mask[var]] += variable_value[var]

    def compute_model(self):
        if self.normalization_model is None:
            self.model = self.additive_model + self.external_model
        else:
            self.model = self.additive_model + \
                self.external_model + \
                (1. + self.unitary_model)*self.normalization_model

    def compute_model_from_arbitrary_datasets(self, additive_model, unitary_model, normalization_model, external_model):
        if normalization_model is None or unitary_model is None:
            return additive_model + external_model
        else:
            return additive_model + external_model + (1. + unitary_model)*normalization_model

    def compute_residuals(self):
        self.residuals = self.y - self.model

    def model_logchi2(self):
        env = 1.0 / (self.e ** 2.0 + self.jitter ** 2.0)

        #chi2 = -0.5 * (self.n * np.log(2 * np.pi) +
        #               np.sum(self.residuals ** 2 * env - np.log(env)))
        #print('{0:25s} {1:12f} {2:12f} \n'.format(self.name_ref, chi2, np.std(self.residuals)))

        return -0.5 * (self.n * np.log(2 * np.pi) +
                       np.sum(self.residuals ** 2 * env - np.log(env)))

    def update_bounds_spaces_priors_starts(self):

        for var_generic in list(set(self.bounds) & set(self.variable_compressed)):
            for var in self.variable_compressed[var_generic]:
                self.bounds[var] = self.bounds[var_generic]

        for var_generic in list(set(self.spaces) & set(self.variable_compressed)):
            for var in self.variable_compressed[var_generic]:
                self.spaces[var] = self.spaces[var_generic]

        for var_generic in list(set(self.prior_pams) & set(self.variable_compressed)):
            for var in self.variable_compressed[var_generic]:
                self.prior_pams[var] = self.prior_pams[var_generic]
                self.prior_kind[var] = self.prior_kind[var_generic]

        for var_generic in list(set(self.starts) & set(self.variable_compressed)):
            for var in self.variable_compressed[var_generic]:
                self.starts[var] = self.starts[var_generic]

    def has_jitter(self):
        for var in self.list_pams:
            if self.variable_expanded[var] == 'jitter':
                return True
        return False
