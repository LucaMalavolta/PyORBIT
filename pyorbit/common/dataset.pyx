from pyorbit.subroutines.common import np
from pyorbit.common.abstract_common import AbstractCommon
from pyorbit.model_definitions import datatype_definition



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

        self.generic_list_pams = {'jitter', 'offset'}

        self.generic_default_priors = {
            'jitter': ['Uniform', []],
            'offset': ['Uniform', []]}

        self.generic_default_spaces = {
            'jitter': 'Linear',
            'offset': 'Linear'}

        if self.kind == 'Tcent':
            self.generic_default_spaces['jitter'] = 'Logarithmic'
            self.generic_default_priors['jitter'] = ['Uniform', []]

        if self.kind == 'CCFs':
            self.Tref = 0.000

        self.variable_compressed = {}
        self.variable_expanded = {}

        self.model_class = 'dataset'
        self.list_pams = set()
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

    def convert_ancillary_from_file(self, input_file):
        """ Function to read an ancillary file containing a dataset required to
        model some effect (usually, instrumental) but not object of modelling
        itself (i.e., the dataset does not enter in the log-likelihhod)
        Added in PyORBIT 9.0

        Args:
            input_file (string): name of the file, absolute or relative path
        """
        # read ancillary data from txt file
        self.ancillary = np.genfromtxt(input_file, names=True)

    def convert_dataset_from_file(self, input_file):
        data = np.atleast_2d(np.loadtxt(input_file))

        data_input = np.zeros([np.size(data, axis=0), 6], dtype=np.double) - 1.
        data_input[:, :np.size(data, axis=1)] = data[:, :]
        return data_input

    def define_dataset_base(self, data_input, update=False,
                            flag_shutdown_jitter=False):
        # Add a flag to save internally the input dataset

        if flag_shutdown_jitter:
            data_input[:, 4] = -1

        if not self.models:
            data_input[:, 3:] = -1

        if self.kind == 'Tcent':
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

        if not update:
            if self.Tref is None:
                self.Tref = np.mean(self.x, dtype=np.double)

            """Default boundaries are defined according to the characteristic
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
            self.submodel_flag = np.int(np.amax(data_input[:, 5])) + 1
            self.submodel_id = data_input[:, 5]
        else:
            self.submodel_flag = None

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
        for ii in range(0, n_sys):
            var = var_generic + '_' + repr(ii)
            self.mask[var] = np.zeros(self.n, dtype=bool)
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

    def shutdown_offset(self):
        self._delete_systematic_dictionaries_mask('offset')

    def common_Tref(self, Tref_in):
        if self.kind == 'CCFs':
            self.Tref = 0.000
        else:
            self.Tref = Tref_in
        self.x0 = self.x - self.Tref
        return

    def model_reset(self):
        self.residuals = None
        self.model = None
        self.additive_model = np.zeros(self.n, dtype=np.double)
        self.unitary_model = np.zeros(self.n, dtype=np.double)
        self.normalization_model = None
        self.external_model = np.zeros(self.n, dtype=np.double)
        self.jitter = np.zeros(self.n, dtype=np.double)
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
