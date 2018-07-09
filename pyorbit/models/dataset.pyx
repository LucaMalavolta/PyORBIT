from ..classes.common import *
from abstract_common import AbstractCommon


class Dataset(AbstractCommon):

    def __init__(self, model_name, kind, models):

        super(self.__class__, self).__init__(None)

        self.kind = kind
        # model kind:  'RV', 'PHOT', 'ACT'...
        self.models = models
        self.name_ref = model_name

        self.dynamical = False
        self.planet_name = None

        self.generic_list_pams = {'jitter': 'U', 'offset': 'U', 'linear': 'U'}
        #self.generic_list_pams = {'jitter': 'LU', 'offset': 'U', 'linear': 'U'}

        self.generic_default_priors = {
            #'jitter': ['Jeffreys', []],
            'jitter': ['Uniform', []],
            'offset': ['Uniform', []],
            'linear': ['Uniform', []]}

        if self.kind == 'Tcent':
            self.generic_list_pams['jitter'] = 'LU'
            self.generic_default_priors['jitter'] = ['Jeffreys', []]

        self.variable_compressed = {}
        self.variable_expanded = {}

        self.model_class = 'dataset'
        self.list_pams = {}
        self.default_bounds = {}
        self.default_priors = {}
        self.recenter_pams = {}

        self.Tref = None
        self.model= None
        self.jitter = None
        self.mask = {}

    def convert_dataset_from_file(self, input_file):
        print 'Opening: ', input_file
        data = np.atleast_2d(np.loadtxt(input_file))

        data_input = np.zeros([np.size(data, axis=0), 6], dtype=np.double) - 1.
        data_input[:, :np.size(data, axis=1)] = data[:, :]
        return data_input

    def define_dataset_base(self, data_input, update=False, flag_shutdown_jitter=False):
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

            """Default boundaries are defined according to the characteristic of the dataset"""
            self.generic_default_bounds = {'offset': [np.min(self.y) - 1000., np.max(self.y) + 100.],
                                           'jitter': [np.min(self.e)/100., 100 * np.max(self.e)],
                                           'linear': [-1., 1.]}

            self.create_systematic_dictionaries('jitter', data_input[:, 3])
            self.create_systematic_dictionaries('offset', data_input[:, 4])
            self.create_systematic_dictionaries('linear', data_input[:, 5])

        self.x0 = self.x - self.Tref
        self.create_systematic_mask('jitter', data_input[:, 3])
        self.create_systematic_mask('offset', data_input[:, 4])
        self.create_systematic_mask('linear', data_input[:, 5])

        self.model = np.zeros(self.n, dtype=np.double)
        self.jitter = np.zeros(self.n, dtype=np.double)

    def create_systematic_dictionaries(self, var_generic, dataset_vals):
        n_sys = np.max(dataset_vals.astype(np.int64)) + 1
        self.variable_compressed[var_generic] = {}
        for ii in xrange(0, n_sys):
            var = var_generic + '_' + repr(ii)
            self.list_pams[var] = self.generic_list_pams[var_generic]
            self.default_bounds[var] = self.generic_default_bounds[var_generic]
            self.default_priors[var] = self.generic_default_priors[var_generic]
            self.variable_compressed[var_generic][var] = None
            self.variable_expanded[var] = var_generic

    def create_systematic_mask(self, var_generic, dataset_vals):
        n_sys = np.max(dataset_vals.astype(np.int64)) + 1
        for ii in xrange(0, n_sys):
            var = var_generic + '_' + repr(ii)
            self.mask[var] = np.zeros(self.n, dtype=bool)
            self.mask[var][(abs(dataset_vals - ii) < 0.1)] = True

    def delete_systematic_dictionaries_mask(self, var_generic):
        if var_generic not in self.variable_compressed: return
        for var in self.variable_compressed[var_generic]:
            self.list_pams.pop(var, None)
            self.default_bounds.pop(var, None)
            self.default_priors.pop(var, None)
            self.variable_expanded.pop(var, None)
            self.mask.pop(var)
        self.variable_compressed[var_generic] = {}

    def shutdown_jitter(self):
        self.delete_systematic_dictionaries_mask('jitter')

    def shutdown_offset(self):
        self.delete_systematic_dictionaries_mask('offset')

    def shutdown_linear(self):
        self.delete_systematic_dictionaries_mask('linear')

    def common_Tref(self, Tref_in):
        self.Tref = Tref_in
        self.x0 = self.x - self.Tref
        return

    def model_reset(self):
        self.model = np.zeros(self.n, dtype=np.double)
        self.jitter = np.zeros(self.n, dtype=np.double)
        return

    def compute(self, variable_value):
        for var in self.list_pams:
            if self.variable_expanded[var] == 'jitter':
                self.jitter[self.mask[var]] += variable_value[var]
            else:
                self.model[self.mask[var]] += variable_value[var]

    def model_logchi2(self):
        env = 1.0 / (self.e ** 2.0 + self.jitter ** 2.0)
        return -0.5 * (np.sum((self.y - self.model) ** 2 * env - np.log(env)))

    def update_priors_starts_bounds(self):

        for var_generic in list(set(self.prior_pams) & set(self.variable_compressed)):
            for var in self.variable_compressed[var_generic]:
                self.prior_pams[var] = self.prior_pams[var_generic]
                self.prior_kind[var] = self.prior_kind[var_generic]

        for var_generic in list(set(self.starts) & set(self.variable_compressed)):
            for var in self.variable_compressed[var_generic]:
                self.starts[var] = self.starts[var_generic]

        for var_generic in list(set(self.bounds) & set(self.variable_compressed)):
            for var in self.variable_compressed[var_generic]:
                self.bounds[var] = self.bounds[var_generic]

    def has_jitter(self):
        for var in self.list_pams:
            if self.variable_expanded[var] == 'jitter':
                return True
        return False
