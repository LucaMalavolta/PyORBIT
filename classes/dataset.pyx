from common import *
from abstract_common import AbstractCommon


class Dataset(AbstractCommon):

    def __init__(self, kind, input_file, models):

        super(self.__class__, self).__init__(None)

        self.kind = kind
        # model kind:  'RV', 'PHOT', 'ACT'...
        self.models = models
        self.name_ref = input_file

        self.dynamical = False
        self.planet_name = 'aaaaaaa'

        self.generic_list_pams = {'jitter': 'LU', 'offset': 'U', 'linear': 'U'}

        self.variable_compressed = {}
        self.variable_expanded = {}

        self.model_class = 'dataset'
        self.list_pams = {}
        self.default_bounds = {}
        self.recenter_pams = {}

        print 'Opening: ', input_file
        self.data = np.atleast_2d(np.loadtxt(input_file))

        n_cols = np.size(self.data, axis=1)

        if self.kind == 'Tcent':
            """ Special input reading from T0 files """
            self.n_transit = np.asarray(self.data[:, 0], dtype=np.int16)
            self.x = np.asarray(self.data[:, 1], dtype=np.double)
            self.e = np.asarray(self.data[:, 2], dtype=np.double)
            """ copy of self.y added for consistency with the rest of the code
            """
            self.y = self.x
        else:
            self.x = np.asarray(self.data[:, 0], dtype=np.double)
            self.y = np.asarray(self.data[:, 1], dtype=np.double)
            self.e = np.asarray(self.data[:, 2], dtype=np.double)

        self.n = np.size(self.x)

        self.Tref = np.mean(self.x, dtype=np.double)
        self.x0 = self.x - self.Tref

        """Default boundaries are defined according to the characteristic of the dataset"""
        self.generic_default_bounds = {'offset': [np.min(self.y) - 50., np.max(self.y) + 50.],
                               'jitter': [0.0001, 20 * np.max(self.e)],
                               'linear': [-1., 1.]}

        self.sys = {}
        self.mask = {}

        if self.models:
            if n_cols > 3:
                self.create_list_pams('jitter', 3)
            if n_cols > 4:
                self.create_list_pams('offset', 4)
            if n_cols > 5:
                self.create_list_pams('linear', 5)

        self.model = np.zeros(self.n, dtype=np.double)
        self.jitter = np.zeros(self.n, dtype=np.double)

    def create_list_pams(self, var_generic, data_id):
        self.sys[var_generic] = np.asarray(self.data[:, data_id], dtype=np.double)
        n_sys = np.max(self.sys[var_generic].astype(np.int64)) + 1
        self.variable_compressed[var_generic] = {}
        for ii in xrange(0, n_sys):
            var = var_generic + '_' + repr(ii)
            self.list_pams[var] = self.generic_list_pams[var_generic]
            self.default_bounds[var] = self.generic_default_bounds[var_generic]
            self.variable_compressed[var_generic] = var
            self.variable_expanded[var] = var_generic
            self.mask[var] = np.zeros(self.n, dtype=bool)
            self.mask[var][(abs(self.sys[var_generic] - ii) < 0.1)] = True

    def common_Tref(self, Tref_in):
        self.Tref = Tref_in
        self.x0 = self.x - self.Tref
        return

    def model_reset(self):
        self.model[:] = 0.0
        self.jitter[:] = 0.0
        return

    def compute(self, variable_value):
        for var in self.list_pams:
            if self.variable_expanded[var] is 'jitter':
                self.jitter[self.mask[var]] += variable_value[var]
            else:
                self.model[self.mask[var]] += variable_value[var]

    def model_logchi2(self):
        env = 1.0 / (self.e ** 2.0 + self.jitter ** 2.0)
        return -0.5 * (np.sum((self.y - self.model) ** 2 * env - np.log(env)))

    def update_priors_starts_bounds(self):

        for var_generic in list(set(self.prior_pams) & set(self.variable_compressed)):
            for var in self.variable_compressed[var_generic]:
                self.prior_pams[var] =  self.prior_pams[var_generic]
                self.prior_kind[var] =  self.prior_kind[var_generic]

        for var_generic in list(set(self.starts) & set(self.variable_compressed)):
            for var in self.variable_compressed[var_generic]:
                self.starts[var] =  self.starts[var_generic]

        for var_generic in list(set(self.bounds) & set(self.variable_compressed)):
            for var in self.variable_compressed[var_generic]:
                self.bounds[var] =  self.bounds[var_generic]

    #def shutdown_jitter(self):
    #    self.n_sys['jitter'] = 0

