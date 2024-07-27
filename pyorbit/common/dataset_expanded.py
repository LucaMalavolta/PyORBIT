from pyorbit.subroutines.common import np, OrderedSet
from pyorbit.common.abstract_common import AbstractCommon
from pyorbit.common.dataset import Dataset
from pyorbit.model_definitions import datatype_definition

from numpy.lib.recfunctions import append_fields, drop_fields

try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle
import numpy as np


class DatasetExpanded(Dataset):

    def __init__(self,  model_name, kind, models):
        AbstractCommon.__init__(self, None)

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
            'jitter': ['Uniform', []],
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

        self.compute_plot = False



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
            input_array = pickle.load( open(input_file, "rb" ) )

        """ NOTE: string are only supported in ancillary data files
        """

        if self.ancillary:
            # Data ancillary has been already defined when reading the main files
            # we take the keywords from the file and add them to the existing
            for key, dict_vals in input_array.items():
                self.ancillary[key] = dict_vals
        else:
            self.ancillary = input_array.copy()

            for iname, name in enumerate(input_array.dtype.names):
                self.ancillary_str_index[name] = iname

        self.ancillary_str = self.ancillary.copy()


    def convert_dataset_from_file(self, input_file):

        """ Robust data reading, now encompassing the case of ancillary data
        columns embedded in the main file"""
        data_dictionary =  pickle.load( open(input_file, "rb" ) )

        """ By default, the dataset is copied into the ancillary file  """
        self.ancillary = data_dictionary.copy()

        for key in self.generic_list_pams:
            self.ancillary.pop(key, None)

            if key not in data_dictionary:
                data_dictionary[key] = np.asarray(-1)
            if type(data_dictionary[key]) == bool:
                data_dictionary[key] = np.asarray(data_dictionary[key]).astype(np.int64) -1
            if type(data_dictionary[key]) == int:
                data_dictionary[key] = np.asarray(data_dictionary[key])


        return data_dictionary

    def define_dataset_base(self, data_dictionary, update=False,
                            flag_shutdown_jitter=False):
        # Add a flag to save internally the input dataset

        if flag_shutdown_jitter:
            data_dictionary['jitter'] = -1

        if not self.models:
            data_dictionary['jitter'] = -1
            data_dictionary['offset'] = -1

        self.x = np.asarray(data_dictionary['bjd'], dtype=np.double)
        try:
            self.x1 = np.asarray(data_dictionary['2nd_axis'], dtype=np.double)
        except:
            self.x1 = None
        self.y = np.asarray(data_dictionary['data'], dtype=np.double)
        self.e = np.asarray(data_dictionary['errs'], dtype=np.double)

        self.n = np.size(self.y)
        self.n_shape = np.shape(self.y)

        if not update:
            if self.Tref is None:
                self.Tref = np.mean(self.x, dtype=np.double)

            """Default boundaries are defined according to the characteristic
               of the dataset. They must be large enough to allow for most of
               the anomalous situations
            """
            data_range = np.max(self.y) - np.min(self.y)

            self.generic_default_bounds = {'offset': [np.min(self.y) - 10.*data_range,
                                                        np.max(self.y) + 10.*data_range],
                                           'jitter': [0., 100 * np.max(self.e)]}

            self._setup_systematic_dictionaries('jitter', data_dictionary['jitter'])
            self._setup_systematic_dictionaries('offset', data_dictionary['offset'])

        self.x0 = self.x - self.Tref
        self._setup_systematic_mask('jitter', data_dictionary['jitter'])
        self._setup_systematic_mask('offset', data_dictionary['offset'])

        if np.amax(data_dictionary['subset']) > 0:
            sel = data_dictionary['subset'] >= -0.5
            self.submodel_minflag = np.int64(np.amin(data_dictionary['subset']))
            self.submodel_maxflag = np.int64(np.amax(data_dictionary['subset'])) + 1
            self.submodel_flag = np.int64(np.amax(data_dictionary['subset'])) + 1
            self.submodel_id = data_dictionary['subset']
        else:
            self.submodel_minflag = None
            self.submodel_maxflag = None
            self.submodel_flag = None


        self.model_reset()


