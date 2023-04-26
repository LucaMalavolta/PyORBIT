from pyorbit.subroutines.common import np
from pyorbit.common.abstract_common import AbstractCommon
from pyorbit.common.dataset import Dataset
from pyorbit.model_definitions import datatype_definition

from numpy.lib.recfunctions import append_fields, drop_fields

try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle
import numpy as np


class DatasetExpanded(AbstractCommon, Dataset):

    def __init__(self, *args, **kwargs):
        # this calls all constructors up to AbstractModel
        super().__init__(*args, **kwargs)
        super(AbstractCommon, self).__init__(*args, **kwargs)

        self.compute_plot = False

    def convert_dataset_from_file(self, input_file):

        """ Robust data reading, now encompassing the case of ancillary data
        columns embedded in the main file"""
        data_dictionary =  pickle.load(input_file, "rb")

        for key in self.generic_list_pams:
            if not data_dictionary.get(key, False):
                data_dictionary[key] = np.asarray(-1)
            if type(data_dictionary[key]) == bool:
                data_dictionary[key] = np.asarray(data_dictionary[key]).astype(np.int64) -1

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
        self.x1 = np.asarray(data_dictionary['2nd_axis'], dtype=np.double)
        self.y = np.asarray(data_dictionary['data'], dtype=np.double)
        self.e = np.asarray(data_dictionary['errs'], dtype=np.double)

        self.n = np.size(self.y)
        self.n_shape = np.shape(self.x)

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
            self.submodel_flag = np.int(np.amax(data_dictionary['subset'])) + 1
            self.submodel_id = data_dictionary['subset']
        else:
            self.submodel_flag = None

        self.model_reset()


