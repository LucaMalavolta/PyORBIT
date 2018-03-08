from abstract_common import *
from abstract_model import *


class CommonOffset(AbstractCommon):
    ''' Common offset for datasets in different files
    '''

    model_class = 'common_offset'

    list_pams = {
        'offset': 'U',  # order 1
    }

    default_bounds = {}

    recenter_pams = {}

    def common_initialization_with_dataset(self, dataset):
        if not self.default_bounds:
            min_offset = np.min(dataset.y) - 100
            max_offset = np.max(dataset.y) + 100.
        else:
            min_offset = min(self.default_bounds['offset'][0], np.min(dataset.e) - 100.)
            max_offset = max(self.default_bounds['offset'][1], np.max(dataset.e) + 100.)
        self.default_bounds['offset'] = [min_offset, max_offset]
        dataset.shutdown_offset()
        return


class Offset(AbstractModel):

    model_class = 'common_offset'
    list_pams_common = {'offset': 'U'}
    list_pams_dataset = {}

    recenter_pams_dataset = {}

    def compute(self, variable_value, dataset, x0_input=None):
        return variable_value['offset']



