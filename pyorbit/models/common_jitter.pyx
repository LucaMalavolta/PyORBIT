from pyorbit.classes.common import *
from pyorbit.models.abstract_common import *
from pyorbit.models.abstract_model import *


class CommonJitter(AbstractCommon):
    ''' Common offset for datasets in different files
    '''

    model_class = 'common_jitter'

    list_pams = {
        'jitter'
    }

    default_bounds = {}
    default_spaces = {'jitter': 'Linear'}
    default_priors = {'jitter': ['Uniform', []]}
    default_fixed = {}

    recenter_pams = {}

    def common_initialization_with_dataset(self, dataset):
        if not self.default_bounds:
            min_jitter = np.min(dataset.e) / 100.
            max_jitter = np.max(dataset.e) * 100.
        else:
            min_jitter = min(self.default_bounds['jitter'][0], np.min(dataset.e) / 100.)
            max_jitter = max(self.default_bounds['jitter'][1], np.max(dataset.e) * 100.)
        self.default_bounds['jitter']= [min_jitter, max_jitter]
        dataset.shutdown_jitter()
        return


class Jitter(AbstractModel):

    model_class = 'common_jitter'
    jitter_model = True

    list_pams_common = {'jitter'}
    list_pams_dataset = {}

    recenter_pams_dataset = {}

    def compute(self, variable_value, dataset, x0_input=None):
        return variable_value['jitter']



