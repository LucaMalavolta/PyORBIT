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


class Jitter(AbstractModel):

    model_class = 'common_jitter'
    jitter_model = True

    list_pams_common = {'jitter'}
    list_pams_dataset = {}

    recenter_pams_dataset = {}

    def __init__(self, *args, **kwargs):
        super(Jitter, self).__init__(*args, **kwargs)

        self.common_jitter_ref = None

    def initialize_model(self, mc, **kwargs):

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == 'common_jitter':
                self.common_jitter_ref = common_ref
                break

    def setup_dataset(self, mc, dataset, **kwargs):

        if not mc.common_models[self.common_jitter_ref].default_bounds:
            min_jitter = np.min(dataset.e) / 100.
            max_jitter = np.max(dataset.e) * 100.
        else:
            min_jitter = min(mc.common_models[self.common_jitter_ref].default_bounds['jitter'][0],
                             np.min(dataset.e) / 100.)
            max_jitter = max(mc.common_models[self.common_jitter_ref].default_bounds['jitter'][1],
                             np.max(dataset.e) * 100.)
        mc.common_models[self.common_jitter_ref].default_bounds['jitter'] = [min_jitter, max_jitter]
        dataset.shutdown_jitter()
        return

    def compute(self, variable_value, dataset, x0_input=None):
        return variable_value['jitter']



