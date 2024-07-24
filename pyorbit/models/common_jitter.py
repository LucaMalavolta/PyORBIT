from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *


class Jitter(AbstractModel):

    default_common = 'common_jitter'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'common_jitter'
        self.jitter_model = True

        self.list_pams_common = OrderedSet(['jitter'])

        self.common_jitter_ref = None

    def initialize_model(self, mc, **kwargs):

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == 'common_jitter':
                self.common_jitter_ref = common_ref
                break

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        try:
            min_jitter = min(mc.common_models[self.common_jitter_ref].default_bounds['jitter'][0],
                             np.min(dataset.e) / 100.)
            max_jitter = max(mc.common_models[self.common_jitter_ref].default_bounds['jitter'][1],
                             np.max(dataset.e) * 100.)
        except KeyError:
            min_jitter = np.min(dataset.e) / 100.
            max_jitter = np.max(dataset.e) * 100.

        mc.common_models[self.common_jitter_ref].default_bounds['jitter'] = [
            min_jitter, max_jitter]
        dataset.shutdown_jitter()
        return

    def compute(self, parameter_values, dataset, x0_input=None):
        return parameter_values['jitter']
