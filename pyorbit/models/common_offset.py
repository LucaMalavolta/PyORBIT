from pyorbit.models.abstract_model import *


class Offset(AbstractModel):

    default_common = 'common_offset'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'common_offset'
        self.systematic_model = True

        self.list_pams_common = OrderedSet(['offset'])

        self.common_offset_ref = None

    def initialize_model(self, mc, **kwargs):

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == 'common_offset':
                self.common_offset_ref = common_ref
                break

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        try:
            min_offset = min(mc.common_models[self.common_offset_ref].default_bounds['offset'][0],
                             np.min(dataset.y) - 100.)
            max_offset = max(mc.common_models[self.common_offset_ref].default_bounds['offset'][1],
                             np.max(dataset.y) + 100.)
        except KeyError:
            min_offset = np.min(dataset.y) - 100.
            max_offset = np.max(dataset.y) + 100.

        mc.common_models[self.common_offset_ref].default_bounds['offset'] = [
            min_offset, max_offset]
        dataset.shutdown_offset()
        return

    def compute(self, parameter_values, dataset, x0_input=None):
        return parameter_values['offset']
