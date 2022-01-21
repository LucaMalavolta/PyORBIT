from pyorbit.models.abstract_model import *


class Offset(AbstractModel):

    model_class = 'common_offset'
    systematic_model = True

    list_pams_common = {'offset': 'U'}
    list_pams_dataset = set()

    def __init__(self, *args, **kwargs):
        super(Offset, self).__init__(*args, **kwargs)

        self.common_offset_ref = None

    def initialize_model(self, mc, **kwargs):

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == 'common_offset':
                self.common_offset_ref = common_ref
                break

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        if not mc.common_models[self.common_offset_ref].default_bounds:
            min_offset = np.min(dataset.y) - 100.
            max_offset = np.max(dataset.y) + 100.
        else:
            min_offset = min(mc.common_models[self.common_offset_ref].default_bounds['offset'][0],
                             np.min(dataset.y) - 100.)
            max_offset = max(mc.common_models[self.common_offset_ref].default_bounds['offset'][1],
                             np.max(dataset.y) + 100.)
        mc.common_models[self.common_offset_ref].default_bounds['offset'] = [
            min_offset, max_offset]
        dataset.shutdown_offset()
        return

    def compute(self, variable_value, dataset, x0_input=None):
        return variable_value['offset']
