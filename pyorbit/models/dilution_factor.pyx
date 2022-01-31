from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *


class DilutionFactor(AbstractModel):

    default_common = 'dilution_factor'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'dilution_factor'
        self.unitary_model = True

        self.list_pams_common = {
            'd_factor',  # Dilution factor, expressed as all_other_stars / star_A ratio of flux
        }

    def compute(self, variable_value, dataset, x0_input=None):
        """

        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """
        return np.asarray(variable_value['d_factor'])


class LocalDilutionFactor(AbstractModel):

    default_common = 'dilution_factor'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'local_dilution_factor'
        self.unitary_model = True

        self.list_pams_dataset = {
            'd_factor',  # Dilution factor, expressed as all_other_stars / star_A ratio of flux
        }

    def compute(self, variable_value, dataset, x0_input=None):
        """

        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """
        return np.asarray(variable_value['d_factor'])
