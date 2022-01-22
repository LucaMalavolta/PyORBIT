from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *

class NormalizationFactor(AbstractModel):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'normalization_factor'
        self.unitary_model = False
        self.normalization_model = True

        self.list_pams_common = {
            'n_factor',  # Normalization factor, expressed as all_other_stars / star_A ratio of flux
            }

    def compute(self, variable_value, dataset, x0_input=None):

        """

        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """
        return np.asarray(variable_value['n_factor'])

class LocalNormalizationFactor(AbstractModel):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'normalization_factor'
        self.unitary_model = False
        self.normalization_model = True

        self.list_pams_dataset = {
            'n_factor',  # Normalization factor, expressed as all_other_stars / star_A ratio of flux
            }

    def compute(self, variable_value, dataset, x0_input=None):

        """

        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """
        return np.asarray(variable_value['n_factor'])
