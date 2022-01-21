from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *

class NormalizationFactor(AbstractModel):

    model_class = 'normalization_factor'
    unitary_model = False
    normalization_model = True

    list_pams_common = {
        'n_factor',  # Diluition factor, expressed as all_other_stars / star_A ratio of flux
        }
    list_pams_dataset = set()

    default_bounds = {}
    default_spaces = {}
    default_priors = {}

    def compute(self, variable_value, dataset, x0_input=None):

        """

        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """
        return np.asarray(variable_value['n_factor'])
