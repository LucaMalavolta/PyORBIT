from pyorbit.classes.common import *
from pyorbit.models.abstract_common import *
from pyorbit.models.abstract_model import *


class CommonNormalizationFactor(AbstractCommon):
    ''' This class must be created for each planet in the system
            model_name is the way the planet is identified
    '''
    model_class = 'normalization_factor'
    unitary_model = False
    normalization_model = True

    list_pams = {
        'n_factor',  # normalization factor
    }

    default_bounds = {
        'n_factor': [0.0000010, 100000.0000]
    }

    """ Must be the same parameters as in list_pams, because priors are applied only to _physical_ parameters """
    default_priors = {
        'n_factor': ['Uniform', []]
    }

    default_spaces = {
        'n_factor': 'Logarithmic'
    }

    default_fixed = {
        'n_factor': 1.0000
    }

    recenter_pams = {}


class NormalizationFactor(AbstractModel):

    model_class = 'normalization_factor'
    unitary_model = False
    normalization_model = True

    list_pams_common = {
        'n_factor',  # Diluition factor, expressed as all_other_stars / star_A ratio of flux
    }
    list_pams_dataset = {}

    default_bounds = {}
    default_spaces = {}
    default_priors = {}

    recenter_pams_dataset = {}


    #def setup_dataset(self, dataset, **kwargs):
    #    self.unitary_array[dataset.name_ref] = np.ones(dataset.n, dtype=np.double)

    def compute(self, variable_value, dataset, x0_input=None):

        """

        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """
        return np.asarray(variable_value['n_factor'])
