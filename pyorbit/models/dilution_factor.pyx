from ..classes.common import *
from abstract_common import *
from abstract_model import *

class CommonDilutionFactor(AbstractCommon):
    ''' This class must be created for each planet in the system
            model_name is the way the planet is identified
    '''
    model_class = 'dilution_factor'
    unitary_model = True

    list_pams = {
        'd_factor',  # diluition factor
    }

    default_bounds = {
        'd_factor': [0.0000, 1.0000]
    }

    """ Must be the same parameters as in list_pams, because priors are applied only to _physical_ parameters """
    default_priors = {
        'd_factor': ['Uniform', []]
    }

    default_spaces = {
        'd_factor': 'Linear'
    }

    default_fixed = {
        'd_factor': 0.0000
    }

    recenter_pams = {}

class DilutionFactor(AbstractModel):

    model_class = 'dilution_factor'
    unitary_model = True

    list_pams_common = {
        'd_factor',  # Diluition factor, expressed as all_other_stars / star_A ratio of flux
    }
    list_pams_dataset = {}

    default_bounds = {}
    default_spaces = {}
    default_priors = {}

    recenter_pams_dataset = {}

    unitary_array = {}

    #def setup_dataset(self, dataset, **kwargs):
    #    self.unitary_array[dataset.name_ref] = np.ones(dataset.n, dtype=np.double)

    def compute(self, variable_value, dataset, x0_input=None):

        """

        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """
        return np.asarray(variable_value['d_factor'])
