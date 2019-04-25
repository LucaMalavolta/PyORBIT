from pyorbit.classes.common import *
from pyorbit.models.abstract_common import *
from pyorbit.models.abstract_model import *


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

    recenter_pams_dataset = {}

    # TODO: check if this dictionary is still required
    unitary_array = {}

    def compute(self, variable_value, dataset, x0_input=None):
        """

        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """
        return np.asarray(variable_value['d_factor'])


class LocalDilutionFactor(AbstractModel):
    model_class = 'dilution_factor'
    unitary_model = True

    list_pams_common = {
        'd_factor',  # Diluition factor, expressed as all_other_stars / star_A ratio of flux
    }
    list_pams_dataset = {'d_factor'}

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

    recenter_pams_dataset = {}

    # TODO: check if this dictionary is still required
    unitary_array = {}

    def compute(self, variable_value, dataset, x0_input=None):
        """

        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """
        return np.asarray(variable_value['d_factor'])
