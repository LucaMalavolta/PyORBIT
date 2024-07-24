from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *

class NormalizationFactor(AbstractModel):

    default_common = 'normalization_factor'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'normalization_factor'
        self.unitary_model = False
        self.normalization_model = True

        self.list_pams_common = OrderedSet([
            'n_factor',  # Normalization factor, expressed as all_other_stars / star_A ratio of flux
        ])

    def compute(self, parameter_values, dataset, x0_input=None):

        """

        :param parameter_values:
        :param dataset:
        :param x0_input:
        :return:
        """
        return np.asarray(parameter_values['n_factor'])

class LocalNormalizationFactor(AbstractModel):

    default_common = 'normalization_factor'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'normalization_factor'
        self.unitary_model = False
        self.normalization_model = True

        self.list_pams_dataset = OrderedSet([
            'n_factor',  # Normalization factor, expressed as all_other_stars / star_A ratio of flux
        ])

    def compute(self, parameter_values, dataset, x0_input=None):

        """

        :param parameter_values:
        :param dataset:
        :param x0_input:
        :return:
        """
        return np.asarray(parameter_values['n_factor'])


class SubsetNormalizationFactor(AbstractModel):

    default_common = 'normalization_factor'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'subset_normalization_factor'
        self.unitary_model = False
        self.normalization_model = True

        self.list_pams_dataset = OrderedSet()

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        if not dataset.submodel_flag:
            return

        for i_sub in range(0,dataset.submodel_flag):
            par_original = 'n_factor'
            par_subset = 'n_factor_sub'+repr(i_sub)
            self.transfer_parameter_properties(mc, dataset, par_original, par_addition, dataset_pam=True)

    def compute(self, parameter_values, dataset, x0_input=None):

        """

        :param parameter_values:
        :param dataset:
        :param x0_input:
        :return:
        """

        if x0_input is None:
            y_output = np.zeros(dataset.n)
            one_array = np.ones(dataset.n)

            for i_sub in range(0,dataset.submodel_flag):
                par = 'n_factor_sub'+repr(i_sub)
                sel_data = (dataset.submodel_id==i_sub)
                y_output[sel_data] = one_array[sel_data] * parameter_values[par]

            return y_output
        else:
            return x0_input*1.