from pyorbit.subroutines.common import np, OrderedSet
from pyorbit.common.abstract_common import AbstractCommon
from pyorbit.common.dataset import Dataset
from pyorbit.datatype_definitions import datatype_definition


class DatasetExternal(Dataset):
    """DatasetExternal Special class to handle external datasets. 
    It just loads the dataset required by the external code and nothing else. All flags must be disabled 

    :param Dataset: _description_
    :type Dataset: _type_
    """

    def __init__(self, model_name, kind, models):
        AbstractCommon.__init__(self, None)

        self.kind = kind
        for kind_name, kind_list in datatype_definition.items():
            if kind in kind_list:
                self.kind = kind_name
                break

        # model kind:  'RV', 'PHOT', 'ACT'...
        self.models = models
        self.name_ref = model_name

        self.dynamical = False
        self.planet_name = None

        self.generic_list_pams = OrderedSet()
        """ No parameters associated to this dataset"""

        self.generic_default_priors = {}
        self.generic_default_spaces = {}
        self.generic_default_fixed = {}

        """ transit_duration files have four data columns before flags and ancillary data
            For now, all the other datasets foresee three data columns before flags and ancillary data
        """

        self.model_class = 'dataset'
        self.list_pams = OrderedSet()

        self.Tref = None
        self.residuals = None
        self.residuals_for_regression = None
        self.model = None
        # this model is compute externally and passed to the compute subroutine of the model
        self.external_model = None
        self.buffer_model = None
        # this model is used to compute the additive part of the model
        self.additive_model = None
        self.unitary_model = None
        self.normalization_model = None
        self.jitter = None
        self.mask = {}

        self.ancillary = None

        self.compute_plot = False
        self.n = 0


    def append_ancillary(self, input_file, input_array=False, input_array_str=False):
        """ This method must do nothing, as ancillary data are not supported in orbitize datasets"""
        pass 


    def convert_dataset_from_file(self, input_file):
        """ This method must do nothing, as data files are read by orbitize internal functions"""

        self.input_file = input_file
        

        return None

    def define_dataset_base(self, data_input, update=False,
                            flag_shutdown_jitter=False):
        # Add a flag to save internally the input dataset
        """This method does nothing as there no dataset parameters for orbitize."""
        self.model_reset()

    def _setup_systematic_dictionaries(self, var_generic, dataset_vals):
        """This method does nothing as there no dataset parameters for orbitize."""
        pass

    def _setup_systematic_mask(self, var_generic, dataset_vals):
        """This method does nothing as there no dataset parameters for orbitize."""
        pass

    def _delete_systematic_dictionaries_mask(self, var_generic):
        """This method does nothing as there no dataset parameters for orbitize."""
        return 

    def shutdown_jitter(self):
        """This method does nothing as jitter is dealt by orbitize."""
        pass

    def shutdown_offset(self):
        """This method does nothing as offset is dealt by orbitize."""
        pass

    def common_Tref(self, Tref_in):
        self.Tref = Tref_in
        return

    def model_reset(self):
        self.log_likelihood = 0.000
        return

    def compute(self, variable_value):
        """This method does nothing as internal parameters are dealt by orbitize."""

    def compute_model(self):
        """This method does nothing as modelling is dealt by orbitize."""

    def compute_model_from_arbitrary_datasets(self, additive_model, unitary_model, normalization_model, external_model):
        """This method does nothing as modelling is dealt by orbitize."""
        return None

    def compute_residuals(self):
        """This method does nothing as modelling is dealt by orbitize."""

    def model_logchi2(self):
        return self.log_likelihood

    def update_bounds_spaces_priors_starts(self):
        """This method does nothing as there are no internarl parameters in orbitize datasets."""

    def has_jitter(self):
        return False
