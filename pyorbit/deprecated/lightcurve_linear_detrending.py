from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
import pyorbit.subroutines.constants as constants
from numpy.lib.recfunctions import append_fields

class LightcurveLinearDetrending(AbstractModel):

    default_common = 'lightcurve_detrending'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'lightcurve_linear_detrending'
        self.unitary_model = True
        self.normalization_model = False
        self.multiplicative_model = False
        self.time_independent_model = True


        self.list_pams_common = OrderedSet()
        self.list_pams_dataset = OrderedSet()

        self.ancillary_skip_first_column = True

        self.ancillary_names = {}

        self.baseline_value = 0.0000


    def initialize_model(self, mc, **kwargs):

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == self.default_common:
                self.common_lcd_ref = common_ref
                break

        if kwargs.get('include_zero_point', False) or self.normalization_model:
            self.starting_order = 0
            par_original = 'coeff_c0'
            par_addition = 'lcld_c0'
            mc.common_models[self.common_lcd_ref]._transfer_priors(mc, par_original, par_addition)
            self.list_pams_common.update([par_addition])

        """ WARNING: in order to use global coefficients for individual datasets
            it is compulsory that all the datasets are listed in the keyword sections
            and that the zero point is specified
        """
        for data_name, keywords in kwargs['ancillary_sets'].items():
            self.ancillary_names[data_name] = None
            par_original = 'coeff_linear'
            par_addition = 'lcld_coeff_'+data_name

            mc.common_models[self.common_lcd_ref]._transfer_priors(mc, par_original, par_addition)
            self.list_pams_common.update([par_addition])

            par_original = 'x_zero'
            par_addition = 'x_zero_' + data_name

            mc.common_models[self.common_lcd_ref].fix_list[par_addition] = np.asarray([keywords['x_zero'], 0.0000])
            mc.common_models[self.common_lcd_ref]._transfer_priors(mc, par_original, par_addition)
            self.list_pams_common.update([par_addition])

        if self.normalization_model:
            self.baseline_value = 1.0000

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        for data_name in self.ancillary_names:

            if data_name not in dataset.ancillary.dtype.names:
                dataset.ancillary = append_fields(dataset.ancillary, data_name, np.zeros(dataset.n))

    def compute(self, parameter_values, dataset, x0_input=None):

        if x0_input is None:

            trend = np.ones(dataset.n) * parameter_values.get('lcld_c0', self.baseline_value)

            for data_name in self.ancillary_names:
                trend += parameter_values['lcld_coeff_'+data_name] * (dataset.ancillary[data_name]-parameter_values['x_zero_'+data_name])


            return trend

        else:
            return x0_input * 0.0


class LocalLightcurveLinearDetrending(AbstractModel):

    default_common = 'lightcurve_detrending'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'local_lightcurve_linear_detrending'
        self.unitary_model = True
        self.normalization_model = False
        self.time_independent_model = True


        self.list_pams_common = OrderedSet()
        self.list_pams_dataset = OrderedSet()

        self.ancillary_skip_first_column = True

        self.lightcurve_ancillary = {}

        self.baseline_value = 1.0000

    def initialize_model(self, mc, **kwargs):

        if kwargs.get('ancillary_skip_first_column', False):
            self.ancillary_skip_first_column = kwargs['ancillary_skip_first_column']

        if self.normalization_model:
            self.baseline_value = 1.0000

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        self.lightcurve_ancillary[dataset.name_ref] = {}

        if self.ancillary_skip_first_column:
            skip_name = dataset.ancillary.dtype.names[0]
        else:
            skip_name = None

        if kwargs.get('include_zero_point', False) or self.normalization_model:
            par_original = 'coeff_c0'
            par_addition = 'lcld_c0'
            self.transfer_parameter_properties(mc, dataset, par_original, par_addition, dataset_pam=True)

        for data_name in dataset.ancillary.dtype.names:
            if data_name == skip_name: continue

            par_original = 'coeff_linear'
            par_addition = 'lcld_coeff_'+data_name

            self.transfer_parameter_properties(mc, dataset, par_original, par_addition, dataset_pam=True)


            par_original = 'x_zero'
            par_addition = 'x_zero_' + data_name

            try:
                x_zero = kwargs[data_name]['x_zero']
            except (KeyError, ValueError):
                x_zero = np.average(dataset.ancillary[data_name])

            self.fix_list[dataset.name_ref][par_addition] = np.asarray([x_zero, 0.0000])
            self.transfer_parameter_properties(mc, dataset, par_original, par_addition, dataset_pam=True)

            self.lightcurve_ancillary[dataset.name_ref][data_name] = dataset.ancillary[data_name] - x_zero


    def compute(self, parameter_values, dataset, x0_input=None):

        if x0_input is None:

            trend = np.ones(dataset.n) * parameter_values.get('lcld_c0', self.baseline_value)

            for data_name, data_vals in self.lightcurve_ancillary[dataset.name_ref].items():
                trend += parameter_values['lcld_coeff_'+data_name] * data_vals


            return trend

        else:
            return x0_input * 0.0

