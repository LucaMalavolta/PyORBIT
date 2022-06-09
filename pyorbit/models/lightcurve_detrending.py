from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
import pyorbit.subroutines.constants as constants
from numpy.polynomial import polynomial
from numpy.lib.recfunctions import append_fields

class LightcurveDetrending(AbstractModel):

    default_common = 'lightcurve_detrending'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'lightcurve_detrending'
        self.unitary_model = True
        self.normalization_model = False
        self.multiplicative_model = False
        self.time_independent_model = True


        self.list_pams_common = set()
        self.list_pams_dataset = set()

        self.ancillary_skip_first_column = True

        self.lightcurve_ancillary = {}
        self.ancillary_order = {}

        self.order = 1
        self.starting_order = 1
        self.threshold = 0.001
        self.baseline_value = 0.00000

    def initialize_model(self, mc, **kwargs):

        self.order = kwargs.get('order', 1)

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == self.default_common:
                self.common_lcpd_ref = common_ref
                break

        if self.multiplicative_model:
            self.baseline_value = 1.0000

        if self.normalization_model:
            self.starting_order = 0

        if kwargs.get('include_zero_point', False):
            self.starting_order = 0

        """ The user may decide to include the 0th order anyway -
            be aware of correlations with dataset offset!"""
        if self.starting_order == 0 :
            var_original = 'coeff_c0'
            var_addition = 'lcd_c0'
            mc.common_models[self.common_lcpd_ref]._transfer_priors(mc, var_original, var_addition)
            self.list_pams_common.update([var_addition])


        """ WARNING: in order to use global coefficients for individual datasets
            it is compulsory that all the datasets are listed in the keyword sections
            and that the zero point is specified
        """
        for data_name, keywords in kwargs['ancillary_sets'].items():
            self.ancillary_order[data_name] = keywords['order']

            var_original = 'x_zero'
            var_addition = 'x_zero_' + data_name

            mc.common_models[self.common_lcpd_ref].fix_list[var_addition] = np.asarray([keywords['x_zero'], 0.0000])
            mc.common_models[self.common_lcpd_ref]._transfer_priors(mc, var_original, var_addition)
            self.list_pams_common.update([var_addition])

            for i_order in range(1, self.ancillary_order[data_name]+1):
                var_original = 'coeff_poly'
                var_addition = 'lcd_' + data_name + '_c'+repr(i_order)

                mc.common_models[self.common_lcpd_ref]._transfer_priors(mc, var_original, var_addition)
                self.list_pams_common.update([var_addition])

        if self.normalization_model:
            self.baseline_value = 1.0000

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        for data_name in self.ancillary_order:

            if data_name not in dataset.ancillary.dtype.names:
                dataset.ancillary = append_fields(dataset.ancillary, data_name, np.zeros(dataset.n))


    def compute(self, variable_value, dataset, x0_input=None):

        if x0_input is None:

            trend = np.ones(dataset.n) * variable_value.get('lcd_c0', self.baseline_value)

            for data_name, data_order in self.ancillary_order.items():

                coeff = np.zeros(data_order+1)
                for i_order in range(1, data_order+1):
                    var_addition = 'lcd_' + data_name + '_c'+repr(i_order)
                    coeff[i_order] = variable_value[var_addition]

                trend += polynomial.polyval(dataset.ancillary[data_name]-variable_value['x_zero_'+data_name], coeff)

            return trend

        else:
            return x0_input * 0.0



class LocalLightcurveDetrending(AbstractModel):

    default_common = 'lightcurve_detrending'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'local_lightcurve_detrending'
        self.unitary_model = True
        self.normalization_model = False
        self.multiplicative_model = False
        self.time_independent_model = True


        self.list_pams_common = set()
        self.list_pams_dataset = set()

        self.ancillary_skip_first_column = True

        self.lightcurve_ancillary = {}
        self.ancillary_order = {}

        self.order = 1
        self.starting_order = 1
        self.threshold = 0.001
        self.baseline_value = 0.0000

    def initialize_model(self, mc, **kwargs):

        self.order = kwargs.get('order', 1)

        """ The user may decide to include the 0th order anyway -
            be aware of correlations with dataset offset!"""


        if kwargs.get('ancillary_skip_first_column', False):
            self.ancillary_skip_first_column = kwargs['ancillary_skip_first_column']

        if self.multiplicative_model:
            self.baseline_value = 1.0000

        if self.normalization_model:
            self.starting_order = 0

        if kwargs.get('include_zero_point', False):
            self.starting_order = 0

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        self.lightcurve_ancillary[dataset.name_ref] = {}

        if self.ancillary_skip_first_column:
            skip_name = dataset.ancillary.dtype.names[0]
        else:
            skip_name = None

        if self.starting_order == 0 :
            var_original = 'coeff_c0'
            var_addition = 'lcd_c0'

            self._subset_transfer_priors(mc, dataset, var_original, var_addition)

        for data_name in dataset.ancillary.dtype.names:
            if data_name == skip_name: continue


            if kwargs.get('ancillary_sets', False):
                try:
                    x_zero = kwargs['ancillary_sets'][data_name]['x_zero']
                except (KeyError, ValueError):
                    x_zero = np.average(dataset.ancillary[data_name])

                try:
                    self.ancillary_order[data_name] = kwargs['ancillary_sets'][data_name]['order']
                except (KeyError, ValueError):
                    self.ancillary_order[data_name] = self.order

            elif kwargs.get(data_name, False):
                try:
                    x_zero = kwargs[data_name]['x_zero']
                except (KeyError, ValueError):
                    x_zero = np.average(dataset.ancillary[data_name])

                try:
                    self.ancillary_order[data_name] = kwargs[data_name]['order']
                except (KeyError, ValueError):
                    self.ancillary_order[data_name] = self.order

            elif kwargs.get(data_name + '_order', False):
                x_zero = np.average(dataset.ancillary[data_name])
                self.ancillary_order[data_name] = kwargs[data_name + '_order']
            else:
                x_zero = np.average(dataset.ancillary[data_name])
                self.ancillary_order[data_name] = self.order

            var_original = 'x_zero'
            var_addition = 'x_zero_' + data_name

            self.fix_list[dataset.name_ref][var_addition] = np.asarray([x_zero, 0.0000])
            self._subset_transfer_priors(mc, dataset, var_original, var_addition)

            for i_order in range(1, self.ancillary_order[data_name]+1):

                var_original = 'coeff_poly'
                var_addition = 'lcd_' + data_name + '_c'+repr(i_order)

                self._subset_transfer_priors(mc, dataset, var_original, var_addition)

            self.lightcurve_ancillary[dataset.name_ref][data_name] = dataset.ancillary[data_name] - x_zero


    def compute(self, variable_value, dataset, x0_input=None):

        if x0_input is None:

            trend = np.ones(dataset.n) * variable_value.get('lcd_c0', self.baseline_value)

            for data_name, data_vals in self.lightcurve_ancillary[dataset.name_ref].items():

                coeff = np.zeros(self.ancillary_order[data_name]+1)
                for i_order in range(1, self.ancillary_order[data_name]+1):
                    var_addition = 'lcd_' + data_name + '_c'+repr(i_order)
                    coeff[i_order] = variable_value[var_addition]

                trend += polynomial.polyval(data_vals, coeff)


            return trend

        else:
            return x0_input * 0.0

