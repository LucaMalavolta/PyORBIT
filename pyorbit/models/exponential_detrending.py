from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
import pyorbit.subroutines.constants as constants
from numpy.polynomial import polynomial
from numpy.lib.recfunctions import append_fields

class LightcurveDetrending(AbstractModel):

    default_common = 'lightcurve_detrending'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'exponential_detrending'
        self.unitary_model = False
        self.normalization_model = False
        self.multiplicative_model = True
        self.time_independent_model = True

        self.natural_base = False
        self.local_model = True

        self.list_pams_common = set()
        self.list_pams_dataset = set()

        self.ancillary_skip_first_column = True

        self.lightcurve_ancillary = {}
        self.ancillary_order = {}

        self.order = 1
        self.starting_order = 1
        self.threshold = 0.001
        self.baseline_value = 0.00000

        self.pams_order = {}
        self.pams_name  = {}

        self.initialized = False

    def initialize_model(self, mc, **kwargs):

        self.order = kwargs.get('order', self.order)
        self.local_model = kwargs.get('local_model', self.local_model)

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == self.default_common:
                self.common_lcpd_ref = common_ref
                break

        if self.multiplicative_model:
            self.baseline_value = 0.000

        if self.normalization_model:
            self.starting_order = 0

        if kwargs.get('include_zero_point', False):
            self.starting_order = 0

        """ The user may decide to include the 0th order anyway -
            be aware of correlations with dataset offset!"""
        if self.starting_order == 0 :
            par_original = 'coeff_c0'
            par_addition = 'expdet_c0'

            if self.local_model:
                self.transfer_parameter_properties(mc, dataset, par_original, par_addition, common_pam=True)
            else:
                self.transfer_parameter_properties(mc, dataset, par_original, par_addition, dataset_pam=True)

        """ WARNING: in order to use global coefficients for individual datasets
            it is compulsory that all the datasets are listed in the keyword sections
            and that the zero point is specified
        """
        for data_name, keywords in kwargs['ancillary_sets'].items():
            self.ancillary_order[data_name] = keywords['order']

            par_original = 'x_zero'
            par_addition = 'x_zero_' + data_name

            mc.common_models[self.common_lcpd_ref].fix_list[par_addition] = np.asarray([keywords['x_zero'], 0.0000])
            mc.common_models[self.common_lcpd_ref]._transfer_priors(mc, par_original, par_addition)
            self.list_pams_common.update([par_addition])

            for i_order in range(1, self.ancillary_order[data_name]+1):
                par_original = 'coeff_poly'
                par_addition = 'expdet_' + data_name + '_c'+repr(i_order)

                mc.common_models[self.common_lcpd_ref]._transfer_priors(mc, par_original, par_addition)
                self.list_pams_common.update([par_addition])

        if self.normalization_model:
            self.baseline_value = 0.00

        self.natural_base = kwargs.get('natural_base', self.natural_base )
        if self.natural_base:
            print(' Exponential detrending: using natural base')
        else:
            print(' Exponential detrending: using base 10')

        """BUG  incorrect procedure """

        self.initialized = False

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        """
        for i_sub in range(0, dataset.submodel_flag):
            par_original = 'Tc'
            par_subset = 'Tc_'+repr(i_sub)
            
            self.Tc_names[dataset.name_ref].append(par_subset)
            self.transfer_parameter_properties(mc, dataset, par_original, par_subset, dataset_pam=True)

            sub_dataset = dataset.x[(dataset.submodel_id == i_sub)]
            if kwargs[dataset.name_ref].get('boundaries', False):
                par_update = kwargs[dataset.name_ref]['boundaries'].get(
                    par_subset, [min(sub_dataset), max(sub_dataset)])
            elif kwargs.get('boundaries', False):
                par_update = kwargs['boundaries'].get(par_subset, [min(sub_dataset), max(sub_dataset)])
            else:
                par_update = [min(sub_dataset), max(sub_dataset)]

            self.bounds[dataset.name_ref].update({par_subset: par_update})
        """

        """ There is only one zeroth-order coefficient for the full trend model, otherwise
            parameter degenaries will arise
        """
        if self.starting_order == 0 and self.initialized == False:

            par_original = 'coeff_c0'
            par_addition = 'exp_c0'

            if self.local_model:
                self.transfer_parameter_properties(mc, dataset, par_original, par_addition, dataset_pam=True)
            else:
                self.transfer_parameter_properties(mc, dataset, par_original, par_addition, common_pam=True)

            """ No need to redo everything is the parameter has been aready initialized during the first instance"""
            self.initialized = True


        """ the name and the polynomial order for eac dataset and each set is preserved
            We assume that the same order are used for all the datasets that are using this model
            even if the coefficient for each dataset may vary
        """
    
        self.pams_order = {}
        self.pams_name[dataset.name_ref]  = {}

        for data_name in kwargs['detrending_variables']:

            self.pams_order[data_name] = {}


            if data_name not in dataset.ancillary.dtype.names:
                dataset.ancillary = append_fields(dataset.ancillary, data_name, np.zeros(dataset.n))
                self.pams_order[data_name] = 0
                continue


            nested_dict = kwargs['detrending_variables'][data_name]
            try:
                self.pams_order[data_name] = nested_dict.get('order', self.order)
            except:
                self.pams_order[data_name] = self.order

            for i_order in range(1, self.pams_order[data_name]+1):

                par_original = 'coeff_poly'
                par_addition = 'exp_' + data_name + '_c'+repr(i_order)

                if self.local_model:
                    self.transfer_parameter_properties(mc, dataset, par_original, par_addition, dataset_pam=True)
                else:
                    self.transfer_parameter_properties(mc, dataset, par_original, par_addition, common_pam=True)







    def compute(self, parameter_values, dataset, x0_input=None):

        if x0_input is None:

            trend = np.ones(dataset.n) * parameter_values.get('expdet_c0', self.baseline_value)

            for data_name, data_order in self.ancillary_order.items():

                coeff = np.zeros(data_order+1)
                for i_order in range(1, data_order+1):
                    par_addition = 'expdet_' + data_name + '_c'+repr(i_order)
                    coeff[i_order] = parameter_values[par_addition]

                trend += polynomial.polyval(dataset.ancillary[data_name]-parameter_values['x_zero_'+data_name], coeff)
    
            if self.natural_base:
                return np.exp(trend)
            else:
                return 10**(trend)


        else:
            return x0_input * 0.0



class LocalLightcurveDetrending(AbstractModel):

    default_common = 'lightcurve_detrending'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'local_lightcurve_detrending'
        self.unitary_model = False
        self.normalization_model = False
        self.multiplicative_model = True
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
            par_original = 'coeff_c0'
            par_addition = 'lcd_c0'

            self._subset_transfer_priors(mc, dataset, par_original, par_addition)

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

            par_original = 'x_zero'
            par_addition = 'x_zero_' + data_name

            self.fix_list[dataset.name_ref][par_addition] = np.asarray([x_zero, 0.0000])
            self._subset_transfer_priors(mc, dataset, par_original, par_addition)

            for i_order in range(1, self.ancillary_order[data_name]+1):

                par_original = 'coeff_poly'
                par_addition = 'lcd_' + data_name + '_c'+repr(i_order)

                self._subset_transfer_priors(mc, dataset, par_original, par_addition)

            self.lightcurve_ancillary[dataset.name_ref][data_name] = dataset.ancillary[data_name] - x_zero


    def compute(self, parameter_values, dataset, x0_input=None):

        if x0_input is None:

            trend = np.ones(dataset.n) * parameter_values.get('lcd_c0', self.baseline_value)

            for data_name, data_vals in self.lightcurve_ancillary[dataset.name_ref].items():

                coeff = np.zeros(self.ancillary_order[data_name]+1)
                for i_order in range(1, self.ancillary_order[data_name]+1):
                    par_addition = 'lcd_' + data_name + '_c'+repr(i_order)
                    coeff[i_order] = parameter_values[par_addition]

                trend += polynomial.polyval(data_vals, coeff)

            if self.natural_base:
                return np.exp(trend)
            else:
                return 10**(trend)

        else:
            return x0_input * 0.0

