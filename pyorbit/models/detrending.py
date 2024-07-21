from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
import pyorbit.subroutines.constants as constants
from numpy.polynomial import polynomial
from numpy.lib.recfunctions import append_fields
from scipy.interpolate import interp1d

class Detrending(AbstractModel):

    default_common = 'detrending'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        import warnings
        warnings.filterwarnings("ignore", category=RuntimeWarning)

        self.model_class = 'detrending'
        self.unitary_model = False
        self.normalization_model = False

        self.exponential_detrending = False
        self.natural_base = False
        self.local_model = True
        self.use_median_xzero = True

        self.list_pams_common = OrderedSet()
        self.list_pams_dataset = OrderedSet()

        self.order = 1
        self.x_zero = 0.0

        self.starting_order = 1
        self.baseline_value = 0.000

        self.pams_order = {}

        self.interpolated = {}

    def initialize_model(self, mc, **kwargs):

        self.order = kwargs.get('order', self.order)
        self.x_zero = kwargs.get('x_zero', self.x_zero)

        self.local_model = kwargs.get('local_model', self.local_model)

        self.common_model = ~self.local_model
        for keyword in ['use_common_parameters', 'common_parameters', 'use_common_model', 'common_model']:
            self.common_model = kwargs.get(keyword, self.common_model)
        self.local_model = ~self.common_model

        self.use_median_xzero = kwargs.get('use_median_xzero', self.use_median_xzero)
        self.exponential_detrending = kwargs.get('exponential_detrending', self.exponential_detrending)

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == self.default_common:
                self.common_detrending = common_ref
                break

        """ If the polynomial is used as normalization factor, the first order must be included"""
        if self.normalization_model:
            self.starting_order = 0
            if self.exponential_detrending:
                self.baseline_value = kwargs.get('baseline_value', 0.0)
            else:
                self.baseline_value = kwargs.get('baseline_value', 1.0)
        else:
            self.starting_order = 1
            self.baseline_value = kwargs.get('baseline_value', 0.0)

        """ The user may decide to include the 0th order anyway -
            be aware of correlations with dataset offset!"""
        if kwargs.get('include_zero_point', self.include_zero_point):
            self.starting_order = 0

        if kwargs.get('exclude_zero_point', self.exclude_zero_point):
            self.starting_order = 1

        self.natural_base = kwargs.get('natural_base', self.natural_base )

        print()
        print('     model:              ', self.model_name)
        print('     unitary_model:      ', self.unitary_model)
        print('     normalization_model:', self.normalization_model)

        print('     local_model:        ', self.local_model)
        print('     use_median_xzero:   ', self.use_median_xzero)

        print('     order:              ', self.order)
        print('     starting_order:     ', self.starting_order)
        print('     baseline_value:     ', self.baseline_value)

        if self.exponential_detrending:
            if self.natural_base:
                print('     Exponential detrending: using natural base')
            else:
                print('     Exponential detrending: using base 10')
        else:
            print('     Polynomial detrending')
        print()


    def initialize_model_dataset(self, mc, dataset, **kwargs):

        """ There is only one zeroth-order coefficient for the full trend model, otherwise
            parameter degeneracies will arise
            Be aware of possible correlations with other parameters!
        """
        if self.starting_order == 0:
            par_original = 'det_c0'
            par_addition = 'det_c0'
            if self.local_model:
                self.transfer_parameter_properties(mc, dataset, par_original, par_addition, keywords=kwargs, dataset_pam=True)
            else:
                self.transfer_parameter_properties(mc, dataset, par_original, par_addition, keywords=kwargs, common_pam=True)

        """ the name and the polynomial order for each dataset and each set is preserved
            We assume that the same order are used for all the datasets that are using this model
            even if the coefficient for each dataset may vary
        """

        self.pams_order = {}

        """ For plotting purposes only """
        self.interpolated[dataset.name_ref] = {}


        for data_name in kwargs['detrending_variables']:

            self.pams_order[data_name] = {}

            if data_name not in dataset.ancillary.dtype.names:
                dataset.ancillary = append_fields(dataset.ancillary, data_name, np.zeros(dataset.n))
                self.pams_order[data_name] = 0
                continue

            if self.use_median_xzero:
                x_zero = np.median(dataset.ancillary[data_name])
            else:
                x_zero = self.x_zero

            try:
                nested_dict = kwargs['detrending_variables'][data_name]
                self.pams_order[data_name] = nested_dict.get('order', self.order)
                x_zero =  nested_dict.get('x_zero', x_zero)
            except AttributeError:
                self.pams_order[data_name] = kwargs['detrending_variables'][data_name]
            except TypeError:
                self.pams_order[data_name] = self.order


            for i_order in range(1, self.pams_order[data_name]+1):

                par_original = 'det_c'+repr(i_order)
                par_addition = 'det_' + data_name + '_c'+repr(i_order)

                if self.local_model:
                    self.transfer_parameter_properties(mc, dataset, par_original, par_addition, keywords=kwargs, dataset_pam=True)
                else:
                    self.transfer_parameter_properties(mc, dataset, par_original, par_addition, keywords=kwargs, common_pam=True)

            par_original = 'x_zero'
            par_addition = 'x_zero_' + data_name

            if self.local_model:
                self.fix_list[dataset.name_ref][par_addition] = np.asarray([x_zero, 0.0])
                self.transfer_parameter_properties(mc, dataset, par_original, par_addition, keywords=kwargs, dataset_pam=True)
            else:
                mc.common_models[self.common_detrending].fix_list[par_addition] = np.asarray([self.x_zero, 0.0000])
                self.transfer_parameter_properties(mc, dataset, par_original, par_addition, keywords=kwargs, common_pam=True)

            self.interpolated[dataset.name_ref][data_name]=interp1d(
            dataset.x0,
            dataset.ancillary[data_name],
            bounds_error=False,
            fill_value=(np.amin(dataset.ancillary[data_name]), np.amax(dataset.ancillary[data_name])))


    def compute(self, parameter_values, dataset, x0_input=None):

        if x0_input is None:

            trend = np.ones(dataset.n) * parameter_values.get('det_c0', self.baseline_value)

            for data_name, data_order in self.pams_order.items():

                coeff = np.zeros(data_order+1)
                for i_order in range(1, data_order+1):
                    par_addition = 'det_' + data_name + '_c'+repr(i_order)
                    coeff[i_order] = parameter_values[par_addition]

                trend += polynomial.polyval(dataset.ancillary[data_name]-parameter_values['x_zero_'+data_name], coeff)
                #if dataset.name_ref == 'GROND_r':
                #    print(self.model_name, data_name, trend[20:23], coeff)
                #    print()
            if self.exponential_detrending:
                if self.natural_base:
                    return np.exp(trend)
                else:
                    return 10**(trend)
            else:
                return trend

        else:

            trend = np.ones_like(x0_input) * parameter_values.get('det_c0', self.baseline_value)

            for data_name, data_order in self.pams_order.items():

                coeff = np.zeros(data_order+1)
                for i_order in range(1, data_order+1):
                    par_addition = 'det_' + data_name + '_c'+repr(i_order)
                    coeff[i_order] = parameter_values[par_addition]

                trend += polynomial.polyval(self.interpolated[dataset.name_ref][data_name](x0_input)-parameter_values['x_zero_'+data_name], coeff)

            if self.exponential_detrending:
                if self.natural_base:
                    return np.exp(trend)
                else:
                    return 10**(trend)
            else:
                return trend


class PolynomialDetrending(Detrending):

    def __init__(self, *args, **kwargs):
        AbstractModel.__init__(self, *args, **kwargs)
        Detrending.__init__(self, *args, **kwargs)
        self.exponential_detrending = False

class ExponentialDetrending(Detrending):

    def __init__(self, *args, **kwargs):
        AbstractModel.__init__(self, *args, **kwargs)
        Detrending.__init__(self, *args, **kwargs)
        self.exponential_detrending = True