from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
from numpy.polynomial import polynomial

class PolynomialTrend(AbstractModel):

    default_common = 'polynomial_trend'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'polynomial_trend'

        self.list_pams_common = {'x_zero'}

        self.order = 1
        self.starting_order = 1

        """
        The x-intercept must be defined within the interval of at least one dataset,
        otherwise there will be a degeneracy between the offset parameter and the coefficients
        of the polynomial
        """
        self.x_zero = None
        self.common_poly_ref = None

        self.time_interval = 1.000000000


    def initialize_model(self, mc, **kwargs):

        self.order = kwargs.get('order', 1)

        """ The user may decide to include the 0th order anyway - be aware of correlations with dataset offset!"""
        if kwargs.get('include_zero_point', False):
            self.starting_order = 0

        """ The user may decide to compute the polynomial parameters over a different time interval
            useful for long-term with very slow variations over a single day
        """
        self.time_interval = kwargs.get('time_interval', 1.000000000)

        """ If the polynomial is used as normalization factor, the first order must be included"""
        if self.normalization_model:
            self.starting_order = 0

        for i_order in range(self.starting_order, self.order+1):
            var = 'poly_c'+repr(i_order)
            self.list_pams_common.update([var])

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == 'polynomial_trend':
                self.common_poly_ref = common_ref
                break

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        try:
            mc.common_models[self.common_poly_ref].fix_list['x_zero'] = np.asarray([kwargs['x_zero'], 0.0000], dtype=np.double)
        except (KeyError, ValueError):
            if np.amin(dataset.x) < mc.Tref < np.amax(dataset.x):
                self.x_zero = mc.Tref
            elif not self.x_zero:
                self.x_zero = np.average(dataset.x)
            mc.common_models[self.common_poly_ref].fix_list['x_zero'] = np.asarray([self.x_zero, 0.0000])

    def compute(self, variable_value, dataset, x0_input=None):

        coeff = np.zeros(self.order+1)
        for i_order in range(self.starting_order, self.order+1):
            var = 'poly_c'+repr(i_order)
            coeff[i_order] = variable_value[var]

        """ In our array, coefficient are sorted from the lowest degree to the higher
        Numpy Polynomials requires the inverse order (from high to small) as input"""

        if x0_input is None:
            return polynomial.polyval((dataset.x-variable_value['x_zero'])/self.time_interval, coeff)
        else:
            return polynomial.polyval((x0_input+dataset.Tref-variable_value['x_zero'])/self.time_interval, coeff)


class LocalPolynomialTrend(AbstractModel):

    default_common = 'polynomial_trend'

    def __init__(self, *args, **kwargs):
        super(LocalPolynomialTrend, self).__init__(*args, **kwargs)

        self.model_class = 'local_polynomial_trend'

        self.list_pams_common = set()

        self.list_pams_dataset = {'x_zero'}

        self.order = 1
        self.starting_order = 1

        """
        The x-intercept must be defined within the interval of at least one dataset,
        otherwise there will be a degeneracy between the offset parameter and the coefficients
        of the polynomial
        """
        self.x_zero = {}

        self.time_interval = 1.000000000

    def initialize_model(self, mc, **kwargs):

        self.order = kwargs.get('order', 1)

        """ The user may decide to include the 0th order anyway - be aware of correlations with dataset offset!"""
        if kwargs.get('include_zero_point', False):
            self.starting_order = 0

        """ The user may decide to compute the polynomial parameters over a different time interval
            useful for long-term with very slow variations over a single day
        """
        self.time_interval = kwargs.get('time_interval', 1.000000000)

        """ If the polynomial is used as normalization factor, the first order must be included"""
        if self.normalization_model:
            self.starting_order = 0

        for i_order in range(self.starting_order, self.order+1):
            var = 'poly_c'+repr(i_order)
            self.list_pams_dataset.update([var])

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        try:
            self.fix_list[dataset.name_ref]['x_zero'] = np.asarray([kwargs['x_zero'], 0.0000], dtype=np.double)
        except (KeyError, ValueError):
            if np.amin(dataset.x) < dataset.Tref < np.amax(dataset.x):
                self.x_zero[dataset.name_ref] = dataset.Tref
            else:
                self.x_zero[dataset.name_ref] = np.average(dataset.x)
            self.fix_list[dataset.name_ref]['x_zero'] = np.asarray([self.x_zero[dataset.name_ref], 0.0000])

    def compute(self, variable_value, dataset, x0_input=None):

        coeff = np.zeros(self.order+1)
        for i_order in range(self.starting_order, self.order+1):
            var = 'poly_c'+repr(i_order)
            coeff[i_order] = variable_value[var]

        """ In our array, coefficient are sorted from the lowest degree to the highest """

        if x0_input is None:
            return polynomial.polyval((dataset.x-variable_value['x_zero'])/self.time_interval, coeff)
        else:
            return polynomial.polyval((x0_input+dataset.Tref-variable_value['x_zero'])/self.time_interval, coeff)

    def compute_alt(self, variable_value, dataset, x0_input=None):

        coeff = np.zeros(self.order+1)
        for i_order in range(self.starting_order, self.order+1):
            var = 'poly_c'+repr(i_order)
            coeff[-1-i_order] = variable_value[var]

        """ In our array, coefficient are sorted from the lowest degree to the highest """

        if x0_input is None:
            return self._polyval(coeff, (dataset.x-variable_value['x_zero'])/self.time_interval)
        else:
            return self._polyval(coeff, (x0_input+dataset.Tref-variable_value['x_zero'])/self.time_interval)

    def _polyval(p, x):
        y = np.zeros(x.shape, dtype=float)
        for i, v in enumerate(p):
            y *= x
            y += v
        return y

class SubsetPolynomialTrend(AbstractModel):

    default_common = 'polynomial_trend'

    def __init__(self, *args, **kwargs):
        super(SubsetPolynomialTrend, self).__init__(*args, **kwargs)

        self.model_class = 'subset_polynomial_trend'

        self.list_pams_common = set()

        self.list_pams_dataset = set()

        self.order = 1
        self.starting_order = 1

        """
        The x-intercept must be defined within the interval of at least one dataset,
        otherwise there will be a degeneracy between the offset parameter and the coefficients
        of the polynomial
        """
        self.x_zero = {}

        self.time_interval = 1.000000000

    def initialize_model(self, mc, **kwargs):

        self.order = kwargs.get('order', 1)

        """ The user may decide to include the 0th order anyway - be aware of correlations with dataset offset!"""
        if kwargs.get('include_zero_point', False):
            self.starting_order = 0

        """ The user may decide to compute the polynomial parameters over a different time interval
            useful for leng-term with very slow variations over a single day
        """
        self.time_interval = kwargs.get('time_interval', 1.000000000)

        """ If the polynomial is used as normalization factor, the first order must be included"""
        if self.normalization_model:
            self.starting_order = 0

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        if not dataset.submodel_flag:
            return

        for i_sub in range(0, dataset.submodel_flag):

            for i_order in range(self.starting_order, self.order+1):
                var_original = 'poly_c'+repr(i_order)
                var_subset = 'poly_sub'+repr(i_sub)+'_c'+repr(i_order)

                self._subset_transfer_priors(mc, dataset, var_original, var_subset)

            sub_dataset = dataset.x[(dataset.submodel_id==i_sub)]

            var_original = 'x_zero'
            var_subset = 'x_zero_sub'+repr(i_sub)
            try:
                xzero_ref = kwargs[var_subset] * 1.
                self.fix_list[dataset.name_ref][var_subset] = np.asarray([kwargs[var_subset], 0.0000], dtype=np.double)
            except (KeyError, ValueError):
                xzero_ref = np.average(sub_dataset)

            self.fix_list[dataset.name_ref][var_subset] = np.asarray([xzero_ref, 0.0000])

            self._subset_transfer_priors(mc, dataset, var_original, var_subset)


    def compute(self, variable_value, dataset, x0_input=None):

        if x0_input is None:
            y_output = np.zeros(dataset.n)
            xd_input = dataset.x
        else:
            y_output = x0_input * 0.
            xd_input = x0_input+dataset.Tref

        for i_sub in range(0,dataset.submodel_flag):

            x_zero_var = 'x_zero_sub'+repr(i_sub)
            x_input = xd_input-variable_value[x_zero_var]

            coeff = np.zeros(self.order+1)
            """ In our array, coefficient are sorted from the lowest degree to the highest """
            for i_order in range(self.starting_order, self.order+1):
                var = 'poly_sub'+repr(i_sub)+'_c'+repr(i_order)
                coeff[i_order] = variable_value[var]

            if x0_input is None:
                sel_data = (dataset.submodel_id==i_sub)
            else:
                original_dataset = dataset.x[(dataset.submodel_id==i_sub)] -variable_value[x_zero_var]
                sel_data = (x_input >= np.amin(original_dataset)) &  (x_input <= np.amax(original_dataset))
            y_output[sel_data] = polynomial.polyval(x_input[sel_data]/self.time_interval, coeff)

        return y_output
