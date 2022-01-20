from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
import numpy.polynomial.polynomial

class PolynomialTrend(AbstractModel):

    model_class = 'polynomial_trend'

    def __init__(self, *args, **kwargs):
        super(PolynomialTrend, self).__init__(*args, **kwargs)

        self.list_pams_common = {'x_zero'}

        """ in this model x_zero will be common to all the models, however it is treated
            as a specific
        """
        self.list_pams_dataset = set()


        self.recenter_pams_dataset = set()

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
            useful for leng-term with very slow variations over a single day
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
            return numpy.polynomial.polynomial.polyval((dataset.x-variable_value['x_zero'])/self.time_interval, coeff)
        else:
            return numpy.polynomial.polynomial.polyval((x0_input+dataset.Tref-variable_value['x_zero'])/self.time_interval, coeff)


class LocalPolynomialTrend(AbstractModel):

    model_class = 'local_polynomial_trend'

    def __init__(self, *args, **kwargs):
        super(LocalPolynomialTrend, self).__init__(*args, **kwargs)

        self.list_pams_common = set()

        self.list_pams_dataset = {'x_zero'}
        self.default_bounds = {'x_zero': [-10**6, 10**6]}
        self.default_spaces = {'x_zero': 'Linear'}
        self.default_priors = {'x_zero': ['Uniform', []]}

        self.recenter_pams_dataset = {}

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
            return numpy.polynomial.polynomial.polyval((dataset.x-variable_value['x_zero'])/self.time_interval, coeff)
        else:
            return numpy.polynomial.polynomial.polyval((x0_input+dataset.Tref-variable_value['x_zero'])/self.time_interval, coeff)

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
