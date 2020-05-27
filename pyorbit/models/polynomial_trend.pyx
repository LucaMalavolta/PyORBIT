from pyorbit.classes.common import *
from pyorbit.models.abstract_common import *
from pyorbit.models.abstract_model import *
import numpy.polynomial.polynomial


class CommonPolynomialTrend(AbstractCommon):

    model_class = 'polynomial_trend'

    "polynomial trend up to 10th order"
    list_pams = {
        'x_zero',
        'poly_c0',  # order 0
        'poly_c1',  # order 1
        'poly_c2',  # order 2
        'poly_c3',  # order 3
        'poly_c4',  # order 4
        'poly_c5',  # order 5
        'poly_c6',  # order 6
        'poly_c7',  # order 7
        'poly_c8',  # order 8
        'poly_c9',  # order 9
    }

    """These default boundaries are used when the user does not define them in the yaml file"""
    default_bounds = {
        'x_zero': [-10 ** 6, 10 ** 6],
        'poly_c0': [-10 ** 5, 10 ** 6],
        'poly_c1': [-10 ** 5, 10 ** 6],
        'poly_c2': [-10 ** 6, 10 ** 6],
        'poly_c3': [-10 ** 6, 10 ** 6],
        'poly_c4': [-10 ** 6, 10 ** 6],
        'poly_c5': [-10 ** 6, 10 ** 6],
        'poly_c6': [-10 ** 6, 10 ** 6],
        'poly_c7': [-10 ** 6, 10 ** 6],
        'poly_c8': [-10 ** 6, 10 ** 6],
        'poly_c9': [-10 ** 6, 10 ** 6]
    }

    default_spaces = {
        'x_zero': 'Linear',
        'poly_c0': 'Linear',  # order 1
        'poly_c1': 'Linear',  # order 1
        'poly_c2': 'Linear',  # order 2
        'poly_c3': 'Linear',  # order 3
        'poly_c4': 'Linear',  # order 4
        'poly_c5': 'Linear',  # order 5
        'poly_c6': 'Linear',  # order 6
        'poly_c7': 'Linear',  # order 7
        'poly_c8': 'Linear',  # order 8
        'poly_c9': 'Linear',  # order 9
    }

    default_priors = {
        'x_zero': ['Uniform', []],
        'poly_c0': ['Uniform', []],
        'poly_c1': ['Uniform', []],
        'poly_c2': ['Uniform', []],
        'poly_c3': ['Uniform', []],
        'poly_c4': ['Uniform', []],
        'poly_c5': ['Uniform', []],
        'poly_c6': ['Uniform', []],
        'poly_c7': ['Uniform', []],
        'poly_c8': ['Uniform', []],
        'poly_c9': ['Uniform', []]
    }

    default_fixed = {}

    recenter_pams = {}


class PolynomialTrend(AbstractModel):

    model_class = 'polynomial_trend'

    def __init__(self, *args, **kwargs):
        super(PolynomialTrend, self).__init__(*args, **kwargs)

        self.list_pams_common = {'x_zero': None}

        """ in this model x_zero will be common to all the models, however it is treated
            as a specific
        """
        self.list_pams_dataset = {}


        self.recenter_pams_dataset = {}

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

        if 'order' in kwargs:
            self.order = kwargs['order']

        """ The user may decide to include the 0th order anyway - be aware of correlations with dataset offset!"""
        try:
            if kwargs['include_zero_point']:
                self.starting_order = 0
        except:
            self.starting_order = 1

        """ The user may decide to compute the polynomial parameters over a different time interval
            useful for leng-term with very slow variations over a single day
        """
        try:
            if kwargs['time_interval']:
                self.time_interval = kwargs['time_interval']
            print('TIME INTERVAL: ', self.time_interval)
        except:
            self.time_interval = 1.

        """ If the polynomial is used as normalization factor, the first order must be included"""
        if self.normalization_model:
            self.starting_order = 0

        for i_order in range(self.starting_order, self.order+1):
            var = 'poly_c'+repr(i_order)
            self.list_pams_common.update({var: None})

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == 'polynomial_trend':
                self.common_poly_ref = common_ref
                break

    def setup_dataset(self, mc, dataset, **kwargs):

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

        self.list_pams_common = {}

        self.list_pams_dataset = {'x_zero': None}
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

        if 'order' in kwargs:
            self.order = kwargs['order']

        """ The user may decide to include the 0th order anyway - be aware of correlations with dataset offset!"""
        try:
            if kwargs['include_zero_point']:
                self.starting_order = 0
            self.starting_order
        except:
            self.starting_order = 1

        """ The user may decide to compute the polynomial parameters over a different time interval
            useful for leng-term with very slow variations over a single day
        """
        try:
            if kwargs['time_interval']:
                self.time_interval = kwargs['time_interval']
            self.time_interval
        except:
            self.time_interval = 1.

        """ If the polynomial is used as normalization factor, the first order must be included"""
        if self.normalization_model:
            self.starting_order = 0

        for i_order in range(self.starting_order, self.order+1):
            var = 'poly_c'+repr(i_order)
            self.list_pams_dataset.update({var: None})

    def setup_dataset(self, mc, dataset, **kwargs):

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

        """ In our array, coefficient are sorted from the lowest degree to the highest
        Numpy Polynomials requires the inverse order (from high to small) as input"""

        if x0_input is None:
            return numpy.polynomial.polynomial.polyval((dataset.x-variable_value['x_zero'])/self.time_interval, coeff)
        else:
            return numpy.polynomial.polynomial.polyval((x0_input+dataset.Tref-variable_value['x_zero'])/self.time_interval, coeff)
