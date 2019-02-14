from abstract_common import *
from abstract_model import *
import numpy.polynomial.polynomial


class CommonPolynomialTrend(AbstractCommon):
    ''' This class must be created for each planet in the system
            model_name is the way the planet is identified
    '''

    model_class = 'polynomial_trend'

    "polynomial trend up to 10th order"
    list_pams = {
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
        'poly_c0': [-1.0, 1.0],
        'poly_c1': [-1.0, 1.0],
        'poly_c2': [-1.0, 1.0],
        'poly_c3': [-1.0, 1.0],
        'poly_c4': [-1.0, 1.0],
        'poly_c5': [-1.0, 1.0],
        'poly_c6': [-1.0, 1.0],
        'poly_c7': [-1.0, 1.0],
        'poly_c8': [-1.0, 1.0],
        'poly_c9': [-1.0, 1.0]
    }

    default_spaces = {
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
        'poly_c0': ['Uniform', []], # 10 m/s/day would be already an unbelievable value
        'poly_c1': ['Uniform', []], # 10 m/s/day would be already an unbelievable value
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
    list_pams_common = {}
    list_pams_dataset = {}

    recenter_pams_dataset = {}

    order = 1
    starting_order = 1

    def initialize_model(self, mc, **kwargs):
        """ A special kind of initialization is required for this module, since it has to take a second dataset
        and check the correspondence with the points

        """
        if 'order' in kwargs:
            self.order = kwargs['order']

        """ The user may decide to include the 0th order anyway - be aware of correlations with dataset offset!"""
        try:
            if kwargs['include_zero_point']:
                self.starting_order = 0
        except:
            pass

        """ If the polynomial is used as normalization factor, the first order must be included"""
        if self.normalization_model:
            self.starting_order = 0

        for i_order in xrange(self.starting_order, self.order+1):
            var = 'poly_c'+repr(i_order)
            self.list_pams_common.update({var: None})

    def compute(self, variable_value, dataset, x0_input=None):

        coeff = np.zeros(self.order+1)
        for i_order in xrange(self.starting_order, self.order+1):
            var = 'poly_c'+repr(i_order)
            coeff[i_order] = variable_value[var]

        """ In our array, coefficient are sorted from the lowest degree to the highr
        Numpy Polinomials requires the inverse order (from high to small) as input"""

        if x0_input is None:
            return numpy.polynomial.polynomial.polyval(dataset.x0, coeff)
        else:
            return numpy.polynomial.polynomial.polyval(x0_input, coeff)


class LocalPolynomialTrend(AbstractModel):

    model_class = 'local_polynomial_trend'
    list_pams_common = {}
    list_pams_dataset = {}

    recenter_pams_dataset = {}

    order = 1
    starting_order = 1
    x_zero = {}

    def initialize_model(self, mc, **kwargs):
        """ A special kind of initialization is required for this module, since it has to take a second dataset
        and check the correspondence with the points

        """
        if 'order' in kwargs:
            self.order = kwargs['order']

        """ The user may decide to include the 0th order anyway - be aware of correlations with dataset offset!"""
        try:
            if kwargs['include_zero_point']:
                self.starting_order = 0
        except:
            pass

        """ If the polynomial is used as normalization factor, the first order must be included"""
        if self.normalization_model:
            self.starting_order = 0

        for i_order in xrange(self.starting_order, self.order+1):
            var = 'poly_c'+repr(i_order)
            self.list_pams_dataset.update({var: None})

    def setup_dataset(self, dataset, **kwargs):

        if 'x_zero' in kwargs:
            self.x_zero[dataset.name_ref] = kwargs['x_zero']
        else:
            self.x_zero[dataset.name_ref] = np.median(dataset.x0)

    def compute(self, variable_value, dataset, x0_input=None):

        coeff = np.zeros(self.order+1)
        for i_order in xrange(self.starting_order, self.order+1):
            var = 'poly_c'+repr(i_order)
            coeff[i_order] = variable_value[var]

        """ In our array, coefficient are sorted from the lowest degree to the highest
        Numpy Polinomials requires the inverse order (from high to small) as input"""

        if x0_input is None:
            return numpy.polynomial.polynomial.polyval(dataset.x0-self.x_zero[dataset.name_ref], coeff)
        else:
            return numpy.polynomial.polynomial.polyval(x0_input-self.x_zero[dataset.name_ref], coeff)
