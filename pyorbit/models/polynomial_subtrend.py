from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
import numpy.polynomial.polynomial


class LocalPolynomialTrend(AbstractModel):

    model_class = 'local_polynomial_trend'

    def __init__(self, *args, **kwargs):
        super(LocalPolynomialTrend, self).__init__(*args, **kwargs)

        self.list_pams_common = set()

        self.list_pams_dataset = {'x_zero'}
        self.default_bounds = {'x_zero': [-10**6, 10**6]}
        self.default_spaces = {'x_zero': 'Linear'}
        self.default_priors = {'x_zero': ['Uniform', []]}

        self.order = 1
        self.starting_order = 1

        """
        The x-intercept must be defined within the interval of at least one dataset,
        otherwise there will be a degeneracy between the offset parameter and the coefficients
        of the polynomial
        """
        self.x_zero = {}
        self.selection_matrix = {}

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

        self.selection_matrix[]

        for i_sub in range(0,dataset.submodel_flag):

            self.fix_list[dataset.name_ref]['sub'+repr(i_order)] = {}

            for i_order in range(self.starting_order, self.order+1):
                var = 'poly_sub'+repr(i_sub)+'_c'+repr(i_order)
                self.list_pams_dataset.update([var])

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


        for i_num in range(0,10):
            output[(sel_array==i_num)] 

        if x0_input is None:
            return numpy.polynomial.polynomial.polyval((dataset.x-variable_value['x_zero'])/self.time_interval, coeff)
        else:
            return numpy.polynomial.polynomial.polyval((x0_input+dataset.Tref-variable_value['x_zero'])/self.time_interval, coeff)

    def _polyval(p, x):
        y = np.zeros(x.shape, dtype=float)
        for i, v in enumerate(p):
            y *= x
            y += v
        return y
