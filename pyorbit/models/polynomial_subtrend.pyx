from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
from numpy.polynomial import polynomial


class SubPolynomialTrend(AbstractModel):

    model_class = 'sub_polynomial_trend'

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

        for i_sub in range(0,dataset.submodel_flag):

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

        if x0_input is None:
            y_output = np.zeros(dataset.n)
            x_input = dataset.x-variable_value['x_zero']
        else:
            y_output = x0_input * 0.
            x_input = x0_input+dataset.Tref-variable_value['x_zero']

        for i_sub in range(0,dataset.submodel_flag):

            coeff = np.zeros(self.order+1)
            """ In our array, coefficient are sorted from the lowest degree to the highest """
            for i_order in range(self.starting_order, self.order+1):
                var = 'poly_sub'+repr(i_sub)+'_c'+repr(i_order)
                coeff[i_order] = variable_value[var]

            if x0_input is None:
                sel_data = (dataset.submodel_id==i_sub)
            else:
                original_dataset = dataset.x[(dataset.submodel_id==i_sub)] -variable_value['x_zero']
                sel_data = (x_input >= np.amin(original_dataset)) &  (x_input <= np.amax(original_dataset))
            y_output[sel_data] = polynomial.polyval(x_input[sel_data]/self.time_interval, coeff)

        return y_output
