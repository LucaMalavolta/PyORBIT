from pyorbit.models.abstract_model import *
from numpy.polynomial import polynomial
from numpy.polynomial import Polynomial
from scipy.optimize import minimize

class ComplexCorrelation(AbstractModel):

    default_common = 'complex_correlation'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'complex_correlation'
        self.time_independent_model = True
        self.residuals_analysis = True

        self.list_pams_common = OrderedSet(['x_zero'])
        self.list_pams_dataset = OrderedSet()

        self.order = 1
        self.starting_order = 1
        #self.threshold = 0.001
        self.common_correlation = None

    def initialize_model(self, mc, **kwargs):

        # Threshold value to check if the datasets have the same dimensions
        # Well it looks like I have assumed that the user is careful enough
        #self.threshold = kwargs.get('threshold', 0.001)
        self.order = kwargs.get('order', 1)

        """ The user may decide to include the 0th order anyway -
            be aware of correlations with dataset offset!"""
        if kwargs.get('include_zero_point', False):
            self.starting_order = 0


        for i_order in range(self.starting_order, self.order+1):
            par = 'corr_c'+repr(i_order)
            self.list_pams_common.update([par])

        try:
            self.x_dataset = kwargs.get('dataset_x')
            self.y_dataset = kwargs.get('dataset_y')
        except:
            print('x_dataset and y_dataset must be specified')

        self.gp_before_correlation = kwargs.get('gp_before_correlation', False)

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == 'complex_correlation':
                self.common_correlation = common_ref
                break

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        try:
             mc.common_models[self.common_correlation].fix_list['x_zero'] = \
                np.asarray([kwargs['x_zero'], 0.0000], dtype=np.double)
        except (KeyError, ValueError):
            if dataset.name_ref == self.x_dataset:
                mc.common_models[self.common_correlation].fix_list['x_zero'] = np.asarray([np.median(dataset.y), 0.0000])


    def _compute_distance(self, x_comp):
            y_comp =self.temp_poly(x_comp)
            return (self.temp_xx-x_comp)**2 + (self.temp_yy-y_comp)**2

    def compute(self, parameter_values, x_dataset, y_dataset, x0_input=None):

        coeff = np.zeros(self.order+1)
        for i_order in range(self.starting_order, self.order+1):
            par = 'corr_c'+repr(i_order)
            coeff[i_order] = parameter_values[par]

        """ In our array, coefficient are sorted from the lowest degree to the higher
        This is the order accepted by NumPy.polynomial.polynomial.polyval ,
        which is reversed with respect to from NumPy.polyval
        """
        if x0_input is None:

            self.temp_poly = Polynomial(coeff)
            out_xx = np.zeros_like(x_dataset.residuals)

            for ii in range(0, x_dataset.n):
                self.temp_xx = x_dataset.residuals[ii] - parameter_values['x_zero']
                self.temp_yy = y_dataset.residuals[ii]

                out_xx[ii] = minimize(self._compute_distance, self.temp_xx)['x']
            out_yy = self.temp_poly(out_xx- parameter_values['x_zero'])

            return out_xx, out_yy

        else:
            return 0.00


