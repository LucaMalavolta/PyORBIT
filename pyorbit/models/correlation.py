from pyorbit.models.abstract_model import *
from numpy.polynomial import polynomial
from scipy.interpolate import interp1d

class LocalCorrelation(AbstractModel):

    default_common = 'correlation'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'correlation'
        self.unitary_model = False
        self.normalization_model = False

        self.list_pams_common = OrderedSet()
        self.list_pams_dataset = OrderedSet(['x_zero'])

        self.use_median_xzero = True

        self.order = 1
        self.starting_order = 1
        self.x_zero = 0.00
        self.baseline_value = 0.000

        self.interpolated = {}

    def initialize_model(self, mc, **kwargs):

        # Threshold value to check if the datasets have the same dimensions
        # Well it looks like I have assumed that the user is careful enough
        self.order = kwargs.get('order', self.order)
        self.x_zero = kwargs.get('x_zero', self.x_zero)
        self.use_median_xzero = kwargs.get('use_median_xzero', self.use_median_xzero)

        self.correlated_val = kwargs.get('correlated_variable', 'corr')

        """ If the polynomial is used as normalization factor, the first order must be included"""
        if self.normalization_model:
            self.starting_order = 0
            self.baseline_value = kwargs.get('baseline_value', 1.0)
        else:
            self.starting_order = 1
            self.baseline_value = kwargs.get('baseline_value', 0.0)


        """ The user may decide to include the 0th order anyway -
            be aware of correlations with dataset offset!"""
        if kwargs.get('include_zero_point', self.include_zero_point):
            self.starting_order = 0

        """ zero point should be excluded  """
        if kwargs.get('exclude_zero_point', self.exclude_zero_point):
            self.starting_order = 1


        for i_order in range(self.starting_order, self.order+1):
            par = 'corr_c'+repr(i_order)
            self.list_pams_dataset.update([par])

    def initialize_model_dataset(self, mc, dataset, **kwargs):


        if self.use_median_xzero:
            x_zero = np.median(dataset.ancillary[self.correlated_val])
        else:
            x_zero = kwargs.get('x_zero', self.x_zero)

        self.fix_list[dataset.name_ref]['x_zero'] = np.asarray([x_zero, 0.0000], dtype=np.double)

        self.interpolated[dataset.name_ref]=interp1d(
                dataset.ancillary['time'],
                dataset.ancillary[self.correlated_val],
                bounds_error=False,
                fill_value=(np.amin(dataset.ancillary[self.correlated_val]), np.amax(dataset.ancillary[self.correlated_val])))

    def compute(self, parameter_values, dataset, x0_input=None):


        coeff = np.zeros(self.order+1)
        """ starting from the first order coefficient, as the constant may be added as baseline"""
        for i_order in range(1, self.order+1):
            par = 'corr_c'+repr(i_order)
            coeff[i_order] = parameter_values[par]

        """ In our array, coefficient are sorted from the lowest degree to the higher
        This is the order accepted by NumPy.polynomial.polynomial.polyval ,
        which is reversed with respect to from NumPy.polyval
        """
        if x0_input is None:
            return polynomial.polyval(
                dataset.ancillary[self.correlated_val] - parameter_values['x_zero'], coeff) + parameter_values.get('corr_c0', self.baseline_value)
        else:
            t = x0_input + dataset.Tref

            return polynomial.polyval(
                self.interpolated[dataset.name_ref](t) - parameter_values['x_zero'], coeff)+ parameter_values.get('corr_c0', self.baseline_value)
