from pyorbit.models.abstract_model import *
from numpy.polynomial import polynomial

class LocalCorrelation(AbstractModel):

    default_common = 'correlation'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'correlation'
        self.time_independent_model = True

        self.list_pams_common = set()
        self.list_pams_dataset = {'x_zero'}

        self.order = 1
        self.starting_order = 1
        self.threshold = 0.001

    def initialize_model(self, mc, **kwargs):

        # Threshold value to check if the datasets have the same dimensions
        # Well it looks like I have assumed that the user is careful enough
        self.threshold = kwargs.get('threshold', 0.001)
        self.order = kwargs.get('order', 1)
        self.correlated_val = kwargs.get('correlated_val', 'flux')

        """ The user may decide to include the 0th order anyway -
            be aware of correlations with dataset offset!"""
        if kwargs.get('include_zero_point', False):
            self.starting_order = 0


        for i_order in range(self.starting_order, self.order+1):
            var = 'corr_c'+repr(i_order)
            self.list_pams_dataset.update([var])


    def initialize_model_dataset(self, mc, dataset, **kwargs):

        try:
            self.fix_list[dataset.name_ref]['x_zero'] = np.asarray(
                [kwargs['x_zero'], 0.0000], dtype=np.double)
        except (KeyError, ValueError):
            self.fix_list[dataset.name_ref]['x_zero'] = np.asarray(
                [np.median(dataset.ancillary[self.correlated_val]), 0.0000], dtype=np.double)

    def compute(self, variable_value, dataset, x0_input=None):

        coeff = np.zeros(self.order+1)
        for i_order in range(self.starting_order, self.order+1):
            var = 'corr_c'+repr(i_order)
            coeff[i_order] = variable_value[var]

        """ In our array, coefficient are sorted from the lowest degree to the higher
        This is the order accepted by NumPy.polynomial.polynomial.polyval ,
        which is reversed with respect to from NumPy.polyval
        """
        if x0_input is None:
            return polynomial.polyval(
                dataset.ancillary[self.correlated_val] - variable_value['x_zero'],coeff)
        else:
            return 0.00
