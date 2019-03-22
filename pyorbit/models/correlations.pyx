import numpy.polynomial.polynomial
from pyorbit.models.abstract_common import *
from pyorbit.models.abstract_model import *

class Correlation_SingleDataset(AbstractModel):

    model_class = 'correlations'
    time_independent_model = True

    def __init__(self, *args, **kwargs):
        super(Correlation_SingleDataset, self).__init__(*args, **kwargs)

        self.list_pams_common = {}
        self.list_pams_dataset = {}
        self.default_bounds = {}
        self.default_spaces = {}
        self.default_priors = {}

        self.recenter_pams_dataset = {}

        self.order = 1
        self.x_vals = None
        self.x_mask = None
        self.x_zero = 0.000
        self.threshold = 0.001

    def initialize_model(self, mc, **kwargs):
        """ A special kind of initialization is required for this module, since it has to take a second dataset
        and check the correspondence with the points

        """
        dataset_ref = mc.dataset_dict[kwargs['reference']]
        dataset_asc = mc.dataset_dict[kwargs['associated']]

        if 'threshold' in kwargs:
            self.threshold = kwargs['threshold']
        if 'order' in kwargs:
            self.order = kwargs['order']

        for i_order in xrange(1, self.order+1):
            var = 'c'+repr(i_order)
            self.list_pams_dataset.update({var: None})
            self.default_bounds.update({var: [-10**9, 10**9]})
            self.default_spaces.update({var: 'Linear'})
            self.default_priors.update({var: ['Uniform', []]})

        self.x_vals = np.zeros(dataset_ref.n, dtype=np.double)
        self.x_mask = np.zeros(dataset_ref.n, dtype=bool)

        """ HERE: we must associated the data from name_asc dataset to the one from name_ref
                remove that part from dataset.pyx
                Add a None option for the dataset
                Fix input_parser to accomodate the new changes
                Jitter must not be included in the analysis, but how about the offset? 
                Or maybe I should just leave the zero point of the polynomial fit free?
        """

        for i_date, v_date in enumerate(dataset_asc.x):
            match = np.where(np.abs(v_date-dataset_ref.x) < self.threshold)[0]
            self.x_vals[match] = dataset_asc.y[i_date]
            self.x_mask[match] = True

        if 'x_zero' in kwargs:
            self.x_zero = kwargs['x_zero']
        else:
            self.x_zero = np.median(self.x_vals[self.x_mask])

    def compute(self, variable_value, dataset, x0_input=None):

        coeff = np.zeros(self.order+1)
        for i_order in xrange(1, self.order+1):
            var = 'c'+repr(i_order)
            coeff[i_order] = variable_value[var]
        """ In our array, coefficient are sorted from the lowest degree to the higher
        Numpy Polinomials requires the inverse order (from high to small) as input"""
        if x0_input is None:
            return np.where(self.x_mask, numpy.polynomial.polynomial.polyval(self.x_vals-self.x_zero, coeff), 0.0)
        else:
            return numpy.polynomial.polynomial.polyval(x0_input-self.x_zero, coeff)
