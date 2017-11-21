from abstract_model import *

class Correlation_SingleDataset(AbstractModel):

    model_class = 'correlations'
    list_pams_common = {}
    list_pams_dataset = {
        'c1': 'U',
        'c2': 'U',
        'c3': 'U',
        'c4': 'U',
        'c5': 'U',
        'c6': 'U',
        'c7': 'U',
        'c8': 'U',
        'c9': 'U'
    }

    default_bounds = {
        'c1': [-10**9,10**9],
        'c2': [-10**9,10**9],
        'c3': [-10**9,10**9],
        'c4': [-10**9,10**9],
        'c5': [-10**9,10**9],
        'c6': [-10**9,10**9],
        'c7': [-10**9,10**9],
        'c8': [-10**9,10**9],
        'c9': [-10**9,10**9]
    }

    order = 1
    x_vals = None
    x_mask = None
    x_zero = 0.000
    threshold = 0.001

    def initialize_model(self, dataset_ref, dataset_asc, **kwargs):
        """ A special kind of initialization is required for this module, since it has to take a second dataset
        and check the corrispondence with the points

        """
        if 'threshold' in kwargs:
            self.threshold = kwargs['threshold']
        if 'order' in kwargs:
            self.order = kwargs['order']

        self.fix_list[dataset_ref.name_ref] = {}
        for ii in xrange(self.order+1, 10):
            self.fix_list[dataset_ref.name_ref]['c'+repr(ii)] = 0.0000

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
        for ii in xrange(1,10):
            var = 'c'+repr(ii)
            coeff[ii] = variable_value[var]

        """ In our array, coefficient are sorted from the lowest degree to the highr
        Numpy Polinomials requires the inverse order (from high to small) as input"""
        return np.where(self.x_mask, np.polynomial.polynomial.polyval(self.x_vals-self.x_zero, coeff), 0.0)
