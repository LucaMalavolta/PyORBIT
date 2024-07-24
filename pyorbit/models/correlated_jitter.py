import numpy.polynomial.polynomial
from pyorbit.models.abstract_model import *


class LocalCorrelatedJitter(AbstractModel):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'local_correlated_jitter'
        self.jitter_model = True

        self.list_pams_common = OrderedSet()
        self.list_pams_dataset = OrderedSet(['x_zero'])
        self.default_bounds = {'x_zero': [-10**5, 10**5]}
        self.default_spaces = {'x_zero': 'Linear'}
        self.default_priors = {'x_zero': ['Uniform', []]}

        self.recenter_pams_dataset = OrderedSet()

        self.order = 1
        self.x_vals = None
        self.x_mask = None
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

        """ HERE: we must associated the data from name_asc dataset to the one from name_ref
                remove that part from dataset.pyx
                Add a None option for the dataset
                Fix input_parser to accommodate the new changes
                Jitter must not be included in the analysis, but how about the offset?
                Or maybe I should just leave the zero point of the polynomial fit free?
        """

        self.x_vals = np.zeros(dataset_ref.n, dtype=np.double)
        self.x_mask = np.zeros(dataset_ref.n, dtype=bool)

        for i_date, v_date in enumerate(dataset_asc.x):
            match = np.where(np.abs(v_date-dataset_ref.x) < self.threshold)[0]
            self.x_vals[match] = dataset_asc.y[i_date]
            self.x_mask[match] = True

        print()
        print('Correlation model')
        print('Reference dataset:  ', kwargs['reference'])
        print('Associated dataset: ', kwargs['associated'])
        print('Cross-match between datasets: {0:d} out of {1:d} '.format(int(np.sum(self.x_mask)), dataset_ref.n))

        try:
            self.fix_list[dataset_asc.name_ref]['x_zero'] = np.asarray([kwargs['x_zero'], 0.0000], dtype=np.double)
        except (KeyError, ValueError):
            self.fix_list[dataset_asc.name_ref]['x_zero'] = np.asarray([np.amin(self.x_vals[self.x_mask]), 0.0000], dtype=np.double)

        for i_order in range(1, self.order+1):
            par = 'c'+repr(i_order)
            self.list_pams_dataset.update(par)
            self.default_bounds.update({par: [0.0, 10**6]})
            self.default_spaces.update({par: 'Linear'})
            self.default_priors.update({par: ['Uniform', []]})

    def compute(self, parameter_values, dataset, x0_input=None):

        coeff = np.zeros(self.order+1)
        for i_order in range(1, self.order+1):
            par = 'c'+repr(i_order)
            coeff[i_order] = parameter_values[par]
        x_zero = parameter_values['x_zero']

        """ In our array, coefficient are sorted from the lowest degree to the higher
        This is the order accepted by NumPy.polynomial.polynomial.polyval ,
        which is reversed with respect to from NumPy.polyval
        """
        if x0_input is None:
            return np.where(self.x_mask, numpy.polynomial.polynomial.polyval(self.x_vals-x_zero, coeff), 0.0)
        else:
            return numpy.polynomial.polynomial.polyval(x0_input-x_zero, coeff)

