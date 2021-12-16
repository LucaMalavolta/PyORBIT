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


def _polyval(p, x):
    y = np.zeros(x.shape, dtype=float)
    for i, v in enumerate(p):
        y *= x
        y += v
    return y


class Harmonics(AbstractModel):

    model_class = 'harmonics'

    def __init__(self, *args, **kwargs):
        super(Harmonics, self).__init__(*args, **kwargs)

        self.list_pams_common = {}
        self.list_pams_dataset = {}

        self.recenter_pams_dataset = {}

        #self.order = 1
        #self.starting_order = 1

        """
        The x-intercept must be defined within the interval of at least one dataset,
        otherwise there will be a degeneracy between the offset parameter and the coefficients
        of the polynomial
        """
        #self.x_zero = None
        #self.common_poly_ref = None

        #self.time_interval = 1.000000000


        self.sine_harmonics = [1, 2]
        self.cosine_harmonics = [1]
        self.use_t0 = False
        self.use_common_period = True

    def initialize_model(self, mc, **kwargs):

        if 'sine_harmonics' in kwargs:
            self.sine_harmonics = np.arange(1, kwargs['sine_harmonics']+1, dtpye=np.int16)
        if 'sine_harmonics_selection' in kwargs:
            self.sine_harmonics = kwargs['sine_harmonics_selection']

        for i_order in self.sine_harmonics:
            var = 'amp_S'+repr(i_order)
            self.list_pams_dataset.update({var: None})

        if 'cosine_harmonics' in kwargs:
            self.cosine_harmonics = np.arange(1, kwargs['cosine_harmonics']+1, dtpye=np.int16)
        if 'cosine_harmonics_selection' in kwargs:
            self.cosine_harmonics = kwargs['cosine_harmonics_selection']

        for i_order in self.cosine_harmonics:
            var = 'amp_C'+repr(i_order)
            self.list_pams_dataset.update({var: None})

        if kwargs.get('use_common_independent_phases', False):
            if kwargs.get('use_T0', False):
                print("Harmonics model: independent phases and T0 are not compatible options") 
            for i_order in self.sine_harmonics:
                var = 'pha_S'+repr(i_order)
                self.list_pams_common.update({var: None})
            for i_order in self.cosine_harmonics:
                var = 'pha_C'+repr(i_order)
                self.list_pams_common.update({var: None})
        elif kwargs.get('use_independent_phases', False):
            if kwargs.get('use_T0', False) or kwargs.get('use_common_T0', False):
                print("Harmonics model: independent phases and T0 are not compatible options") 
            for i_order in self.sine_harmonics:
                var = 'pha_S'+repr(i_order)
                self.list_pams_dataset.update({var: None})
            for i_order in self.cosine_harmonics:
                var = 'pha_C'+repr(i_order)
                self.list_pams_dataset.update({var: None})
        elif kwargs.get('use_common_T0', False):
            self.use_t0 = True
            self.list_pams_common.update({'T0': None})
        elif kwargs.get('use_T0', False):
            self.use_t0 = True
            self.list_pams_dataset.update({'T0': None})
        elif kwargs.get('use_common_phase', False):
            self.list_pams_common.update({'phase': None})
        else:
            self.list_pams_dataset.update({'phase': None})

        if kwargs.get('use_common_period', self.use_common_period):
            self.list_pams_common.update({'P': None})
        else:
            self.list_pams_dataset.update({'P': None})




    def compute(self, variable_value, dataset, x0_input=None):

        if self.use_t0:
            if x0_input is None:
                phi_value = (dataset.x - variable_value['T0']) / variable_value['P'] * 2 * np.pi
            else:
                phi_value = (x0_input + dataset.Tref - variable_value['T0']) / variable_value['P'] * 2 * np.pi
        else:
            if x0_input is None:
                phi_value = dataset.x0 / variable_value['P'] * 2 * np.pi
            else:
                phi_value = x0_input / variable_value['P'] * 2 * np.pi

        harmonics_output = 0.000 * phi_value

        phase = variable_value.get('phase', 0.00)

        for i_order in self.sine_harmonics:
            var = 'amp_S'+repr(i_order)
            pha = variable_value.get('pha_S'+repr(i_order),  phase)
            harmonics_output += variable_value[var] * np.sin(phi_value * i_order + pha)

        for i_order in self.cosine_harmonics:
            var = 'amp_C'+repr(i_order)
            pha = variable_value.get('pha_C'+repr(i_order),  phase)
            harmonics_output += variable_value[var] * np.cos(phi_value * i_order + pha)

        return harmonics_output

