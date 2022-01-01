from pyorbit.classes.common import *
from pyorbit.models.abstract_model import *


class Harmonics(AbstractModel):

    model_class = 'harmonics'

    def __init__(self, *args, **kwargs):
        super(Harmonics, self).__init__(*args, **kwargs)

        self.list_pams_common = {}
        self.list_pams_dataset = {}

        self.recenter_pams_dataset = {}

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

