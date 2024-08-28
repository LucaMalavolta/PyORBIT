from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *


class Harmonics(AbstractModel):

    def __init__(self, *args, **kwargs):
        super(Harmonics, self).__init__(*args, **kwargs)

        self.model_class = 'harmonics'

        self.list_pams_common = OrderedSet()
        self.list_pams_dataset = OrderedSet()

        self.sine_harmonics = [1, 2]
        self.cosine_harmonics = [1]
        self.use_t0 = False
        self.use_common_period = True

    def initialize_model(self, mc, **kwargs):


        if kwargs.get('sine_harmonics_selection', False):
            self.sine_harmonics = kwargs['sine_harmonics_selection']
        else:
            n_sines = kwargs.get('sine_harmonics', 2)
            self.sine_harmonics = np.arange(1, n_sines+1, dtype=np.int16)
        if self.sine_harmonics is None or self.sine_harmonics == 'None':
            self.sine_harmonics = []

        for i_order in self.sine_harmonics:
            if i_order == 0: continue
            par = 'amp_S'+repr(i_order)
            self.list_pams_dataset.update([par])

        if kwargs.get('cosine_harmonics_selection', False):
            self.cosine_harmonics = kwargs['cosine_harmonics_selection']
        else:
            n_cosines = kwargs.get('cosine_harmonics', 1)
            self.cosine_harmonics = np.arange(1, n_cosines+1, dtype=np.int16)
        if self.cosine_harmonics is None or self.cosine_harmonics == 'None':
            self.cosine_harmonics = []

        for i_order in self.cosine_harmonics:
            if i_order == 0: continue
            par = 'amp_C'+repr(i_order)
            self.list_pams_dataset.update([par])



        if kwargs.get('use_common_independent_phases', False):
            if kwargs.get('use_T0', False):
                print(
                    "Harmonics model: independent phases and T0s are not compatible options")
            for i_order in self.sine_harmonics:
                par = 'pha_S'+repr(i_order)
                self.list_pams_common.update([par])
            for i_order in self.cosine_harmonics:
                par = 'pha_C'+repr(i_order)
                self.list_pams_common.update([par])
        elif kwargs.get('use_independent_phases', False):
            if kwargs.get('use_T0', False) or kwargs.get('use_common_T0', False):
                print(
                    "Harmonics model: independent phases and T0s are not compatible options")
            for i_order in self.sine_harmonics:
                par = 'pha_S'+repr(i_order)
                self.list_pams_dataset.update([par])
            for i_order in self.cosine_harmonics:
                par = 'pha_C'+repr(i_order)
                self.list_pams_dataset.update([par])
        elif kwargs.get('use_common_T0', False):
            self.use_t0 = True
            self.list_pams_common.update(['T0'])
        elif kwargs.get('use_T0', False):
            self.use_t0 = True
            self.list_pams_dataset.update(['T0'])
        elif kwargs.get('use_common_phase', False):
            self.list_pams_common.update(['phase'])
        else:
            self.list_pams_dataset.update(['phase'])

        if kwargs.get('use_common_period', self.use_common_period):
            self.list_pams_common.update(['P'])
        else:
            self.list_pams_dataset.update(['P'])

    def compute(self, parameter_values, dataset, x0_input=None):

        if self.use_t0:
            if x0_input is None:
                phi_value = (
                    dataset.x - parameter_values['T0']) / parameter_values['P'] * 2 * np.pi
            else:
                phi_value = (x0_input + dataset.Tref -
                             parameter_values['T0']) / parameter_values['P'] * 2 * np.pi
        else:
            if x0_input is None:
                phi_value = dataset.x0 / parameter_values['P'] * 2 * np.pi
            else:
                phi_value = x0_input / parameter_values['P'] * 2 * np.pi

        harmonics_output = 0.000 * phi_value

        phase = parameter_values.get('phase', 0.00)

        for i_order in self.sine_harmonics:
            par = 'amp_S'+repr(i_order)
            pha = parameter_values.get('pha_S'+repr(i_order),  phase)
            harmonics_output += parameter_values[par] * \
                np.sin(phi_value * i_order + pha/180.*np.pi)

        for i_order in self.cosine_harmonics:
            par = 'amp_C'+repr(i_order)
            pha = parameter_values.get('pha_C'+repr(i_order),  phase)
            harmonics_output += parameter_values[par] * \
                np.cos(phi_value * i_order + pha/180.*np.pi)

        return harmonics_output
