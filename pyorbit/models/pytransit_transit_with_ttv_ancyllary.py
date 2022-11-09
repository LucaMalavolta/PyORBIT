
from pyorbit.subroutines.common import np, convert_rho_to_a, convert_b_to_i
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import AbstractTransit

try:
    from pytransit import QuadraticModel
    from pytransit import RoadRunnerModel
    from pytransit import QPower2Model
except ImportError:
    pass


class PyTransit_Transit_With_TTV_Ancillary(AbstractModel, AbstractTransit):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            from pytransit import QuadraticModel
        except ImportError:
            print("ERROR: PyTransit not installed, this will not work")
            quit()

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = {
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'omega',  # argument of pericenter (in radians)
            'R_Rs',  # planet radius (in units of stellar radii)
        }

        self.pytransit_models = {}
        self.pytransit_plot = {}

        """ Dataset-specific time of transit boundaries are stored here"""
        self.transit_time_boundaries = {}


    def initialize_model(self, mc, **kwargs):
        """ Force the use of the central time of transit"""
        self.use_time_of_transit = True

        self._prepare_planetary_parameters(mc, **kwargs)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

        # Qui: leggere file con i tempi di transito e gli id

        self.list_pams_common.discard('Tc')
        self.list_pams_dataset.update(['Tc'])

    def initialize_model_dataset(self, mc, dataset, **kwargs):
        self._prepare_dataset_options(mc, dataset, **kwargs)

        if self.limb_darkening_model == 'quadratic':
            self.pytransit_models[dataset.name_ref] = QuadraticModel()
            self.pytransit_plot[dataset.name_ref] = QuadraticModel()

        self.pytransit_models[dataset.name_ref].set_data(dataset.x0,
                                                         exptimes=self.code_options[dataset.name_ref]['exp_time'],
                                                         nsamples=self.code_options[dataset.name_ref]['sample_factor'])

    def define_special_variable_properties(self,
                                           ndim,
                                           output_lists,
                                           dataset_name,
                                           var):

        if var == 'Tc' and (var not in self.bounds[dataset_name]):
            self.bounds[dataset_name][var] = self.code_options[dataset_name]['Tc_boundaries']
        return ndim, output_lists, False

    def compute(self, variable_value, dataset, x0_input=None):
        """
        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """

        pams_a, pams_i = self.retrieve_ai(variable_value)
        pams_t0 = self.retrieve_t0(variable_value, dataset.Tref)
        omega_rad = variable_value['omega'] / 180. * np.pi

        for var, i_var in self.ldvars.items():
            self.ld_vars[i_var] = variable_value[var]

        if x0_input is None:
            return self.pytransit_models[dataset.name_ref].evaluate_ps(
                variable_value['R_Rs'],
                self.ld_vars,
                pams_t0, variable_value['P'], pams_a, pams_i, variable_value['e'], omega_rad) - 1.

        else:
            self.pytransit_plot[dataset.name_ref].set_data(x0_input,
                                                           exptimes=self.code_options[dataset.name_ref]['exp_time'],
                                                           nsamples=self.code_options[dataset.name_ref]['sample_factor'])

            return self.pytransit_plot[dataset.name_ref].evaluate_ps(
                variable_value['R_Rs'],
                self.ld_vars,
                pams_t0, variable_value['P'], pams_a, pams_i, variable_value['e'], omega_rad) - 1.
