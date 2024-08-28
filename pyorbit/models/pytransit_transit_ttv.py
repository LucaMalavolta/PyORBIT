
from pyorbit.subroutines.common import np, constants, OrderedSet
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import AbstractTransit

try:
    from pytransit import QuadraticModel
    from pytransit import RoadRunnerModel
    from pytransit import QPower2Model
except (ModuleNotFoundError,ImportError):
    pass


class PyTransit_Transit_TTV(AbstractModel, AbstractTransit):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            from pytransit import QuadraticModel
        except (ModuleNotFoundError,ImportError):
            print("ERROR: PyTransit not installed, this will not work")
            quit()

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = OrderedSet([
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'omega',  # argument of pericenter (in radians)
            'R_Rs',  # planet radius (in units of stellar radii)
        ])

        self.pytransit_models = {}
        self.pytransit_plot = {}

    def initialize_model(self, mc, **kwargs):
        """ Force the use of the time of inferior conjunction"""
        mc.common_models[self.planet_ref].use_time_inferior_conjunction = True

        self.use_roadrunner = kwargs.get('use_roadrunner', True)
        if self.use_roadrunner:
            print('Using RoadRunner Model from PyTransit')

        self._prepare_planetary_parameters(mc, **kwargs)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

        self.list_pams_common.discard('Tc')
        self.list_pams_dataset.update(['Tc'])

    def initialize_model_dataset(self, mc, dataset, **kwargs):
        self._prepare_dataset_options(mc, dataset, **kwargs)

        if self.use_roadrunner:
            self.pytransit_models[dataset.name_ref] = RoadRunnerModel(self.limb_darkening_model)
            self.pytransit_plot[dataset.name_ref] = RoadRunnerModel(self.limb_darkening_model)
        elif self.limb_darkening_model == 'quadratic':
            self.pytransit_models[dataset.name_ref] = QuadraticModel()
            self.pytransit_plot[dataset.name_ref] = QuadraticModel()

        self.pytransit_models[dataset.name_ref].set_data(dataset.x0,
                                                            exptimes=self.code_options[dataset.name_ref]['exp_time'],
                                                            nsamples=self.code_options[dataset.name_ref]['sample_factor'])

    def compute(self, parameter_values, dataset, x0_input=None):
        """
        :param parameter_values:
        :param dataset:
        :param x0_input:
        :return:
        """
        self.update_parameter_values(parameter_values, dataset.Tref)

        for par, i_par in self.ldvars.items():
            self.ld_vars[i_par] = parameter_values[par]

        if x0_input is None:
            return self.pytransit_models[dataset.name_ref].evaluate(
                parameter_values['R_Rs'],
                self.ld_vars,
                parameter_values['Tc'] - dataset.Tref,
                parameter_values['P'],
                parameter_values['a_Rs'],
                parameter_values['i'] * constants.deg2rad,
                parameter_values['e'],
                parameter_values['omega'] * constants.deg2rad) - 1.

        else:
            self.pytransit_plot[dataset.name_ref].set_data(x0_input,
                                                            exptimes=self.code_options[dataset.name_ref]['exp_time'],
                                                            nsamples=self.code_options[dataset.name_ref]['sample_factor'])

            return self.pytransit_plot[dataset.name_ref].evaluate(
                parameter_values['R_Rs'],
                self.ld_vars,
                parameter_values['Tc'] - dataset.Tref,
                parameter_values['P'],
                parameter_values['a_Rs'],
                parameter_values['i'],
                parameter_values['e'],
                parameter_values['omega'] * constants.deg2rad) - 1.