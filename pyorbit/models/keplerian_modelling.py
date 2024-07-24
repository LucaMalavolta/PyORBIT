from __future__ import print_function
from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import AbstractModel

import pyorbit.subroutines.kepler_exo as kepler_exo


class RVkeplerian(AbstractModel):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'radial_velocities'

        self.list_pams_common = OrderedSet([
            'P',  # Period
            'K',  # RV semi-amplitude
            'e',  # eccentricity, uniform prior - to be fixed
            'omega'  # argument of pericenter
        ])

        self.use_time_inferior_conjunction = False
        self.use_mass = False

    def initialize_model(self, mc, **kwargs):

        if mc.common_models[self.planet_ref].parametrization[:8] == 'Ford2006' \
            and mc.common_models[self.planet_ref].orbit != 'circular':
            self.list_pams_common.discard('e')
            self.list_pams_common.discard('omega')

            self.list_pams_common.update(['e_coso'])
            self.list_pams_common.update(['e_sino'])

        elif mc.common_models[self.planet_ref].parametrization[:8] != 'Standard' \
            and mc.common_models[self.planet_ref].orbit != 'circular':
                # 'Eastman2013' is the standard choice
            self.list_pams_common.discard('e')
            self.list_pams_common.discard('omega')

            self.list_pams_common.update(['sre_coso'])
            self.list_pams_common.update(['sre_sino'])

        if mc.common_models[self.planet_ref].use_time_inferior_conjunction:
            self.list_pams_common.update(['Tc'])
            self.use_time_inferior_conjunction = True
            # Copying the property to the class for faster access
        else:
            self.list_pams_common.update(['mean_long'])

        if mc.common_models[self.planet_ref].use_mass:
            self.list_pams_common.update(['M_Me'])
            self.list_pams_common.update(['mass'])
            self.use_mass = True
        else:
            self.list_pams_common.update(['K'])

    def compute(self, parameter_values, dataset, x0_input=None):

        #if parameter_values['P']<10:
        #    print(parameter_values['P'], parameter_values['K'], parameter_values['e'], parameter_values['omega'])

        if self.use_time_inferior_conjunction:
            mean_long = kepler_exo.kepler_Tc2phase_Tref(parameter_values['P'],
                                                parameter_values['Tc'] - dataset.Tref,
                                                parameter_values['e'],
                                                parameter_values['omega'])
        else:
            mean_long = parameter_values['mean_long']

        if self.use_mass:

            K = kepler_exo.kepler_K1(parameter_values['mass'],
                                     parameter_values['M_Me'] / constants.Msear, parameter_values['P'], parameter_values['i'],
                                     parameter_values['e'])
        else:
            K = parameter_values['K']

        if x0_input is None:
            return kepler_exo.kepler_RV_T0P(dataset.x0,
                                            mean_long,
                                            parameter_values['P'],
                                            K,
                                            parameter_values['e'],
                                            parameter_values['omega'])
        else:
            return kepler_exo.kepler_RV_T0P(x0_input,
                                            mean_long,
                                            parameter_values['P'],
                                            K,
                                            parameter_values['e'],
                                            parameter_values['omega'])


class TransitTimeKeplerian(AbstractModel):
    model_class = 'transit_times'

    def __init__(self, *args, **kwargs):
        super(TransitTimeKeplerian, self).__init__(*args, **kwargs)
        self.use_time_inferior_conjunction = False

        self.list_pams_common = OrderedSet(['P'])  # Period

        self.list_pams_dataset = OrderedSet()

    def initialize_model(self, mc, **kwargs):

        if mc.common_models[self.planet_ref].use_time_inferior_conjunction:
            self.list_pams_common.update(['Tc'])
            self.use_time_inferior_conjunction = True
            # Copying the property to the class for faster access
        else:
            self.list_pams_common.update(['mean_long'])
            self.list_pams_common.update(['e'])
            self.list_pams_common.update(['omega'])
            # mean longitude = argument of pericenter + mean anomaly at Tref

    def compute(self, parameter_values, dataset, x0_input=None):

        if self.use_time_inferior_conjunction:
            delta_T = parameter_values['Tc'] - \
                np.floor((parameter_values['Tc'] - dataset.Tref) / parameter_values['P']) * parameter_values['P']
        else:
            delta_T = dataset.Tref + \
                      kepler_exo.kepler_phase2Tc_Tref(parameter_values['P'],
                                                      parameter_values['mean_long'],
                                                      parameter_values['e'],
                                                      parameter_values['omega'])

        if x0_input is None:
            return np.floor(dataset.x0 / parameter_values['P']) * parameter_values['P'] + delta_T
        else:
            return np.floor(x0_input / parameter_values['P']) * parameter_values['P'] + delta_T


