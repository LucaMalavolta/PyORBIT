from __future__ import print_function
from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import AbstractModel

import pyorbit.subroutines.kepler_exo as kepler_exo

from scipy.stats import norm

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

        if mc.common_models[self.planet_ref].use_longitude_of_nodes:
            self.list_pams_common.update(['Omega'])

        try:
            self.default_Omega = mc.common_models[self.planet_ref].default_omega
        except AttributeError:
            self.default_Omega = 0.00

    def compute(self, parameter_values, dataset, x0_input=None):

        Omega = parameter_values.get('Omega', self.default_Omega)

        if self.use_time_inferior_conjunction:
            mean_long = kepler_exo.kepler_compute_meanlong_from_deltaTc(parameter_values['P'],
                                                parameter_values['Tc'] - dataset.Tref,
                                                parameter_values['e'],
                                                parameter_values['omega'],
                                                Omega)
        else:
            mean_long = parameter_values['mean_long']

        if self.use_mass:

            rv_semiamplitude = kepler_exo.kepler_compute_rv_semiamplitude(parameter_values['mass'],
                                                                            parameter_values['M_Me'] / constants.Msear,
                                                                            parameter_values['P'], 
                                                                            parameter_values['i'],
                                                                            parameter_values['e'])
        else:
            rv_semiamplitude = parameter_values['K']

        if x0_input is None:
            return kepler_exo.kepler_compute_rv_deltabjd(dataset.x0,
                                            rv_semiamplitude,
                                            parameter_values['P'],
                                            mean_long,
                                            parameter_values['e'],
                                            parameter_values['omega'],
                                            Omega)
        else:
            return kepler_exo.kepler_compute_rv_deltabjd(x0_input,
                                            rv_semiamplitude,
                                            parameter_values['P'],
                                            mean_long,
                                            parameter_values['e'],
                                            parameter_values['omega'],
                                            Omega)



class ApodizedRVkeplerian(RVkeplerian):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.list_pams_common.update(['apo_center'])
        self.list_pams_common.update(['apo_timescale'])

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        try:
            min_apo_center = min(mc.common_models[self.planet_ref].default_bounds['apo_center'][0],
                                    np.min(dataset.x))
            max_apo_center = max(mc.common_models[self.planet_ref].default_bounds['apo_center'][1],
                                    np.max(dataset.x))
        except KeyError:
            min_apo_center = np.min(dataset.x)
            max_apo_center = np.max(dataset.x)

        mc.common_models[self.planet_ref].default_bounds['apo_center'] = [min_apo_center, max_apo_center]
        return

    def compute(self, parameter_values, dataset, x0_input=None):

        Omega = parameter_values.get('Omega', self.default_Omega)

        if self.use_time_inferior_conjunction:
            mean_long = kepler_exo.kepler_compute_meanlong_from_deltaTc(parameter_values['P'],
                                                parameter_values['Tc'] - dataset.Tref,
                                                parameter_values['e'],
                                                parameter_values['omega'],
                                                Omega)
        else:
            mean_long = parameter_values['mean_long']

        if self.use_mass:

            rv_semiamplitude = kepler_exo.kepler_compute_rv_semiamplitude(parameter_values['mass'],
                                                                            parameter_values['M_Me'] / constants.Msear,
                                                                            parameter_values['P'],
                                                                            parameter_values['i'],
                                                                            parameter_values['e'])
        else:
            rv_semiamplitude = parameter_values['K']

        if x0_input is None:

            # Apodization factor
            apo_function = (norm.pdf(dataset.x,
                                        parameter_values['apo_center'],
                                        parameter_values['apo_timescale']) 
                            / norm.pdf(parameter_values['apo_center'], 
                                        parameter_values['apo_center'],
                                        parameter_values['apo_timescale']))

            return apo_function * kepler_exo.kepler_compute_rv_deltabjd(dataset.x0,
                                            rv_semiamplitude,
                                            parameter_values['P'],
                                            mean_long,
                                            parameter_values['e'],
                                            parameter_values['omega'],
                                            Omega)
        else:

            # Apodization factor
            apo_function = (norm.pdf(dataset.x,
                                        parameter_values['apo_center'],
                                        parameter_values['apo_timescale']) 
                            / norm.pdf(parameter_values['apo_center'], 
                                        parameter_values['apo_center'],
                                        parameter_values['apo_timescale']))

            return apo_function * kepler_exo.kepler_compute_rv_deltabjd(dataset.x0,
                                            rv_semiamplitude,
                                            parameter_values['P'],
                                            mean_long,
                                            parameter_values['e'],
                                            parameter_values['omega'],
                                            Omega)


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

        try:
            self.default_Omega = mc.common_models[self.planet_ref].default_omega
        except AttributeError:
            self.default_Omega = 0.00


    def compute(self, parameter_values, dataset, x0_input=None):


        #TODO check this function, I believe it's a bit crazy.... URGENT
        Omega = parameter_values.get('Omega', self.default_Omega)

        if self.use_time_inferior_conjunction:
            delta_T = parameter_values['Tc'] - \
                np.floor((parameter_values['Tc'] - dataset.Tref) / parameter_values['P']) * parameter_values['P']
        else:
            delta_T = dataset.Tref + \
                      kepler_exo.kepler_compute_deltaTc_from_meanlong(parameter_values['P'],
                                                      parameter_values['mean_long'],
                                                      parameter_values['e'],
                                                      parameter_values['omega'],
                                                      Omega)

        if x0_input is None:
            return np.floor(dataset.x0 / parameter_values['P']) * parameter_values['P'] + delta_T
        else:
            return np.floor(x0_input / parameter_values['P']) * parameter_values['P'] + delta_T


