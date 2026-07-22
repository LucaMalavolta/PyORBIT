
from os import system

from pyorbit.subroutines.common import *
import pyorbit.subroutines.kepler_exo as kepler_exo
import pyorbit.subroutines.transformations as transformations
import pyorbit.subroutines.constants as constants
import os   



class AbstractAstrometry(object):

    def __init__(self, *args, **kwargs):

        self.compute_semimajor_axis = True
        self.compute_scaled_semimajor_axis = True
        self.compute_semimajor_axis_from_mass = False

    def _prepare_astrometry_parameters(self, mc, **kwargs):


        print("*** {0:s} global parameters:".format(self.model_name))

        transformations._prepare_planet_parametrization(self, mc, **kwargs)
        transformations._prepare_planet_scaled_semimajor_axis(self, mc, **kwargs)
        transformations._prepare_planet_semimajor_axis(self, mc, **kwargs)
        transformations._prepare_planet_mass(self, mc, **kwargs)
        transformations._prepare_stellar_mass(self, mc, **kwargs)
        transformations._prepare_planet_inclination(self, mc, **kwargs)
        transformations._prepare_planet_time_inferior_conjunction(self, mc, **kwargs)
        print()

    def update_parameter_values(self,
                                parameter_values,
                                input_prepend=''):

        if input_prepend == '':
            prepend = ''
        else:
            prepend = input_prepend + '__'

        if self.compute_planet_mass:
            sin_i = abs(np.sin(parameter_values[prepend+"i"]*constants.deg2rad))
            parameter_values[prepend+'M_Me'] = kepler_exo.kepler_get_planet_mass(
                parameter_values[prepend+'P'],
                parameter_values[prepend+'K']/sin_i,
                parameter_values[prepend+'e'],
                parameter_values['mass'],
                approximation_limit=10,
                verbose=False)

        if self.compute_semimajor_axis_from_mass and self.compute_scaled_semimajor_axis:
            parameter_values[prepend+'a_AU']= convert_PMsMp_to_a(
                    parameter_values[prepend+'P'],
                    parameter_values['mass'],
                    parameter_values[prepend+'M_Me'])
            parameter_values[prepend+'a_Rs'] = convert_a_to_ars(parameter_values[prepend+'a_AU'], parameter_values['radius'])
            if self.compute_inclination:
                parameter_values[prepend+'i'] = convert_b_to_i(
                    parameter_values[prepend+'b'], parameter_values[prepend+'e'], parameter_values[prepend+'omega'], parameter_values[prepend+'a_Rs'])
        else:
            if self.compute_inclination:
                if self.compute_scaled_semimajor_axis:
                    parameter_values[prepend+'a_Rs'] = convert_rho_to_ars(parameter_values[prepend+'P'], parameter_values['density'])
                parameter_values[prepend+'i'] = convert_b_to_i(
                    parameter_values[prepend+'b'], parameter_values[prepend+'e'], parameter_values[prepend+'omega'], parameter_values[prepend+'a_Rs'])

            if self.compute_semimajor_axis:
                if self.compute_scaled_semimajor_axis:
                    parameter_values[prepend+'a_Rs'] = convert_rho_to_ars(parameter_values[prepend+'P'], parameter_values['density'])
                parameter_values[prepend+'a_AU'] = convert_ars_to_a(parameter_values[prepend+'a_Rs'], parameter_values['radius'])

            sin_i = abs(np.sin(parameter_values[prepend+"i"]*constants.deg2rad))
            parameter_values[prepend+'M_Me'] = kepler_exo.kepler_get_planet_mass(
                parameter_values[prepend+'P'],
                parameter_values[prepend+'K']/sin_i,
                parameter_values[prepend+'e'],
                parameter_values['mass'],
                approximation_limit=10,
                verbose=False)

        if self.compute_time_inferior_conjunction:
            parameter_values[prepend+'Tc']= kepler_exo.kepler_compute_deltaTc_from_meanlong(
                parameter_values[prepend+'P'],
                parameter_values[prepend+'mean_long'],
                parameter_values[prepend+'e'],
                parameter_values[prepend+'omega'],
                parameter_values[prepend+'Omega']) + self.Tref

        if self.compute_mean_longitude:
            parameter_values[prepend+'mean_long'] = kepler_exo.kepler_compute_meanlong_from_deltaTc(
                parameter_values[prepend+'P'],
                parameter_values[prepend+'Tc'] - self.Tref,
                parameter_values[prepend+'e'],
                parameter_values[prepend+'omega'],
                parameter_values[prepend+'Omega'])

            parameter_values[prepend+'Tperi'] = kepler_exo.kepler_compute_deltaTperi_from_deltaTc(
                    parameter_values[prepend+'P'],
                    parameter_values[prepend+'Tc'] - self.Tref,
                    parameter_values[prepend+'e'],
                    parameter_values[prepend+'omega'])
        else:
            parameter_values[prepend+'Tperi'] = kepler_exo.kepler_compute_deltaTperi_from_meanlong(
                    parameter_values[prepend+'P'],
                    parameter_values[prepend+'mean_long'],
                    parameter_values[prepend+'e'],
                    parameter_values[prepend+'omega'],
                    parameter_values[prepend+'Omega'])

