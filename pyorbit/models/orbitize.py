
from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import AbstractModel
import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.subroutines.transformations import *

class AbstractAstrometry(object):

    def __init__(self, *args, **kwargs):

        ''' Orbital parameters to be used in the astrometric fit '''
        self.list_pams_common = OrderedSet([
            'P',     # Period in days
            'Omega', # longitude of ascending node
            'e',     # eccentricity, uniform prior - to be fixed
            'R_Rs',  # planet radius (in units of stellar radii)
            'omega', # argument of pericenter
            #'i',     # inclination in degrees
            'mass', #stellar mass
            'parallax', #stellar parallax
        ])

        self.list_pams_dataset = OrderedSet()
        self.warning_given = False


    def _prepare_astrometry_parameters(self, mc, **kwargs):


        print("*** {0:s} global parameters:".format(self.model_name))

        _prepare_planet_parametrization(self, mc, **kwargs)
        _prepare_planet_scaled_semimajor_axis(self, mc, **kwargs)
        _prepare_planet_mass(self, mc, **kwargs)
        _prepare_stellar_mass(self, mc, **kwargs)
        _prepare_planet_inclination(self, mc, **kwargs)
        _prepare_planet_time_inferior_conjunction(self, mc, **kwargs)

    def update_parameter_values_for_astrometry(self, parameter_values, Tref, prepend=''):

        if self.compute_inclination:
            if self.compute_scaled_semimajor_axis:
                parameter_values[prepend+'a_Rs'] = convert_rho_to_ars(parameter_values[prepend+'P'], parameter_values['density'])
            parameter_values[prepend+'i'] = convert_b_to_i(
                parameter_values[prepend+'b'], parameter_values[prepend+'e'], parameter_values[prepend+'omega'], parameter_values[prepend+'a_Rs'])

        if self.compute_semimajor_axis:
            if self.compute_scaled_semimajor_axis:
                parameter_values[prepend+'a_Rs'] = convert_rho_to_ars(parameter_values[prepend+'P'], parameter_values['density'])
            parameter_values[prepend+'a_AU'] = convert_ars_to_a(parameter_values[prepend+'a_Rs'], parameter_values['radius'])

        if self.compute_time_inferior_conjunction:
            parameter_values[prepend+'Tc']= kepler_exo.kepler_compute_deltaTc_from_meanlong(
                parameter_values[prepend+'P'],
                parameter_values[prepend+'mean_long'],
                parameter_values[prepend+'e'],
                parameter_values[prepend+'omega'],
                parameter_values[prepend+'Omega']) + Tref

        if self.compute_mean_longitude:
            parameter_values[prepend+'mean_long'] = kepler_exo.kepler_compute_meanlong_from_deltaTc(
                parameter_values[prepend+'P'],
                parameter_values[prepend+'Tc'] - Tref,
                parameter_values[prepend+'e'],
                parameter_values[prepend+'omega'],
                parameter_values[prepend+'Omega'])


class Orbitize(AbstractModel, AbstractAstrometry):
    model_class = 'orbitize_model'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

    # brainless workaround
    def initialize_model(self, mc, **kwargs):
        self._prepare_astrometry_parameters(mc, **kwargs)

        print(**kwargs)

    def compute(self, mc, theta, x_input=None, *args, **kwargs):
        print(self.multiple_planets)
        for planet_name in self.multiple_planets:
            self.update_parameter_values(parameter_values, dataset.Tref, planet_name+'_' )


class OrbitizeRunner(object):
    def __init__(self):
        self.model_name = 'orbitize_runner'
        self.to_be_initialized = True

        print("    {0:s} WARNING:".format(self.model_name))
        print('        Astrometry modelling requires the use of the stellar mass')
        print('        This may cause a clash with models requiring stellar density and radius, e.g., RM modelling')
        print('        The use of a multivariate approach is strongly suggested')
        print('        You can control the behaviour of mass/radius/density with the specific keywords')
        print('        compute_mass, compute_radius, compute_density')
        print()

    #def compute(self, mc, theta, x_input=None, *args, **kwargs):
        
