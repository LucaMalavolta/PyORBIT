from pyorbit.subroutines.common import np, convert_rho_to_a, convert_b_to_i
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo


class AbstractTransit(object):

    def __init__(self, *args, **kwargs):

        self.model_class = 'transit'
        self.unitary_model = True

        self.ldvars = {}
        self.ld_ncoeff = 2
        self.parametrization = 'Standard'

        self.use_semimajor_axis = False
        self.use_inclination = False
        self.use_time_of_transit = False
        self.use_stellar_radius = False
        self.use_stellar_temperature = False

        self.limb_darkening_model = None

        self.retrieve_ai = None
        self.retrieve_t0 = None
        self.retrieve_radius = None
        self.retrieve_temperature = None

        self.multivariate_mass_radius = False
        self.code_options = {}

    def _prepare_planetary_parameters(self, mc, **kwargs):

        try:
            multivariate_vars = mc.common_models[self.stellar_ref].multivariate_vars
        except AttributeError:
            multivariate_vars = []

        """ Default parametrization uses the stellar density and the impact
            parameter, it is possible to switch back to scaled semi-major axis and
            inclination respectively by activating the proper flag """

        if mc.common_models[self.planet_ref].use_semimajor_axis:
            """ a is the semi-major axis (in units of stellar radii) """
            self.list_pams_common.update({'a': None})
            self.use_semimajor_axis = True
        else:
            if 'mass' in multivariate_vars and 'radius' in multivariate_vars:
                self.list_pams_common.update({'mass': None, 'radius':None})
                self.multivariate_mass_radius = True
            else:
                """ rho is the density of the star (in solar units) """
                self.list_pams_common.update({'rho': None})
                self.multivariate_mass_radius = False

        if mc.common_models[self.planet_ref].use_inclination:
            """ i is the orbital inclination (in degrees) """
            self.list_pams_common.update({'i': None})
            self.use_inclination = True
        else:
            """ b is the impact parameter """
            self.list_pams_common.update({'b': None})

        if mc.common_models[self.planet_ref].use_time_of_transit:
            self.list_pams_common.update({'Tc': None})
            self.use_time_of_transit = True
            # Copying the property to the class for faster access
        else:
            self.list_pams_common.update({'f': None})
            # mean longitude = argument of pericenter + mean anomaly at Tref

        """ The appropriate function for variable conversion is stored internally
        """
        if self.use_semimajor_axis and self.use_inclination:
            self.retrieve_ai = self._internal_transformation_mod03
        elif self.use_semimajor_axis:
            if self.multivariate_mass_radius:
                self.retrieve_ai = self._internal_transformation_mod07
            else:
                self.retrieve_ai = self._internal_transformation_mod02
        elif self.use_inclination:
            self.retrieve_ai = self._internal_transformation_mod01
        else:
            if self.multivariate_mass_radius:
                self.retrieve_ai = self._internal_transformation_mod06
            else:
                self.retrieve_ai = self._internal_transformation_mod00

        if self.use_time_of_transit:
            self.retrieve_t0 = self._internal_transformation_mod04
        else:
            self.retrieve_t0 = self._internal_transformation_mod05

    def _prepare_star_parameters(self, mc, **kwargs):
        """ Additional stellar parameters,
            check if the stellar radius and effect temperature are provided as
            fixed values or not.
            As in version 9., only used by spiderman models
        """
        self.use_stellar_radius = True

        stellarradius_names = [
            'stellar_radius',
            'radius',
            'star_radius'
        ]

        for dict_name in stellarradius_names:
            if kwargs.get(dict_name, False):
                self.code_options['radius'] = kwargs[dict_name]
                self.use_stellar_radius = False
                self.retrieve_radius = _internal_transformation_mod10
        if self.use_stellar_radius:
            self.list_pams_common.update({'radius': None})
            self.retrieve_radius = _internal_transformation_mod11

        self.use_stellar_temperature = True

        effectivetemperature_names = [
            'teff',
            'temperature',
            'eff_temperature'
            'temperature_eff'
        ]
        for dict_name in effectivetemperature_names:
            if kwargs.get(dict_name, False):
                self.code_options['temperature'] = kwargs[dict_name]
                self.use_stellar_temperature = False
                self.retrieve_temperature = _internal_transformation_mod12
        if self.use_stellar_temperature:
            self.list_pams_common.update({'temperature': None})
            self.retrieve_temperature = _internal_transformation_mod13

    def _prepare_limnb_darkening_coefficients(self, mc, **kwargs):
        """ Setting up the limb darkening calculation"""

        self.limb_darkening_model = kwargs['limb_darkening_model']
        self.ld_vars = [0.00] * kwargs['limb_darkening_ncoeff']
        for i_coeff in range(1, kwargs['limb_darkening_ncoeff'] + 1):
            var = 'ld_c' + repr(i_coeff)
            self.ldvars[var] = i_coeff - 1
            self.list_pams_common.update({var: None})

    def _prepare_dataset_options(self, mc, dataset, **kwargs):

        supersample_names = ['supersample_factor',
                             'supersample',
                             'supersampling',
                             'oversample_factor',
                             'oversample',
                             'oversampling',
                             'sample_factor',
                             'sample',
                             'sampling'
                             'nsample_factor',
                             'nsample',
                             'nsampling'
                             ]

        sample_factor = 1
        exposure_time = 0.01

        for dict_name in supersample_names:
            if kwargs[dataset.name_ref].get(dict_name, False):
                sample_factor = kwargs[dataset.name_ref][dict_name]
            elif kwargs.get(dict_name, False):
                sample_factor = kwargs[dict_name]

        exptime_names = ['exposure_time',
                         'exposure',
                         'exp_time',
                         'exptime',
                         'obs_duration',
                         'integration',
                         ]

        for dict_name in exptime_names:
            if kwargs[dataset.name_ref].get(dict_name, False):
                exposure_time = kwargs[dataset.name_ref][dict_name]
            elif kwargs.get(dict_name, False):
                exposure_time = kwargs[dict_name]

        self.code_options[dataset.name_ref] = {
            'sample_factor': sample_factor,
            'exp_time': exposure_time / constants.d2s,
        }

        wavebounds_names = [
            'wavelength_range',
            'wavelength_boundaries',
        ]

        """ Lower and upper wavelength boundaries for the filter, in (nm) """
        for dict_name in wavebounds_names:
            if kwargs[dataset.name_ref].get(dict_name, False):
                self.code_options[dataset.name_ref]['l1'] = kwargs[dataset.name_ref][dict_name][0] / 10**9
                self.code_options[dataset.name_ref]['l2'] = kwargs[dataset.name_ref][dict_name][1] / 10**9
            elif kwargs.get(dict_name, False):
                self.code_options[dataset.name_ref]['l1'] = kwargs[dict_name][0] / 10**9
                self.code_options[dataset.name_ref]['l2'] = kwargs[dict_name][1] / 10**9


    """ function for internal transformation of variables, to avoid if calls"""
    def _internal_transformation_mod00(self, variable_value):
        """ this function transforms b and rho to i and a  """
        a = convert_rho_to_a(variable_value['P'], variable_value['rho'])
        i = convert_b_to_i(
            variable_value['b'], variable_value['e'], variable_value['o'], a)
        return a, i

    def _internal_transformation_mod01(self, variable_value):
        """ this function transforms b to i"""
        i = convert_b_to_i(
            variable_value['b'], variable_value['e'], variable_value['o'], variable_value['a'])
        return variable_value['a'], i

    def _internal_transformation_mod02(self, variable_value):
        """ this function transforms rho to a  """
        a = convert_rho_to_a(variable_value['P'], variable_value['rho'])
        return a, variable_value['i']

    def _internal_transformation_mod03(self, variable_value):
        """ no transformation needed  """
        return variable_value['a'], variable_value['i']

    def _internal_transformation_mod04(self, variable_value, Tref):
        """ this function transforms Tc into Tc- Tref t"""
        return variable_value['Tc'] - Tref

    def _internal_transformation_mod05(self, variable_value, Tref):
        """ this function transforms Tc into Tc- Tref t"""
        return kepler_exo.kepler_phase2Tc_Tref(variable_value['P'],
                                               variable_value['f'],
                                               variable_value['e'],
                                               variable_value['o'])

    def _internal_transformation_mod06(self, variable_value):
        """ this function transforms b, mass, radius to i and a 
            it replaces _internal_transformation_mod00 when mass & radius
            multivariate are used
        """
        rho = variable_value['mass']/variable_value['radius']**3
        a = convert_rho_to_a(variable_value['P'], rho)
        i = convert_b_to_i(
            variable_value['b'], variable_value['e'], variable_value['o'], a)
        return a, i

    def _internal_transformation_mod07(self, variable_value):
        """ this function transforms P,mass, radius to a
            it replaces _internal_transformation_mod02 when mass & radius
            multivariate are used
        """
        rho = variable_value['mass']/variable_value['radius']**3
        a = convert_rho_to_a(variable_value['P'], rho)
        return a, variable_value['i']

    def _internal_transformation_mod10(self, variable_value):
        return self.code_options['radius']

    def _internal_transformation_mod11(self, variable_value):
        return variable_value['radius']

    def _internal_transformation_mod12(self, variable_value):
        return self.code_options['temperature']

    def _internal_transformation_mod13(self, variable_value):
        return variable_value['temperature']