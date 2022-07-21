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
        self.use_stellar_period = False
        self.stellar_period_from_activity = False

        self.limb_darkening_model = None

        self.retrieve_ai = None
        self.retrieve_t0 = None
        self.retrieve_radius = None
        self.retrieve_temperature = None
        self.retrieve_Omega_Istar = None

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
            self.list_pams_common.update(['a_Rs'])
            self.use_semimajor_axis = True
        else:
            if 'mass' in multivariate_vars and 'radius' in multivariate_vars:
                self.list_pams_common.update(['mass'])
                self.list_pams_common.update(['radius'])
                self.multivariate_mass_radius = True
            else:
                """ rho is the density of the star (in solar units) """
                self.list_pams_common.update(['density'])
                self.multivariate_mass_radius = False

        if mc.common_models[self.planet_ref].use_inclination:
            """ i is the orbital inclination (in degrees) """
            self.list_pams_common.update(['i'])
            self.use_inclination = True
        else:
            """ b is the impact parameter """
            self.list_pams_common.update(['b'])

        if mc.common_models[self.planet_ref].use_time_of_transit:
            self.list_pams_common.update(['Tc'])
            self.use_time_of_transit = True
            # Copying the property to the class for faster access
        else:
            self.list_pams_common.update(['mean_long'])
            # mean longitude = argument of pericenter + mean anomaly at Tref

        """ The appropriate function for variable conversion is stored internally
        """
        if self.use_semimajor_axis and self.use_inclination:
            self.retrieve_ai = self._internal_transformation_mod03
        elif self.use_inclination:
            if self.multivariate_mass_radius:
                self.retrieve_ai = self._internal_transformation_mod07
            else:
                self.retrieve_ai = self._internal_transformation_mod02
        elif self.use_semimajor_axis:
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

        stellarradius_names = [
            'stellar_radius',
            'radius',
            'star_radius'
        ]

        for dict_name in stellarradius_names:
            if kwargs.get(dict_name, False):
                self.code_options['radius'] = kwargs[dict_name]
                self.use_stellar_radius = False
                self.retrieve_radius = self._internal_transformation_mod10
        if self.use_stellar_radius:
            self.list_pams_common.update(['radius'])
            self.retrieve_radius = self._internal_transformation_mod11

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
                self.retrieve_temperature = self._internal_transformation_mod12
        if self.use_stellar_temperature:
            self.list_pams_common.update(['temperature'])
            self.retrieve_temperature = self._internal_transformation_mod13


        """ check if the stellar_rotation has to be used as parameter """
        self.use_stellar_period = kwargs.get('use_stellar_period', self.use_stellar_period)


        """ check if the stellar_inclination has to be used as parameter """
        self.use_stellar_inclination = kwargs.get('use_stellar_inclination', self.use_stellar_inclination)


        """ Check if the stellar rotation period is given as a starting value
            If so, the angular rotation of the star and the stellar inclination
            are computed through the rotation period
        """

        rotationperiod_names =[
            'rotation_period',
            'rotational_period',
            'stellar_period',
            'stellar_rotation',
            'Prot',
        ]
        for dict_name in rotationperiod_names:
            if kwargs.get(dict_name, False):
                self.code_options['rotation_period'] = kwargs[dict_name]
                self.use_stellar_period = False
                self.use_stellar_inclination = False
                self.retrieve_Omega_Istar = self._internal_transformation_mod23

        """ Check if the activity model is included in the keyword list
            If so, the angular rotation of the star and the stellar inclination
            are computed through the rotation period
        """
        rotationperiod_from_activity_names = [
            'activity',
            'activity_model',
            'activity_rotation',
            'stellar_rotation_from_activity',
            'stellar_period_from_activity',
            'Prot_from_activity',
            'star_rotation_from_activity',
            'star_period_from_activity',
            'rotation_period_from_activity'
        ]

        for dict_name in rotationperiod_from_activity_names:
            if kwargs.get(dict_name, False):
                self.use_stellar_period = False
                self.use_stellar_inclination = False
                self.list_pams_common.update(['Prot'])
                self.common_ref.append(kwargs[dict_name])
                self.retrieve_Omega_Istar = self._internal_transformation_mod22

        if self.use_stellar_period and self.use_stellar_inclination:
            print('Error in Stellar model  inside Transit model definition ')
            print('It is not possibile to have both stellar_period and stellar_inclination as free parameters')
            quit()

        if self.use_stellar_period:
                self.list_pams_common.update(['rotation_period'])
                self.retrieve_Omega_Istar = self._internal_transformation_mod21

        if self.use_stellar_inclination:
                self.list_pams_common.update(['i_star'])
                self.retrieve_Omega_Istar = self._internal_transformation_mod20


    def _prepare_limb_darkening_coefficients(self, mc, **kwargs):
        """ Setting up the limb darkening calculation"""

        self.limb_darkening_model = kwargs['limb_darkening_model']
        self.ld_vars = [0.00] * kwargs['limb_darkening_ncoeff']
        for i_coeff in range(1, kwargs['limb_darkening_ncoeff'] + 1):
            var = 'ld_c' + repr(i_coeff)
            self.ldvars[var] = i_coeff - 1
            self.list_pams_common.update([var])

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

        """ Keep track of the boundaries of each dataset, so that the user do
            not have to write down the boundaries of each transit in case of TTV fit
        """

        tc_boundaries_names = ['Tc_boundaries',
            'Tc_bounds',
            'T_boundaries',
            'T_bounds',
            'TC_boundaries',
            'TC_bounds',
        ]

        tc_boundaries = [np.amin(dataset.x), np.amax(dataset.x)]
        if kwargs[dataset.name_ref].get('boundaries', False):
            tc_boundaries = kwargs[dataset.name_ref]['boundaries']['Tc']
        else:
            for dict_name in tc_boundaries_names:
                if kwargs[dataset.name_ref].get(dict_name, False):
                    tc_boundaries = kwargs[dataset.name_ref][dict_name]

        self.code_options[dataset.name_ref]['Tc_boundaries'] = tc_boundaries


    """ function for internal transformation of variables, to avoid if calls"""
    @staticmethod
    def _internal_transformation_mod00(variable_value):
        """ this function transforms b and stellar density to i and a  """
        a = convert_rho_to_a(variable_value['P'], variable_value['density'])
        i = convert_b_to_i(
            variable_value['b'], variable_value['e'], variable_value['omega'], a)
        return a, i

    @staticmethod
    def _internal_transformation_mod01(variable_value):
        """ this function transforms b to i"""
        i = convert_b_to_i(
            variable_value['b'], variable_value['e'], variable_value['omega'], variable_value['a_Rs'])
        return variable_value['a_Rs'], i

    @staticmethod
    def _internal_transformation_mod02(variable_value):
        """ this function transforms stellar density to a  """
        a = convert_rho_to_a(variable_value['P'], variable_value['density'])
        return a, variable_value['i']

    @staticmethod
    def _internal_transformation_mod03(variable_value):
        """ no transformation needed  """
        return variable_value['a_Rs'], variable_value['i']

    @staticmethod
    def _internal_transformation_mod04(variable_value, Tref):
        """ this function transforms Tc into Tc- Tref t"""
        return variable_value['Tc'] - Tref

    @staticmethod
    def _internal_transformation_mod05(variable_value, Tref):
        """ this function transforms phase into Tc- Tref t"""
        return kepler_exo.kepler_phase2Tc_Tref(variable_value['P'],
                                               variable_value['mean_long'],
                                               variable_value['e'],
                                               variable_value['omega'])

    @staticmethod
    def _internal_transformation_mod06(variable_value):
        """ this function transforms b, mass, radius to i and a
            it replaces _internal_transformation_mod00 when mass & radius
            multivariate are used
        """
        rho = variable_value['mass']/variable_value['radius']**3
        a = convert_rho_to_a(variable_value['P'], rho)
        i = convert_b_to_i(
            variable_value['b'], variable_value['e'], variable_value['omega'], a)
        return a, i

    @staticmethod
    def _internal_transformation_mod07(variable_value):
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

    @staticmethod
    def _internal_transformation_mod20(variable_value):
        """ This function extracts from the variable_value dictionary
            - 'v_sini': projected tangential velocity of the star
            - 'radius': stellar radius in Solar units
            - 'i_star': inclination of the of stellar rotation axis

            into angular velocity Omega [rad/s] of the star
            it returns the angular velocity Omega and the stellar inclination

        Args:
            variable_value (dict): dictionary with all the parameters of the current model

        Returns:
            float: angular velocity Omega in [rad/s]
            float: inclination of the of stellar rotation axis, [degrees]
        """

        Omega = variable_value['v_sini'] / (variable_value['radius'] * constants.Rsun) / np.sin(variable_value['i_star']/180.*np.pi)

        return Omega, variable_value['i_star']

    @staticmethod
    def _internal_transformation_mod21(variable_value):
        """ This function extracts from the variable_value dictionary
            - 'v_sini': projected tangential velocity of the star
            - 'radius': stellar radius in Solar units
            - 'rotation_period': rotational period of the star from the stellar model

            into angular velocity Omega [rad/s] of the star
            it returns the angular velocity Omega and the stellar inclination

        Args:
            variable_value (dict): dictionary with all the parameters of the current model

        Returns:
            float: angular velocity Omega, [rad/s]
            float: inclination of the of stellar rotation axis, [degrees]
        """

        sin_is =  variable_value['v_sini']  / (variable_value['radius'] * constants.Rsun) * (variable_value['rotation_period'] * constants.d2s)  / (2* np.pi)
        Omega = 2* np.pi / (variable_value['Prot'] * constants.d2s)

        return Omega, np.arcsin(sin_is)/np.pi*180.

    @staticmethod
    def _internal_transformation_mod22(variable_value):
        """ This function extracts from the variable_value dictionary
            - 'v_sini': projected tangential velocity of the star
            - 'radius': stellar radius in Solar units
            - 'Prot': rotational period of the star from the activity model

            into angular velocity Omega [rad/s] of the star
            it returns the angular velocity Omega and the stellar inclination

        Args:
            variable_value (dict): dictionary with all the parameters of the current model

        Returns:
            float: angular velocity Omega, [rad/s]
            float: inclination of the of stellar rotation axis, [degrees]
        """

        sin_is =  variable_value['v_sini']  / (variable_value['radius'] * constants.Rsun) * (variable_value['Prot'] * constants.d2s)  / (2* np.pi)
        Omega = 2* np.pi / (variable_value['Prot'] * constants.d2s)

        return Omega, np.arcsin(sin_is)/np.pi*180.

    def _internal_transformation_mod23(self, variable_value):
        """ This function extracts from the variable_value dictionary
            - 'v_sini': projected tangential velocity of the star
            - 'radius': stellar radius in Solar units
            in addition, it takes from self.code_options
            - 'rotation_period': rotational period of the star from the stellar model

            into angular velocity Omega [rad/s] of the star
            it returns the angular velocity Omega and the stellar inclination

        Args:
            variable_value (dict): dictionary with all the parameters of the current model

        Returns:
            float: angular velocity Omega, [rad/s]
            float: inclination of the of stellar rotation axis, [degrees]
        """

        sin_is =  variable_value['v_sini']  / (variable_value['radius'] * constants.Rsun) * (self.code_options['rotation_period'] * constants.d2s)  / (2* np.pi)
        Omega = 2* np.pi / (variable_value['Prot'] * constants.d2s)

        return Omega, np.arcsin(sin_is)/np.pi*180.