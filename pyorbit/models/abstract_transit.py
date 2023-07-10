from pyorbit.subroutines.common import np, convert_rho_to_ars, convert_b_to_i
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo


class AbstractTransit(object):

    def __init__(self, *args, **kwargs):

        self.model_class = 'transit'
        self.unitary_model = True

        self.ldvars = {}
        self.ld_ncoeff = 2
        self.parametrization = 'Standard'

        """ Keywords inherited from Planet common model (with switched logical sign)"""
        self.compute_semimajor_axis = True
        self.compute_inclination = True
        self.compute_time_of_transit = False

        """ This keywors is specific to Rossiter-McLaughlin analysis"""
        self.compute_Omega_rotation = False
        self.compute_star_inclination = False
        self.compute_vsini = False

        """ This keywors is specific to TTV analysis analysis
        Right now implemented only in pytransit_ttv_ancillary
        """
        self.use_shared_ttvs = False

        """ Keywords inherited from Star_parameter common model (with switched logical sign)"""
        #self.use_stellar_rotation = False
        #self.use_stellar_inclination = False
        #self.use_equatorial_velocity = False
        #self.use_rotation_from_activity = False

        """ Some models just want fixed values for stellar raddi and temperature"""
        self.fixed_stellar_radius = False
        self.fixed_stellar_temperature = False

        self.limb_darkening_model = None

        #self.retrieve_Omega_Istar = None
        #self.retrieve_Istar = None

        self.multivariate_mass_radius = False
        self.code_options = {}

    def _prepare_planetary_parameters(self, mc, **kwargs):

        try:
            multivariate_pams = mc.common_models[self.stellar_ref].multivariate_pams
        except AttributeError:
            multivariate_pams = []

        """ Default parametrization uses the stellar density and the impact
            parameter, it is possible to switch back to scaled semi-major axis and
            inclination respectively by activating the proper flag """

        if mc.common_models[self.planet_ref].use_semimajor_axis:
            """ a is the semi-major axis (in units of stellar radii) """
            self.list_pams_common.update(['a_Rs'])
            self.compute_semimajor_axis = False
        else:
            if 'mass' in multivariate_pams and 'radius' in multivariate_pams:
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
            self.compute_inclination = False
        else:
            """ b is the impact parameter """
            self.list_pams_common.update(['b'])

        if mc.common_models[self.planet_ref].use_time_of_transit:
            self.list_pams_common.update(['Tc'])
        else:
            self.list_pams_common.update(['mean_long'])
            self.compute_time_of_transit = True
            # mean longitude = argument of pericenter + mean anomaly at Tref

        self.use_shared_ttvs = mc.common_models[self.planet_ref].use_shared_ttvs

    def _prepare_star_parameters(self, mc, **kwargs):
        """ Additional stellar parameters
            in 9.2 the possibility of fixing a parameter for a specific model
            may be suppressed
        """

        if mc.common_models[self.stellar_ref].use_stellar_inclination:
            self.list_pams_common.update(['i_star'])
        #else:
        #    self.compute_star_inclination = True

        if mc.common_models[self.stellar_ref].use_equatorial_velocity:
            self.list_pams_common.update(['veq_star'])
        #else:
        #    self.compute_equatorial_velocity = True

        if mc.common_models[self.stellar_ref].use_stellar_rotation:
            self.list_pams_common.update(['rotation_period'])
        #else:
        #    self.compute_rotation_period = True

        if mc.common_models[self.stellar_ref].use_stellar_radius:
            self.list_pams_common.update(['radius'])

        self.convective_order = kwargs.get('convective_order', mc.common_models[self.stellar_ref].convective_order)
        self.mu_step = 0.00001
        self.mu_integral = np.arange(0.,1.+self.mu_step, self.mu_step)
        if self.convective_order == 0:
            self.retrieve_convective_c0 = self._convective_c0_order0
            self.retrieve_convective_rv = self._convective_rv_order0
        elif self.convective_order == 1:
            self.retrieve_convective_c0 = self._convective_c0_order1
            self.retrieve_convective_rv = self._convective_rv_order1
            self.list_pams_common.update(['convective_c1'])
        elif self.convective_order == 2:
            self.retrieve_convective_c0 = self._convective_c0_order2
            self.retrieve_convective_rv = self._convective_rv_order2
            self.list_pams_common.update(['convective_c1'])
            self.list_pams_common.update(['convective_c2'])
        elif self.convective_order == 3:
            self.retrieve_convective_c0 = self._convective_c0_order3
            self.retrieve_convective_rv = self._convective_rv_order3
            self.list_pams_common.update(['convective_c1'])
            self.list_pams_common.update(['convective_c2'])
            self.list_pams_common.update(['convective_c3'])
        else:
            print(' Maximum order for polynomial model for convective RV shift is 3')
            quit()


        stellar_radius_names = [
            'stellar_radius',
            'radius',
            'star_radius'
        ]

        for dict_name in stellar_radius_names:
            if kwargs.get(dict_name, False):
                self.code_options['radius'] = kwargs[dict_name]
                self.list_pams_common.discard('radius')
                self.fixed_stellar_radius = True

        effective_temperature_names = [
            'teff',
            'temperature',
            'eff_temperature'
            'temperature_eff'
        ]
        for dict_name in effective_temperature_names:
            if kwargs.get(dict_name, False):
                self.code_options['temperature'] = kwargs[dict_name]
                self.list_pams_common.discard('temperature')
                self.fixed_stellar_temperature = True

        """ Check if the stellar rotation period is given as a starting value
            If so, the angular rotation of the star and the stellar inclination
            are computed through the rotation period
        """

        self.code_options['rotation_keyword'] = 'rotation_period'
        rotation_period_names =[
            'rotation_period',
            'rotational_period',
            'stellar_period',
            'stellar_rotation',
            'Prot',
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
        for dict_name in rotation_period_names:
            if kwargs.get(dict_name, False):
                self.code_options['rotation_period'] = kwargs[dict_name]
                self.list_pams_common.discard('rotation_period')
                self.fixed_stellar_rotation = True

    def _prepare_limb_darkening_coefficients(self, mc, **kwargs):
        """ Setting up the limb darkening calculation"""

        self.limb_darkening_model = kwargs['limb_darkening_model']
        self.ld_vars = [0.00] * kwargs['limb_darkening_ncoeff']
        for i_coeff in range(1, kwargs['limb_darkening_ncoeff'] + 1):
            par = 'ld_c' + repr(i_coeff)
            self.ldvars[par] = i_coeff - 1
            self.list_pams_common.update([par])

        if self.limb_darkening_model == 'uniform':
            self.compute_limb_darkening = self._limb_darkening_uniform
        elif self.limb_darkening_model == 'linear':
            self.compute_limb_darkening = self._limb_darkening_linear
        elif self.limb_darkening_model == 'quadratic':
            self.compute_limb_darkening = self._limb_darkening_quadratic
        elif self.limb_darkening_model == 'square-root':
            self.compute_limb_darkening = self._limb_darkening_squareroot
        elif self.limb_darkening_model == 'logarithmic':
            self.compute_limb_darkening = self._limb_darkening_logarithmic
        elif self.limb_darkening_model == 'exponential':
            self.compute_limb_darkening = self._limb_darkening_exponential
        else:
            print('ERROR: Selected limb darkening law not implemented')
            quit()

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


    """ function for internal transformation of parameters """

    def update_parameter_values(self, parameter_values, Tref):

        #t1_start = process_time()
        if self.multivariate_mass_radius:
            parameter_values['density'] = parameter_values['mass']/parameter_values['radius']**3

        if self.compute_semimajor_axis:
            parameter_values['a_Rs'] = convert_rho_to_ars(parameter_values['P'], parameter_values['density'])

        if self.compute_inclination:
            parameter_values['i'] = convert_b_to_i(
            parameter_values['b'], parameter_values['e'], parameter_values['omega'], parameter_values['a_Rs'])

        if self.compute_time_of_transit:
            parameter_values['Tc']= kepler_exo.kepler_phase2Tc_Tref(parameter_values['P'],
                                               parameter_values['mean_long'],
                                               parameter_values['e'],
                                               parameter_values['omega']) + Tref

        if self.fixed_stellar_radius:
            parameter_values['radius'] = self.code_options['radius']

        if self.fixed_stellar_temperature:
            parameter_values['temperature'] = self.code_options['temperature']



    def _limb_darkening_coefficients(self, parameter_values):
        ld_par = np.zeros(2)
        for par, i_par in self.ldvars.items():
            ld_par[i_par] = parameter_values[par]
        return ld_par

    def _limb_darkening_uniform(self, ld_par, mu):
        return  1

    def _limb_darkening_linear(self, ld_par, mu):
        return  1 - ld_par[0]*(1. - mu)

    def _limb_darkening_quadratic(self, ld_par, mu):
        return  1 - ld_par[0]*(1. - mu) - ld_par[1]*(1. - mu)**2

    def _limb_darkening_squareroot(self, ld_par, mu):
        return  1 - ld_par[0]*(1. - mu) - ld_par[1]*(1. - np.sqrt(mu))

    def _limb_darkening_logarithmic(self, ld_par, mu):
        return  1 - ld_par[0]*(1. - mu) - ld_par[1]*mu*np.log(mu)

    def _limb_darkening_exponential(self, ld_par, mu):
        return  1 - ld_par[0]*(1. - mu) - ld_par[1]/(1. - np.exp(mu))




    def _convective_c0_order0(self, ld_par, parameter_values):
        return 0

    def _convective_rv_order0(self, mean_mu, parameter_values):
        return 0

    def _convective_c0_order1(self, ld_par, parameter_values):

        I_integral = self.compute_limb_darkening(ld_par, self.mu_integral)

        int_1=parameter_values['convective_c1']*np.sum(I_integral*(self.mu_integral**2.)*self.mu_step)
        dnm=np.sum(I_integral*self.mu_integral*self.mu_step)
        return -(int_1)/dnm

    def _convective_rv_order1(self, mean_mu, parameter_values):

        return parameter_values['convective_c1']*mean_mu

    def _convective_c0_order2(self, ld_par, parameter_values):

        I_integral = self.compute_limb_darkening(ld_par, self.mu_integral)

        int_1=parameter_values['convective_c1']*np.sum(I_integral*(self.mu_integral**2.)*self.mu_step)
        int_2=parameter_values['convective_c2']*np.sum(I_integral*(self.mu_integral**3.)*self.mu_step)
        dnm=np.sum(I_integral*self.mu_integral*self.mu_step)
        return -(int_1+int_2)/dnm

    def _convective_rv_order2(self, mean_mu, parameter_values):
        return parameter_values['convective_c1']*mean_mu+(parameter_values['convective_c2']*mean_mu**2)

    def _convective_c0_order3(self, ld_par, parameter_values):

        I_integral = self.compute_limb_darkening(ld_par, self.mu_integral)

        int_1=parameter_values['convective_c1']*np.sum(I_integral*(self.mu_integral**2.)*self.mu_step)
        int_2=parameter_values['convective_c2']*np.sum(I_integral*(self.mu_integral**3.)*self.mu_step)
        int_3=parameter_values['convective_c3']*np.sum(I_integral*(self.mu_integral**4.)*self.mu_step)
        dnm=np.sum(I_integral*self.mu_integral*self.mu_step)
        return -(int_1+int_2+int_3)/dnm

    def _convective_rv_order3(self, mean_mu, parameter_values):
        return ((parameter_values['convective_c1']*mean_mu)
            + (parameter_values['convective_c2']*mean_mu**2)
            + (parameter_values['convective_c3']*mean_mu**3))

