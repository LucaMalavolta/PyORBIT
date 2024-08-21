from pyorbit.subroutines.common import np, convert_rho_to_ars, convert_b_to_i
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.keywords_definitions import *

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
        self.compute_time_inferior_conjunction = False

        """ These keywords are specific to Rossiter-McLaughlin analysis"""
        self.compute_Omega_rotation = False
        self.compute_Omega_rotation_model = -1

        """ This keywors is specific to TTV analysis analysis
        Right now implemented only in pytransit_ttv_ancillary
        """
        self.use_shared_ttvs = False

        self.limb_darkening_model = None

        self.multivariate_mass_radius = False
        self.code_options = {}

    def _prepare_planetary_parameters(self, mc, **kwargs):

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
            elif mc.common_models[self.stellar_ref].compute_density:
                self.list_pams_common.update(['mass'])
                self.list_pams_common.update(['radius'])
                self.multivariate_mass_radius = True
            else:
                """ this is the density of the star (in solar units) """
                self.list_pams_common.update(['density'])
                self.multivariate_mass_radius = False

        if mc.common_models[self.planet_ref].use_inclination:
            """ i is the orbital inclination (in degrees) """
            self.list_pams_common.update(['i'])
            self.compute_inclination = False
        else:
            """ b is the impact parameter """
            self.list_pams_common.update(['b'])

        if mc.common_models[self.planet_ref].use_time_inferior_conjunction:
            self.list_pams_common.update(['Tc'])
        else:
            self.list_pams_common.update(['mean_long'])
            self.compute_time_inferior_conjunction = True
            # mean longitude = argument of pericenter + mean anomaly at Tref

        self.use_shared_ttvs =  mc.common_models[self.planet_ref].use_shared_ttvs
        for use_shared_ttvs in keywords_shared_ttv:
            self.use_shared_ttvs = kwargs.get(use_shared_ttvs, self.use_shared_ttvs)
            if self.use_shared_ttvs:
                print('Shared transit-specific time of transits (i.e., TTVs): ', True)
                break

    def _prepare_star_parameters(self, mc, **kwargs):
        """ Additional stellar parameters
        """

        self.use_differential_rotation = kwargs.get(mc.common_models[self.stellar_ref].use_differential_rotation, False)
        for keyword in keywords_differential_rotation:
            self.use_differential_rotation = kwargs.get(keyword, self.use_differential_rotation)


        self.use_stellar_rotation_period = kwargs.get(mc.common_models[self.stellar_ref].use_stellar_rotation_period, False)
        for keyword in keywords_stellar_rotation:
            self.use_stellar_rotation_period = kwargs.get(keyword, self.use_stellar_rotation_period)


        """ check if the differential rotation should be included in the model"""
        if self.use_differential_rotation:
            #self.list_pams_common.discard('v_sini')
            self.list_pams_common.update(['alpha_rotation'])

            mc.common_models[self.stellar_ref].use_differential_rotation = True

            if self.use_stellar_rotation_period:
                mc.common_models[self.stellar_ref].use_differential_rotation = True
                mc.common_models[self.stellar_ref].use_stellar_inclination = True
                mc.common_models[self.stellar_ref].use_stellar_rotation_period = True
                mc.common_models[self.stellar_ref].use_stellar_radius = True
                mc.common_models[self.stellar_ref].use_projected_velocity = False
            else:
                mc.common_models[self.stellar_ref].use_equatorial_velocity =  True
                mc.common_models[self.stellar_ref].use_stellar_inclination =  True
                mc.common_models[self.stellar_ref].use_projected_velocity = False

            """ If stellar rotation is provided as a prior or as the outcome of the fit,
                the code will check if the derived stellar radius is consistent with the prior
            """

        if mc.common_models[self.stellar_ref].use_projected_velocity:
            self.list_pams_common.update(['v_sini'])
        else:
            self.list_pams_common.discard('v_sini')

        if mc.common_models[self.stellar_ref].use_stellar_inclination:
            self.list_pams_common.update(['i_star'])

        if mc.common_models[self.stellar_ref].use_cosine_stellar_inclination:
            self.list_pams_common.update(['cosi_star'])
            try:
                self.list_pams_common.discard('i_star')
            except:
                pass

        if mc.common_models[self.stellar_ref].use_equatorial_velocity:
            self.list_pams_common.update(['veq_star'])

        if mc.common_models[self.stellar_ref].use_stellar_rotation_period:
            self.list_pams_common.update(['rotation_period'])

        if mc.common_models[self.stellar_ref].use_stellar_radius:
            self.list_pams_common.update(['radius'])

        if mc.common_models[self.stellar_ref].use_stellar_rotation_period and mc.common_models[self.stellar_ref].use_equatorial_velocity:
            print('   ***  WARNING *** stellar rotation period and equatorial velocity  ')
            print('                    included as independent parameters ')

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


    def _prepare_limb_darkening_coefficients(self, mc, **kwargs):
        """ Setting up the limb darkening calculation"""

        self.limb_darkening_model = kwargs['limb_darkening_model']
        self.ld_vars = [0.00] * kwargs['limb_darkening_ncoeff']

        for common_model in self.common_ref:
            if mc.common_models[common_model].model_class == 'limb_darkening':
                ld_parametrization = getattr(mc.common_models[common_model], 'parametrization', 'Standard')

        if ld_parametrization=='Kipping':
            self.ldvars['ld_c1'] = 0
            self.ldvars['ld_c2'] = 1
            self.list_pams_common.update(['ld_q1'])
            self.list_pams_common.update(['ld_q2'])
        else:
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

        """ TTV MODULE
            Keep track of the boundaries of each dataset, so that the user do
            not have to write down the boundaries of each transit in case of TTV fit
            We also read the priors on Tc, if specified
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
                if kwargs.get(dict_name, False):
                    if kwargs[dict_name].get(dataset.name_ref, False):
                        tc_boundaries = kwargs[dict_name][dataset.name_ref]

        tc_priors_names = ['Tc_priors',
            'T_priors',
            'TC_priors',
        ]

        tc_priors = ['Uniform', 0.00, 0.00]
        if kwargs[dataset.name_ref].get('priors', False):
            tc_priors = np.atleast_1d(kwargs[dataset.name_ref]['priors']['Tc'])
        else:
            for dict_name in tc_priors_names:
                if kwargs[dataset.name_ref].get(dict_name, False):
                    tc_priors = kwargs[dataset.name_ref][dict_name]
                if kwargs.get(dict_name, False):
                    if kwargs[dict_name].get(dataset.name_ref, False):
                        tc_priors = kwargs[dict_name][dataset.name_ref]

        if 'Tc' in self.list_pams_dataset:
            self.bounds[dataset.name_ref]['Tc'] = tc_boundaries
            self.prior_kind[dataset.name_ref]['Tc'] = tc_priors[0]
            try:
                self.prior_pams[dataset.name_ref]['Tc'] = np.asarray(tc_priors[1:], dtype=np.double)
            except:
                self.prior_pams[dataset.name_ref]['Tc'] = np.asarray([0.00], dtype=np.double)


    """ function for internal transformation of parameters """

    def update_parameter_values(self, parameter_values, Tref, prepend=''):

        #t1_start = process_time()
        if self.multivariate_mass_radius:
            parameter_values['density'] = parameter_values['mass']/parameter_values['radius']**3

        if self.compute_semimajor_axis:
            parameter_values[prepend+'a_Rs'] = convert_rho_to_ars(parameter_values[prepend+'P'], parameter_values['density'])

        if self.compute_inclination:
            parameter_values[prepend+'i'] = convert_b_to_i(
            parameter_values[prepend+'b'], parameter_values[prepend+'e'], parameter_values[prepend+'omega'], parameter_values[prepend+'a_Rs'])

        if self.compute_time_inferior_conjunction:
            parameter_values[prepend+'Tc']= kepler_exo.kepler_phase2Tc_Tref(parameter_values[prepend+'P'],
                                               parameter_values[prepend+'mean_long'],
                                               parameter_values[prepend+'e'],
                                               parameter_values[prepend+'omega']) + Tref


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

