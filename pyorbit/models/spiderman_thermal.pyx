
from pyorbit.classes.common import np, convert_rho_to_a, convert_b_to_i
import pyorbit.classes.constants as constants
import pyorbit.classes.kepler_exo as kepler_exo
from pyorbit.models.abstract_model import AbstractModel

#from time import process_time

#try:
#    import spiderman
#except ImportError:
#    pass


class Spiderman_Thermal(AbstractModel):
    model_class = 'eclipse_phasecurve'
    unitary_model = True

    default_bounds = {}
    default_spaces = {}
    default_priors = {}

    recenter_pams_dataset = {}

    def __init__(self, *args, **kwargs):

        super(Spiderman_Thermal, self).__init__(*args, **kwargs)

        try:
            import spiderman
        except ImportError:
            print("ERROR: spiderman not installed, this will not work")
            quit()

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = {
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'o',  # argument of pericenter (in radians)
            'R',  # planet radius (in units of stellar radii)
            'albedo',  # Bond Albedo
            'redist',  # Heat redistribution
            #'insol'
        }
        self.list_pams_dataset = {}

        self.use_semimajor_axis = False
        self.use_inclination = False
        self.use_time_of_transit = False
        self.use_stellar_radius = True
        self.use_stellar_temperature = True

        self.spiderman_params = None
        self.spiderman_models = {}
        self.spiderman_options = {}

    def initialize_model(self, mc, **kwargs):

        """ check if the stellar radius and effect temperature are provided as fixed values or not """
        stellarradius_names = [
            'stellar_radius',
            'radius',
            'star_radius'
        ]
        for dict_name in stellarradius_names:
            if kwargs.get(dict_name, False):
                self.spiderman_options['radius'] = kwargs[dict_name]
                self.use_stellar_radius = False
        if self.use_stellar_radius:
            self.list_pams_common.update({'radius': None})

        effectivetemperature_names = [
            'teff',
            'temperature',
            'eff_temperature'
        ]
        for dict_name in effectivetemperature_names:
            if kwargs.get(dict_name, False):
                self.spiderman_options['temperature'] = kwargs[dict_name]
                self.use_stellar_temperature = False
        if self.use_stellar_temperature:
            self.list_pams_common.update({'temperature': None})

        """ Default parametrization uses the stellar density and the impact
            parameter, it is possible to switch back to scaled semi-major axis and
            inclination respectively by activating the proper flag """

        if mc.common_models[self.planet_ref].use_semimajor_axis:
            """ a is the semi-major axis (in units of stellar radii) """
            self.list_pams_common.update({'a': None})
            self.use_semimajor_axis = True
        else:
            """ rho is the density of the star (in solar units) """
            self.list_pams_common.update({'rho': None})

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

        self.spiderman_params = spiderman.ModelParams(
            brightness_model='Louden', stellar_model="blackbody")

        print('   Warning on Spiderman Thermal model: null limb darkening parameters for the planet')
        print('   Warning on Spiderman Thermal model: null T_int')
        self.spiderman_params.p_u1 = 0.               # Planetary limb darkening parameter
        self.spiderman_params.p_u2 = 0.               # Planetary limb darkening parameter
        self.spiderman_params.T_int = 0.

    def setup_dataset(self, mc, dataset, **kwargs):

        self.spiderman_options[dataset.name_ref] = {}

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

        sample_factor = 0
        exposure_time = 30.

        for dict_name in supersample_names:
            if kwargs[dataset.name_ref].get(dict_name, False):
                sample_factor = kwargs[dataset.name_ref][dict_name]
            elif kwargs.get(dict_name, False):
                sample_factor = kwargs[dict_name]

        exptime_names = [
            'exposure_time',
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

        wavebounds_names = [
            'wavelength_range',
            'wavelength_boundaries',
        ]

        """ Lower and upper wavelength boundaries for the filter, in (nm) """
        for dict_name in wavebounds_names:
            if kwargs[dataset.name_ref].get(dict_name, False):
                self.spiderman_options[dataset.name_ref]['l1'] = kwargs[dataset.name_ref][dict_name][0] / 10**9
                self.spiderman_options[dataset.name_ref]['l2'] = kwargs[dataset.name_ref][dict_name][1] / 10**9
            elif kwargs.get(dict_name, False):
                self.spiderman_options[dataset.name_ref]['l1'] = kwargs[dict_name][0] / 10**9
                self.spiderman_options[dataset.name_ref]['l2'] = kwargs[dict_name][1] / 10**9

        self.spiderman_params.thermal = True
        self.spiderman_options[dataset.name_ref]['sample_factor'] = sample_factor
        self.spiderman_options[dataset.name_ref]['exp_time'] = exposure_time / constants.d2s

        # # OLD code snippet
        # try:
        #    self.batman_options[dataset.name_ref]['sample_factor'] = kwargs[dataset.name_ref]['supersample_factor']
        # except:
        #    self.batman_options[dataset.name_ref]['sample_factor'] = kwargs['supersample_factor']
        #
        # try:
        #    self.batman_options[dataset.name_ref]['exp_time'] = kwargs[dataset.name_ref][
        #        'exposure_time'] / constants.d2s
        # except:
        #    self.batman_options[dataset.name_ref]['exp_time'] = kwargs['exposure_time'] / constants.d2s
        #

    def compute(self, variable_value, dataset, x0_input=None):
        """[summary]

        Args:
            variable_value ([type]): [description]
            dataset ([type]): [description]
            x0_input ([type], optional): [description]. Defaults to None.
        """

        """ retrieve stellar temperature and radius either from the options or
            from the priors
        """
        if self.use_stellar_radius:
            self.spiderman_options['radius'] = variable_value['radius']
        else:
            stellar_radius = self.spiderman_options['radius']

        if self.use_stellar_temperature:
            self.spiderman_params.T_s = variable_value['temperature']
        else:
            self.spiderman_params.T_s = self.spiderman_options['temperature']

        if self.use_semimajor_axis:
            # semi-major axis (in units of stellar radii)
            self.spiderman_params.a = variable_value['a']
        else:
            self.spiderman_params.a = convert_rho_to_a(
                variable_value['P'], variable_value['rho'])

        if self.use_inclination:
            # orbital inclination (in degrees)
            self.spiderman_params.inc = variable_value['i']
        else:
            self.spiderman_params.inc = convert_b_to_i(variable_value['b'],
                                                       variable_value['e'],
                                                       variable_value['o'],
                                                       self.spiderman_params.a)

        if self.use_time_of_transit:
            self.spiderman_params.t0 = variable_value['Tc'] - dataset.Tref
        else:
            self.spiderman_params.t0 = kepler_exo.kepler_phase2Tc_Tref(variable_value['P'],
                                                                       variable_value['f'],
                                                                       variable_value['e'],
                                                                       variable_value['o'])

        self.spiderman_params.per = variable_value['P']  # orbital period
        # planet radius (in units of stellar radii)
        self.spiderman_params.rp = variable_value['R']
        self.spiderman_params.ecc = variable_value['e']  # eccentricity
        # longitude of periastron (in degrees)
        self.spiderman_params.w = variable_value['o'] * (180. / np.pi)

        self.spiderman_params.l1 = self.spiderman_options[dataset.name_ref]['l1']
        self.spiderman_params.l2 = self.spiderman_options[dataset.name_ref]['l2']

        # The absolute value of the semi-major axis [AU]
        self.spiderman_params.a_abs = (self.spiderman_params.a
                                    * self.spiderman_options['radius']
                                    * constants.RsunAU)


        self.spiderman_params.insol = self.spiderman_options['radius']**2 \
            * (self.spiderman_params.T_s/5777.0)**4 \
                * self.spiderman_params.a_abs ** 2 \
                    * 1367 

        self.spiderman_params.albedo = variable_value['albedo']
        self.spiderman_params.redist = variable_value['redist']
        #self.spiderman_params.insol = variable_value['insol']

        if x0_input is None:
            if self.spiderman_params.a < 1.0: return dataset.x0 * 0.

            if self.spiderman_options[dataset.name_ref]['sample_factor'] > 1:

                time_array = np.linspace(-self.spiderman_options[dataset.name_ref]['exp_time']/2.,
                                        self.spiderman_options[dataset.name_ref]['exp_time']/2.,
                                        self.spiderman_options[dataset.name_ref]['sample_factor'])

                lc_out = dataset.x0 * 0.
                for it, tt in enumerate(dataset.x0):
                    lc_out[it] = np.average(self.spiderman_params.lightcurve(tt + time_array))
                return lc_out - 1.
            else:
                return self.spiderman_params.lightcurve(dataset.x0) - 1.


        else:
            if self.spiderman_params.a < 1.0: return x0_input*0.

            if self.spiderman_options[dataset.name_ref]['sample_factor'] > 1:

                time_array = np.linspace(-self.spiderman_options[dataset.name_ref]['exp_time']/2.,
                                        self.spiderman_options[dataset.name_ref]['exp_time']/2.,
                                        self.spiderman_options[dataset.name_ref]['sample_factor'])

                lc_out = x0_input * 0.
                for it, tt in enumerate(x0_input):
                    lc_out[it] = np.average(self.spiderman_params.lightcurve(tt + time_array))
                return lc_out - 1.
            else:
                return self.spiderman_params.lightcurve(x0_input) - 1.
