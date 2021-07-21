
from pyorbit.classes.common import np, convert_rho_to_a, convert_b_to_i
import pyorbit.classes.constants as constants
import pyorbit.classes.kepler_exo as kepler_exo
from pyorbit.models.abstract_model import AbstractModel

#from time import process_time

try:
    import spiderman
except ImportError:
    pass

class Spiderman_Thermal(AbstractModel):
    model_class = 'phasecurve'
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
            'albedo' # Bond Albedo 
            'redist' # Heat redistribution
        }
        self.list_pams_dataset = {}

        self.use_semimajor_axis = False
        self.use_inclination = False
        self.use_time_of_transit = False

        self.spiderman_params = None
        self.spiderman_models = {}
        self.spiderman_options = {}

    def initialize_model(self, mc, **kwargs):

        """ check if the stellar radius and effect temperature are provided as fixed values or not """
        


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


        self.spiderman_params = spiderman.ModelParams(brightness_model='Louden',stellar_model="blackbody" )
        
        print('   Warning on Spiderman Thermal model: null limb darkening parameters for the planet')
        print('   Warning on Spiderman Thermal model: null T_int')
        self.spider_params.p_u1= 0.               # Planetary limb darkening parameter
        self.spider_params.p_u2= 0.               # Planetary limb darkening parameter
        self.T_int = 0.

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

        """ Lower and upper wavelength boundaries for the filter, in (m) """ 
        for dict_name in wavebounds_names:
            if kwargs[dataset.name_ref].get(dict_name, False):
                self.spiderman_options[dataset.name_ref]['l1'] = kwargs[dataset.name_ref][dict_name][0] / 10**10
                self.spiderman_options[dataset.name_ref]['l2'] = kwargs[dataset.name_ref][dict_name][1] / 10**10

            elif kwargs.get(dict_name, False):
                self.spiderman_options[dataset.name_ref]['l1'] = kwargs[dict_name][0] / 10**10
                self.spiderman_options[dataset.name_ref]['l2'] = kwargs[dict_name][1] / 10**10

        stellarradius_names = [
            'stellar_radius',
            'radius',
            'star_radius'
        ]
        
        """ Lower and upper wavelength boundaries for the filter, in (m) """ 
        for dict_name in wavebounds_names:
            if kwargs[dataset.name_ref].get(dict_name, False):
                self.spiderman_options[dataset.name_ref]['l1'] = kwargs[dataset.name_ref][dict_name][0] / 10**10
                self.spiderman_options[dataset.name_ref]['l2'] = kwargs[dataset.name_ref][dict_name][1] / 10**10

            elif kwargs.get(dict_name, False):
                self.spiderman_options[dataset.name_ref]['l1'] = kwargs[dict_name][0] / 10**10
                self.spiderman_options[dataset.name_ref]['l2'] = kwargs[dict_name][1] / 10**10

        self.spider_params.thermal = True
        self.spiderman_options[dataset.name_ref]['sample_factor'] = sample_factor
        self.spiderman_options[dataset.name_ref]['exp_time'] = exposure_time / constants.d2s



        
        # # OLD code snippet 
        #try:
        #    self.batman_options[dataset.name_ref]['sample_factor'] = kwargs[dataset.name_ref]['supersample_factor']
        #except:
        #    self.batman_options[dataset.name_ref]['sample_factor'] = kwargs['supersample_factor']
        #
        #try:
        #    self.batman_options[dataset.name_ref]['exp_time'] = kwargs[dataset.name_ref][
        #        'exposure_time'] / constants.d2s
        #except:
        #    self.batman_options[dataset.name_ref]['exp_time'] = kwargs['exposure_time'] / constants.d2s
        #


    def compute(self, variable_value, dataset, x0_input=None):
        """
        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """
        #t1_start = process_time()


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

        self.spider_params.a_abs= self.spider_params.a * 0.681 * constants.RsunAU        # The absolute value of the semi-major axis [AU]


        """
        print 'a    ', self.batman_params.a
        print 'inc  ', self.batman_params.inc
        print 't0   ', self.batman_params.t0
        print 'per  ', self.batman_params.per
        print 'rp   ', self.batman_params.rp
        print 'ecc  ', self.batman_params.ecc
        print 'w    ', self.batman_params.w
        print 'u    ', self.batman_params.u
        """
        for var, i_var in self.batman_ldvars.items():
            self.batman_params.u[i_var] = variable_value[var]

        """
        From the batman manual:
        Reinitializing the model is by far the slowest component of batman,because it calculates the optimal step size
        for the integration starting from a very small value. 
        -> However, we estimated the optimal step size from random parameters, so at some point we'll need to 
        reinitialize the model so that the correct step size is computed.
        """
        if self.batman_options['initialization_counter'] > 1000:
            self.batman_options['initialization_counter'] = 0
            self.batman_models[dataset.name_ref] = \
                batman.TransitModel(self.batman_params,
                                    dataset.x0,
                                    supersample_factor=self.batman_options[dataset.name_ref]['sample_factor'],
                                    exp_time=self.batman_options[dataset.name_ref][
                                        'exp_time'],
                                    nthreads=self.nthreads)

        else:
            self.batman_options['initialization_counter'] += 1




        if x0_input is None:
            ##model = self.batman_models[dataset.name_ref].light_curve(self.batman_params) - 1.
            ##t1_stop = process_time()
            ##
            ##print("Elapsed time:", t1_stop-t1_start)
            ##return model
            return self.batman_models[dataset.name_ref].light_curve(self.batman_params) - 1.

        else:
            temporary_model = batman.TransitModel(self.batman_params,
                                                  x0_input,
                                                  supersample_factor=self.batman_options[
                                                      dataset.name_ref]['sample_factor'],
                                                  exp_time=self.batman_options[dataset.name_ref]['exp_time'],
                                                  nthreads=self.nthreads)

            return temporary_model.light_curve(self.batman_params) - 1.
