
from pyorbit.subroutines.common import np, convert_rho_to_a, convert_b_to_i
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.models.abstract_model import AbstractModel

try:
    import batman
except ImportError:
    pass


class Batman_Transit(AbstractModel):
    model_class = 'transit'
    unitary_model = True

    def __init__(self, *args, **kwargs):
        super(Batman_Transit, self).__init__(*args, **kwargs)

        try:
            import batman
        except ImportError:
            print("ERROR: batman package not installed, this will not work")
            quit()

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = {
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'o',  # argument of pericenter (in radians)
            'R',  # planet radius (in units of stellar radii)
        }
        self.list_pams_dataset = {}

        self.batman_ldvars = {}
        self.ld_ncoeff = 2
        self.parametrization = 'Standard'

        self.use_semimajor_axis = False
        self.use_inclination = False
        self.use_time_of_transit = False
        self.nthreads = 1

        self.transittype = 'primary'
        self.batman_params = None
        self.batman_models = {}
        self.batman_options = {
            'initialization_counter': 1000000
        }

        self.multivariate_mass_radius = False

    def initialize_model(self, mc, **kwargs):

        try:
            multivariate_vars = mc.common_models[self.stellar_ref].multivariate_vars
        except AttributeError:
            multivariate_vars = []

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

        if hasattr(kwargs, 'nthreads'):
            self.nthreads = kwargs['nthreads']

        self.batman_params = batman.TransitParams()

        """ Initialization with random transit parameters"""
        self.batman_params.t0 = 0.  # time of inferior conjunction
        self.batman_params.per = 1.  # orbital period
        # planet radius (in units of stellar radii)
        self.batman_params.rp = 0.1
        # semi-major axis (in units of stellar radii)
        self.batman_params.a = 15.
        self.batman_params.inc = 87.  # orbital inclination (in degrees)
        self.batman_params.ecc = 0.  # eccentricity
        self.batman_params.w = 90.  # longitude of periastron (in degrees)

        """ Setting up the limb darkening calculation"""

        self.batman_params.limb_dark = kwargs['limb_darkening_model']

        for i_coeff in range(1, kwargs['limb_darkening_ncoeff'] + 1):
            var = 'ld_c' + repr(i_coeff)
            self.batman_ldvars[var] = i_coeff - 1
            self.list_pams_common.update({var: None})

        self.batman_params.u = np.ones(kwargs['limb_darkening_ncoeff'],
                                       dtype=np.double) * 0.1  # limb darkening coefficients

        self.batman_options['initialization_counter'] = 5000

    def setup_dataset(self, mc, dataset, **kwargs):

        self.batman_options[dataset.name_ref] = {}

        sample_factor = 1
        exposure_time = 30.

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

        for dict_name in supersample_names:
            if kwargs[dataset.name_ref].get(dict_name, False):
                sample_factor = kwargs[dataset.name_ref][dict_name]
            elif kwargs[dataset.name_ref].get(dict_name, False):
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
            elif kwargs[dataset.name_ref].get(dict_name, False):
                exposure_time = kwargs[dict_name]

        self.batman_options[dataset.name_ref]['sample_factor'] = sample_factor
        self.batman_options[dataset.name_ref]['exp_time'] = exposure_time / constants.d2s

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


        self.batman_models[dataset.name_ref] = batman.TransitModel(self.batman_params,
                                                                   dataset.x0,
                                                                   supersample_factor=self.batman_options[dataset.name_ref][
                                                                       'sample_factor'],
                                                                   exp_time=self.batman_options[dataset.name_ref][
                                                                       'exp_time'],
                                                                   nthreads=self.nthreads)

    def compute(self, variable_value, dataset, x0_input=None):
        """
        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """
        #t1_start = process_time()

        if getattr(self, 'multivariate_mass_radius', False):
            variable_value['rho'] = variable_value['mass'] / variable_value['radius']**3

        if self.use_semimajor_axis:
            # semi-major axis (in units of stellar radii)
            self.batman_params.a = variable_value['a']
        else:
            self.batman_params.a = convert_rho_to_a(
                variable_value['P'], variable_value['rho'])

        if self.use_inclination:
            # orbital inclination (in degrees)
            self.batman_params.inc = variable_value['i']
        else:
            self.batman_params.inc = convert_b_to_i(variable_value['b'],
                                                    variable_value['e'],
                                                    variable_value['o'],
                                                    self.batman_params.a)

        if self.use_time_of_transit:
            self.batman_params.t0 = variable_value['Tc'] - dataset.Tref
        else:
            self.batman_params.t0 = kepler_exo.kepler_phase2Tc_Tref(variable_value['P'],
                                                                    variable_value['f'],
                                                                    variable_value['e'],
                                                                    variable_value['o'])

        self.batman_params.per = variable_value['P']  # orbital period
        # planet radius (in units of stellar radii)
        self.batman_params.rp = variable_value['R']
        self.batman_params.ecc = variable_value['e']  # eccentricity
        # longitude of periastron (in degrees)
        self.batman_params.w = variable_value['o'] * (180. / np.pi)


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
