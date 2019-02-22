from ..classes.common import *
from abstract_model import *


class Batman_Transit(AbstractModel):
    model_class = 'transit'
    unitary_model = True

    list_pams_common = {
        'P',  # Period, log-uniform prior
        'e',  # eccentricity, uniform prior
        'o',  # argument of pericenter (in radians)
        'R',  # planet radius (in units of stellar radii)
    }
    list_pams_dataset = {}

    default_bounds = {}
    default_spaces = {}
    default_priors = {}

    recenter_pams_dataset = {}

    transittype = 'primary'
    supersample_factor = 1
    exp_time = 0.0

    batman_params = {}
    batman_models = {}
    batman_options = {}

    ld_ncoeff = 0

    use_semimajor_axis = False
    use_inclination = False
    use_time_of_transit = False
    nthreads = 1

    def initialize_model(self, mc, **kwargs):

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

        if hasattr(kwargs, 'nthreads'):
            self.nthreads = kwargs['nthreads']

    def setup_dataset(self, dataset, **kwargs):

        self.batman_params[dataset.name_ref] = batman.TransitParams()
        self.batman_options[dataset.name_ref] = {}

        """ Initialization with random transit parameters"""
        self.batman_params[dataset.name_ref].t0 = 0.  # time of inferior conjunction
        self.batman_params[dataset.name_ref].per = 1.  # orbital period
        self.batman_params[dataset.name_ref].rp = 0.1  # planet radius (in units of stellar radii)
        self.batman_params[dataset.name_ref].a = 15.  # semi-major axis (in units of stellar radii)
        self.batman_params[dataset.name_ref].inc = 87.  # orbital inclination (in degrees)
        self.batman_params[dataset.name_ref].ecc = 0.  # eccentricity
        self.batman_params[dataset.name_ref].w = 90.  # longitude of periastron (in degrees)

        """ Setting up the limb darkening calculation"""
        try:
            self.ld_ncoeff = kwargs[dataset.name_ref]['limb_darkening_ncoeff']
            self.batman_params[dataset.name_ref].limb_dark = kwargs[dataset.name_ref]['limb_darkening_model']
        except:
            self.ld_ncoeff = kwargs['limb_darkening_ncoeff']
            self.batman_params[dataset.name_ref].limb_dark = kwargs['limb_darkening_model']

        for i_coeff in xrange(1, self.ld_ncoeff + 1):
            var = 'ld_c' + repr(i_coeff)
            self.list_pams_common.update({var: None})

        self.batman_params[dataset.name_ref].u = np.ones(self.ld_ncoeff,
                                                         dtype=np.double) * 0.1  # limb darkening coefficients

        try:
            self.batman_options[dataset.name_ref]['sample_factor'] = kwargs[dataset.name_ref]['supersample_factor']
        except:
            self.batman_options[dataset.name_ref]['sample_factor'] = kwargs['supersample_factor']

        try:
            self.batman_options[dataset.name_ref]['exp_time'] = kwargs[dataset.name_ref][
                                                                    'exposure_time'] / constants.d2s
        except:
            self.batman_options[dataset.name_ref]['exp_time'] = kwargs['exposure_time'] / constants.d2s

        self.batman_models[dataset.name_ref] = batman.TransitModel(self.batman_params[dataset.name_ref],
                                                                   dataset.x0,
                                                                   supersample_factor=
                                                                   self.batman_options[dataset.name_ref][
                                                                       'sample_factor'],
                                                                   exp_time=self.batman_options[dataset.name_ref][
                                                                       'exp_time'],
                                                                   nthreads=self.nthreads)
        self.batman_options[dataset.name_ref]['initialization_counter'] = 10000

    def compute(self, variable_value, dataset, x0_input=None):

        """
        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """

        if self.use_semimajor_axis:
            # semi-major axis (in units of stellar radii)
            self.batman_params[dataset.name_ref].a = variable_value['a']
        else:
            self.batman_params[dataset.name_ref].a = convert_rho_to_a(variable_value['P'], variable_value['rho'])

        if self.use_inclination:
            # orbital inclination (in degrees)
            self.batman_params[dataset.name_ref].inc = variable_value['i']
        else:
            self.batman_params[dataset.name_ref].inc = convert_b_to_i(variable_value['b'],
                                                                      variable_value['e'],
                                                                      variable_value['o'],
                                                                      self.batman_params[dataset.name_ref].a)

        if self.use_time_of_transit:
            self.batman_params[dataset.name_ref].t0 = variable_value['Tc'] - dataset.Tref
        else:
            self.batman_params[dataset.name_ref].t0 = kepler_exo.kepler_phase2Tc_Tref(variable_value['P'],
                                                                                      variable_value['f'],
                                                                                      variable_value['e'],
                                                                                      variable_value['o'])

        self.batman_params[dataset.name_ref].per = variable_value['P']  # orbital period
        self.batman_params[dataset.name_ref].rp = variable_value['R']  # planet radius (in units of stellar radii)
        self.batman_params[dataset.name_ref].ecc = variable_value['e']  # eccentricity
        self.batman_params[dataset.name_ref].w = variable_value['o'] * (
                    180. / np.pi)  # longitude of periastron (in degrees)

        """
        print 'a    ', self.batman_params[dataset.name_ref].a
        print 'inc  ', self.batman_params[dataset.name_ref].inc
        print 't0   ', self.batman_params[dataset.name_ref].t0
        print 'per  ', self.batman_params[dataset.name_ref].per
        print 'rp   ', self.batman_params[dataset.name_ref].rp
        print 'ecc  ', self.batman_params[dataset.name_ref].ecc
        print 'w    ', self.batman_params[dataset.name_ref].w
        print 'u    ', self.batman_params[dataset.name_ref].u
        """

        for i_coeff in xrange(1, self.ld_ncoeff + 1):
            var = 'ld_c' + repr(i_coeff)
            self.batman_params[dataset.name_ref].u[i_coeff - 1] = variable_value[var]

        """ 
        From the batman manual:
        Reinitializing the model is by far the slowest component of batman,because it calculates the optimal step size
        for the integration starting from a very small value. 
        -> However, we estimated the optimal step size from random parameters, so at some point we'll need to 
        reinitialize the model so that the correct step size is computed.
        """
        if self.batman_options[dataset.name_ref]['initialization_counter'] > 1000:
            #fac = self.batman_models[dataset.name_ref].fac
            self.batman_options[dataset.name_ref]['initialization_counter'] = 0
            self.batman_models[dataset.name_ref] = batman.TransitModel(self.batman_params[dataset.name_ref],
                                                                       dataset.x0,
                                                                       supersample_factor=
                                                                       self.batman_options[dataset.name_ref][
                                                                           'sample_factor'],
                                                                       exp_time=self.batman_options[dataset.name_ref][
                                                                           'exp_time'],
                                                                       nthreads=self.nthreads),
                                                                       #fac = fac)
        else:
            self.batman_options[dataset.name_ref]['initialization_counter'] += 1

        if x0_input is None:
            return self.batman_models[dataset.name_ref].light_curve(self.batman_params[dataset.name_ref]) - 1.

        else:
            temporary_model = batman.TransitModel(self.batman_params[dataset.name_ref],
                                                  x0_input,
                                                  supersample_factor=self.batman_options[dataset.name_ref]['sample_factor'],
                                                  exp_time=self.batman_options[dataset.name_ref]['exp_time'],
                                                  nthreads=self.nthreads)

            return temporary_model.light_curve(self.batman_params[dataset.name_ref]) - 1.
