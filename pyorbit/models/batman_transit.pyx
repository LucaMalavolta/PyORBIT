from ..classes.common import *
from abstract_model import *

class Batman_Transit(AbstractModel):

    model_class = 'transit'
    multiplicative_model = True

    list_pams_common = {
        'P',  # Period, log-uniform prior
        'f',  # mean longitude = argument of pericenter + mean anomaly at Tref
        'e',  # eccentricity, uniform prior - to be fixed
        'o',  # argument of pericenter (in radians)
        'i',  # orbital inclination (in degrees)
        'R',  # planet radius (in units of stellar radii)
        'a'  # semi-major axis (in units of stellar radii)
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

    #def initialize_model(self, mc, **kwargs):

    def setup_dataset(self, dataset, **kwargs):

        self.batman_params[dataset.name_ref] = batman.TransitParams()
        self.batman_options[dataset.name_ref] = {}

        """ Initialization with random transit parameters"""
        self.batman_params[dataset.name_ref].t0 = 0. #time of inferior conjunction
        self.batman_params[dataset.name_ref].per = 1. #orbital period
        self.batman_params[dataset.name_ref].rp = 0.1 #planet radius (in units of stellar radii)
        self.batman_params[dataset.name_ref].a = 15. #semi-major axis (in units of stellar radii)
        self.batman_params[dataset.name_ref].inc = 87. #orbital inclination (in degrees)
        self.batman_params[dataset.name_ref].ecc = 0. #eccentricity
        self.batman_params[dataset.name_ref].w = 90. #longitude of periastron (in degrees)


        self.batman_params[dataset.name_ref].u = [0.1, 0.3] #limb darkening coefficients



        try:
            self.batman_params[dataset.name_ref].limb_dark = kwargs[dataset.name_ref]['limb_darkening']
        except:
            self.batman_params[dataset.name_ref].limb_dark = kwargs['limb_darkening']

        try:
            self.batman_options[dataset.name_ref]['sample_factor'] = kwargs[dataset.name_ref]['supersample_factor']
        except:
            self.batman_options[dataset.name_ref]['sample_factor'] = kwargs['supersample_factor']

        try:
            self.batman_options[dataset.name_ref]['exp_time'] = kwargs[dataset.name_ref]['exposure_time']
        except:
            self.batman_options[dataset.name_ref]['exp_time'] = kwargs['exposure_time']

        self.batman_models[dataset.name_ref] = batman.TransitModel(self.batman_params[dataset.name_ref],
                                                                   dataset.x0,
                                                                   supersample_factor=self.batman_options[dataset.name_ref]['sample_factor'],
                                                                   exp_time=self.batman_options[dataset.name_ref]['exp_time'])
        self.batman_options[dataset.name_ref]['initialization_counter'] = 10000

    def compute(self, variable_value, dataset, x0_input=None):

        """

        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """

        self.batman_params[dataset.name_ref].t0 = kepler_exo.kepler_Tcent_T0P(variable_value['P'],
                                                     variable_value['f'],
                                                     variable_value['e'],
                                                     variable_value['o'])

        self.batman_params[dataset.name_ref].per = variable_value['P'] #orbital period
        self.batman_params[dataset.name_ref].rp = variable_value['R']  #planet radius (in units of stellar radii)
        self.batman_params[dataset.name_ref].a = variable_value['a']   #semi-major axis (in units of stellar radii)
        self.batman_params[dataset.name_ref].inc = variable_value['i'] #orbital inclination (in degrees)
        self.batman_params[dataset.name_ref].ecc = variable_value['e'] #eccentricity
        self.batman_params[dataset.name_ref].w = variable_value['o'] * (180./np.pi) #longitude of periastron (in degrees)
        self.batman_params[dataset.name_ref].u = [0.1, 0.3] #limb darkening coefficients


        """ 
        From the batman manual:
        Reinitializing the model is by far the slowest component of batman,because it calculates the optimal step size
        for the integration starting from a very small value. 
        -> However, we estimated the optimal step size from random parameters, so at some point we'll need to 
        reinitialize the model so that the correct step size is computed.
        """
        if self.batman_options[dataset.name_ref]['initialization_counter'] > 100:
            self.batman_options[dataset.name_ref]['initialization_counter'] = 0
            self.batman_models[dataset.name_ref] = batman.TransitModel(self.batman_params[dataset.name_ref],
                                                                       dataset.x0,
                                                                       supersample_factor=
                                                                       self.batman_options[dataset.name_ref][
                                                                           'sample_factor'],
                                                                       exp_time=self.batman_options[dataset.name_ref][
                                                                           'exp_time'])
        else:
            self.batman_options[dataset.name_ref]['initialization_counter'] += 1

        if x0_input is None:
            return 1. - self.batman_models[dataset.name_ref].light_curve(self.batman_params[dataset.name_ref])
        else:
            temporary_model = batman.TransitModel(self.batman_params[dataset.name_ref],
                                                                       x0_input,
                                                                       supersample_factor=
                                                                       self.batman_options[dataset.name_ref][
                                                                           'sample_factor'],
                                                                       exp_time=self.batman_options[dataset.name_ref][
                                                                           'exp_time'])
            return 1. - temporary_model.light_curve(self.batman_params[dataset.name_ref])
