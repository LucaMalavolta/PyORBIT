from pyorbit.subroutines.common import np, convert_rho_to_a, convert_b_to_i
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import AbstractTransit


class Batman_Transit_With_TTV(AbstractModel, AbstractTransit):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)  # this calls all constructors up to AbstractModel
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            import batman
        except ImportError:
            print("ERROR: batman not installed, this will not work")
            quit()

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = {
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'o',  # argument of pericenter (in radians)
            'R',  # planet radius (in units of stellar radii)
        }

        self.list_pams_dataset = {
            'deltaT' # time offset for the central time of transit
        }
        self.default_bounds = {
            'deltaT': [-1., 1.]
        }
        self.default_spaces = {
            'deltaT': 'Linear'
        }
        self.default_priors = {
            'deltaT': ['Uniform', []]
        }
        self.default_fixed = {
            'deltaT': 0.00000
        }

        self.batman_params = None
        self.batman_models = {}
        self.code_options = {
            'nthreads': 1,
            'initialization_counter': 5000
        }

    def initialize_model(self, mc, **kwargs):

        self._prepare_planetary_parameters(mc, **kwargs)
        self._prepare_limnb_darkening_coefficients(mc, **kwargs)

        if hasattr(kwargs, 'nthreads'):
            self.code_options['nthreads'] = kwargs['nthreads']

        self.batman_params = batman.TransitParams()

        """ Initialization with random transit parameters"""
        self.batman_params.t0 = 0.  # time of inferior conjunction
        self.batman_params.per = 1.  # orbital period
        self.batman_params.rp = 0.1  # planet radius (in units of stellar radii)
        self.batman_params.a = 15.  # semi-major axis (in units of stellar radii)
        self.batman_params.inc = 87.  # orbital inclination (in degrees)
        self.batman_params.ecc = 0.  # eccentricity
        self.batman_params.w = 90.  # longitude of periastron (in degrees)

        """ Setting up the limb darkening calculation"""

        self.batman_params.limb_dark = kwargs['limb_darkening_model']
        self.batman_params.u = np.ones(kwargs['limb_darkening_ncoeff'],
                                       dtype=np.double) * 0.1  # limb darkening coefficients

        self.code_options['initialization_counter'] = 5000

    def setup_dataset(self, mc, dataset, **kwargs):

        self._prepare_dataset_options(mc, dataset, **kwargs)
        self.batman_models[dataset.name_ref] = batman.TransitModel(self.batman_params,
                                                                   dataset.x0,
                                                                   supersample_factor=
                                                                   self.code_options[dataset.name_ref][
                                                                       'sample_factor'],
                                                                   exp_time=self.code_options[dataset.name_ref][
                                                                       'exp_time'],
                                                                   nthreads=self.code_options['nthreads'])

    def compute(self, variable_value, dataset, x0_input=None):

        """
        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """

        self.batman_params.a, self.batman_params.inc = self.retrieve_ai(variable_value)
        self.batman_params.t0 = self.retrieve_t0(variable_value, dataset.Tref) + variable_value['deltaT']

        self.batman_params.per = variable_value['P']  # orbital period
        self.batman_params.rp = variable_value['R']  # planet radius (in units of stellar radii)
        self.batman_params.ecc = variable_value['e']  # eccentricity
        self.batman_params.w = variable_value['o'] * (180. / np.pi)  # longitude of periastron (in degrees)

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
        for var, i_var in self.ldvars.items():
            self.batman_params.u[i_var] = variable_value[var]

        """
        From the batman manual:
        Reinitializing the model is by far the slowest component of batman,because it calculates the optimal step size
        for the integration starting from a very small value. 
        -> However, we estimated the optimal step size from random parameters, so at some point we'll need to 
        reinitialize the model so that the correct step size is computed.
        """
        if self.code_options['initialization_counter'] > 1000:
            self.code_options['initialization_counter'] = 0
            self.batman_models[dataset.name_ref] = batman.TransitModel(self.batman_params,
                                                                       dataset.x0,
                                                                       supersample_factor=
                                                                       self.code_options[dataset.name_ref][
                                                                           'sample_factor'],
                                                                       exp_time=self.code_options[dataset.name_ref][
                                                                           'exp_time'],
                                                                       nthreads=self.code_options['nthreads'])
        else:
            self.code_options['initialization_counter'] += 1

        if x0_input is None:
            return self.batman_models[dataset.name_ref].light_curve(self.batman_params) - 1.

        else:
            temporary_model = batman.TransitModel(self.batman_params,
                                                  x0_input,
                                                  supersample_factor=self.code_options[dataset.name_ref]['sample_factor'],
                                                  exp_time=self.code_options[dataset.name_ref]['exp_time'],
                                                  nthreads=self.code_options['nthreads'])

            return temporary_model.light_curve(self.batman_params) - 1.
