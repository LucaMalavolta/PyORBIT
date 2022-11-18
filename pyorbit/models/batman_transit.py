
from pyorbit.subroutines.common import np, convert_rho_to_a, convert_b_to_i
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import AbstractTransit

try:
    import batman
except ImportError:
    pass


class Batman_Transit(AbstractModel, AbstractTransit):

    def __init__(self, *args, **kwargs):
        # this calls all constructors up to AbstractModel
        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            import batman
        except ImportError:
            print("ERROR: batman package not installed, this will not work")
            quit()

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = {
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'omega',  # argument of pericenter (in radians)
            'R_Rs',  # planet radius (in units of stellar radii)
        }

        """ Model-specifc variables, not declared in the abstract class """
        self.batman_params = None
        self.batman_models = {}
        self.code_options = {}
        self.transit_time_boundaries = {}

        self.dataset_x0 ={}
        self.dataset_Tref = {}
        #self.copied_dataset = True


    def initialize_model(self, mc, **kwargs):

        self._prepare_planetary_parameters(mc, **kwargs)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

        self.code_options['nthreads'] = kwargs.get('nthreads', 1)

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
        # limb darkening coefficients
        self.batman_params.u = np.ones(
            kwargs['limb_darkening_ncoeff'], dtype=np.double) * 0.1

        self.code_options['initialization_counter'] = 5000

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        self._prepare_dataset_options(mc, dataset, **kwargs)
        #self.batman_models[dataset.name_ref] = \
        #    batman.TransitModel(self.batman_params,
        #                        dataset.x0,
        #                        supersample_factor=self.code_options[dataset.name_ref]['sample_factor'],
        #                        exp_time=self.code_options[dataset.name_ref]['exp_time'],
        #                        nthreads=self.code_options['nthreads'])
        self.batman_models[dataset.name_ref] = \
            batman.TransitModel(self.batman_params,
                                dataset.x0,
                                supersample_factor=self.code_options[dataset.name_ref]['sample_factor'],
                                exp_time=self.code_options[dataset.name_ref]['exp_time'],
                                nthreads=self.code_options['nthreads'])
        """ Keep track of the boundaries of each dataset, so that the user do
        not have to write down the boundaries of each transit in case of TTV fit
        """
        self.transit_time_boundaries[dataset.name_ref] = \
            [np.amin(dataset.x), np.amax(dataset.x)]

    def compute(self, variable_value, dataset, x0_input=None):
        """
        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """
        #t1_start = process_time()

        self.batman_params.a, self.batman_params.inc = self.retrieve_ai(
            variable_value)
        self.batman_params.t0 = self.retrieve_t0(variable_value, dataset.Tref)

        self.batman_params.per = variable_value['P']  # orbital period
        # planet radius (in units of stellar radii)
        self.batman_params.rp = variable_value['R_Rs']
        self.batman_params.ecc = variable_value['e']  # eccentricity
        # longitude of periastron (in degrees)
        self.batman_params.w = variable_value['omega']

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
        Reinitializing the model is by far the slowest component of batman,
        because it calculates the optimal step size for the integration starting
        from a very small value.
        -> However, we estimated the optimal step size from random parameters,
        so at some point we'll need to reinitialize the model so that the
        correct step size is computed.
        """

        random_selector = np.random.randint(1000)
        if random_selector == 50:
            self.batman_models[dataset.name_ref] = \
                batman.TransitModel(self.batman_params,
                                    dataset.x0,
                                    supersample_factor=self.code_options[dataset.name_ref]['sample_factor'],
                                    exp_time=self.code_options[dataset.name_ref]['exp_time'],
                                    nthreads=self.code_options['nthreads'])

        if x0_input is None:
            return self.batman_models[dataset.name_ref].light_curve(self.batman_params) - 1.

        else:
            temporary_model = batman.TransitModel(self.batman_params,
                                                  x0_input,
                                                  supersample_factor=self.code_options[
                                                      dataset.name_ref]['sample_factor'],
                                                  exp_time=self.code_options[dataset.name_ref]['exp_time'],
                                                  nthreads=self.code_options['nthreads'])

            return temporary_model.light_curve(self.batman_params) - 1.
