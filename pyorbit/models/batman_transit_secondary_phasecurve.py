
from pyorbit.subroutines.common import np, convert_rho_to_a, convert_b_to_i
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import AbstractTransit

try:
    import batman
except ImportError:
    pass


class Batman_Transit_Eclipse_PhaseCurve(AbstractModel, AbstractTransit):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)  # this calls all constructors up to AbstractModel
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            import batman
        except ImportError:
            print("ERROR: batman not installed, this will not work")
            quit()

        self.model_class = 'transit_eclipse_phasecurve'
        self.unitary_model = True

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = {
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'omega',  # argument of pericenter (in radians)
            'R_Rs',  # planet radius (in units of stellar radii)
            'phase_off',
        }
        self.list_pams_dataset = {
            'phase_amp',
            'delta_occ',
        }

        self.batman_params = None
        self.batman_transit = {}
        self.batman_eclipse = {}
        self.code_options = {
            'nthreads': 1,
            'initialization_counter': 5000
        }

    def initialize_model(self, mc, **kwargs):

        self._prepare_planetary_parameters(mc, **kwargs)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

        print('Warning: OpenMP computation on batman temporaroly turned off ')
        self.code_options['nthreads'] = 1
        self.batman_params = batman.TransitParams()

        if kwargs.get('nightside_emission', True):
            self.nightside_emission = True
        else:
            self.nightside_emission = False
            self.list_pams_dataset.discard('delta_occ')

        if kwargs.get('phase_offset', True):
            self.phase_offset = True
        else:
            self.phase_offset = False
            self.list_pams_common.discard('phase_off')

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
        self.batman_params.fp = 0.001
        self.batman_params.t_secondary = self.batman_params.per / 2.

        """ Setting up the limb darkening calculation"""

        self.batman_params.limb_dark = kwargs['limb_darkening_model']
        self.batman_params.u = np.ones(kwargs['limb_darkening_ncoeff'],
                                       dtype=np.double) * 0.1  # limb darkening coefficients

        self.code_options['initialization_counter'] = 5000

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        self._prepare_dataset_options(mc, dataset, **kwargs)

        self.batman_params.fp = 0.001
        self.batman_params.t_secondary = self.batman_params.t0 + self.batman_params.per / 2.

        self.batman_transit[dataset.name_ref] = \
            batman.TransitModel(self.batman_params,
                                dataset.x0,
                                supersample_factor=self.code_options[dataset.name_ref]['sample_factor'],
                                exp_time=self.code_options[dataset.name_ref]['exp_time'],
                                nthreads=self.code_options['nthreads'])

        self.batman_eclipse[dataset.name_ref] = \
            batman.TransitModel(self.batman_params,
                                dataset.x0,
                                supersample_factor=self.code_options[dataset.name_ref]['sample_factor'],
                                exp_time=self.code_options[dataset.name_ref]['exp_time'],
                                nthreads=self.code_options['nthreads'],
                                transittype="secondary")

    def compute(self, variable_value, dataset, x0_input=None):
        """
        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """
        #t1_start = process_time()

        pams_a, pams_i = self.retrieve_ai(variable_value)
        pams_t0 = self.retrieve_t0(variable_value, dataset.Tref)

        self.batman_params.a = pams_a
        self.batman_params.inc = pams_i
        self.batman_params.t0 = pams_t0

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

        if self.nightside_emission:
            self.batman_params.fp = variable_value['delta_occ']
        else:
            self.batman_params.fp = variable_value['phase_amp']

        if self.phase_offset:
            phase_offset = variable_value['phase_off']/180.*np.pi
        else:
            phase_offset = 0.000

        self.batman_params.t_secondary = self.batman_params.t0 + self.batman_params.per / 2.
        amplitude_sin = variable_value['phase_amp']

        """
        From the batman manual:
        Reinitializing the model is by far the slowest component of batman,because it calculates the optimal step size
        for the integration starting from a very small value.
        -> However, we estimated the optimal step size from random parameters, so at some point we'll need to
        reinitialize the model so that the correct step size is computed.
        """
        random_selector = np.random.randint(1000)

        if random_selector == 50:

            self.batman_eclipse[dataset.name_ref] = batman.TransitModel(self.batman_params,
                                            dataset.x0,
                                            supersample_factor=self.code_options[dataset.name_ref]['sample_factor'],
                                            exp_time=self.code_options[dataset.name_ref]['exp_time'],
                                            nthreads=self.code_options['nthreads'],
                                            transittype="secondary")

            self.batman_transit[dataset.name_ref] = batman.TransitModel(self.batman_params,
                                                dataset.x0,
                                                supersample_factor=self.code_options[dataset.name_ref]['sample_factor'],
                                                exp_time=self.code_options[dataset.name_ref]['exp_time'],
                                                nthreads=self.code_options['nthreads'])

        if x0_input is None:

            phase_curve =  amplitude_sin/self.batman_params.fp \
                * (np.cos(2*np.pi*(dataset.x0 - self.batman_params.t_secondary)/self.batman_params.per
                          + phase_offset)/2. + 0.5) \
                + (1 - amplitude_sin/self.batman_params.fp)

            return (self.batman_eclipse[dataset.name_ref].light_curve(self.batman_params)-1.) * phase_curve \
                + self.batman_transit[dataset.name_ref].light_curve(self.batman_params) - 1.

        else:

            phase_curve =  amplitude_sin/self.batman_params.fp \
                * (np.cos(2*np.pi*(x0_input - self.batman_params.t_secondary)/self.batman_params.per
                          + phase_offset)/2. + 0.5) \
                + (1 - amplitude_sin/self.batman_params.fp)

            batman_eclipse = batman.TransitModel(self.batman_params,
                                x0_input,
                                supersample_factor=self.code_options[dataset.name_ref]['sample_factor'],
                                exp_time=self.code_options[dataset.name_ref]['exp_time'],
                                nthreads=self.code_options['nthreads'],
                                transittype="secondary")

            batman_transit = batman.TransitModel(self.batman_params,
                                                x0_input,
                                                supersample_factor=self.code_options[dataset.name_ref]['sample_factor'],
                                                exp_time=self.code_options[dataset.name_ref]['exp_time'],
                                                nthreads=self.code_options['nthreads'])

            return (batman_eclipse.light_curve(self.batman_params)-1.) * phase_curve \
                + batman_transit.light_curve(self.batman_params) - 1.
