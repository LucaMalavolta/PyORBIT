from pyorbit.subroutines.common import np, OrderedSet
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import *
from scipy.ndimage import gaussian_filter1d


class RossiterMcLaughlin_Exomoon(AbstractModel, AbstractTransit):
    model_class = 'rossiter_mclaughlin_exomoon'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)  # this calls all constructors up to AbstractModel
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            from lmfit.models import GaussianModel
        except (ModuleNotFoundError,ImportError):
            print("ERROR: lmfit not installed, this will not work")
            quit()

        import warnings
        from scipy.optimize import OptimizeWarning
        warnings.filterwarnings("ignore", category=OptimizeWarning)

        self.unitary_model = False

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = OrderedSet([
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'omega',  # argument of pericenter (in radians)
            'lambda', # Sky-projected angle between stellar rotation axis and normal of orbit plane [deg]
            'R_Rs',  # planet radius (in units of stellar radii)
            'em_P',
            'em_e',  # eccentricity, uniform prior
            'em_omega',  # argument of pericenter (in radians)
            'em_Tc',  # time of inferior conjunction (in days)
            'natural_contrast',  # relative depth of the CCF
            'natural_broadening',  # natural broadening of the CCF (in km/s)
            'instrumental_broadening',  # instrumental broadening of the CCF (km/s)
            'quadratic_ld_c1',  # quadratic limb darkening coefficient c1
            'quadratic_ld_c2',  # quadratic limb darkening coefficient c2
        ])

        self.star_grid = {}   # write an empty dictionary
        self.model_class = 'rossiter_mclaughlin_exomoon'


    def print_info(self):
        print("*** model {0:s} parameters:".format(self.model_name))
        print("    instrument: ", self.spectrograph_ref)
        print("    Note: assumption of quadratic limb darkening for RML computation")
        for key, value in self.ccf_variables.items():
            print(f"        {key}: {value}")
        print()

    def initialize_model(self, mc, **kwargs):

        self._prepare_planet_parameters(mc, **kwargs)
        self._prepare_star_parameters(mc, **kwargs)

        for common_ref in self.common_ref:
            if getattr(mc.common_models[common_ref], 'model_class', None) == 'spectrograph':
                self.spectrograph_ref = common_ref
                break

        self.use_instrument_specific_line = not mc.common_models[self.spectrograph_ref].use_stellar_lines
        if self.use_instrument_specific_line:
            self.list_pams_common.discard('natural_contrast')
            self.list_pams_common.discard('natural_broadening')
            self.list_pams_common.update(['line_broadening'])
            self.list_pams_common.update(['line_contrast'])

        # HARPS-N instrumental braodening: FWHM of 3.1 pixels
        # At 530 nm for example, in the center of the CCD,
        # the scale is 0.001415 nm/pixel, and the spectral
        # resolution taking into account the measured spotsize
        # is computed to about R = 124’000.
        # fwhm = 2.355 sigma
        # fwhm = 0.01415 * 3.1 * 299792.458 / 5300 = 2.48 km/s
        # instrumental_broadening =  fwhm / 2.355 = 1.108 km/s

        """ Values from spectrograph common object"""
        self.ccf_variables = {
            'rv_min': mc.common_models[ self.spectrograph_ref].rv_min,
            'rv_max': mc.common_models[ self.spectrograph_ref].rv_max,
            'rv_step': mc.common_models[ self.spectrograph_ref].rv_step,
        }


    def initialize_model_dataset(self, mc, dataset, **kwargs):
    
        self._prepare_dataset_options(mc, dataset, **kwargs)
        #self.batman_models[dataset.name_ref] = \
        #    batman.TransitModel(self.batman_params,
        #                        dataset.x0,
        #                        supersample_factor=self.code_options[dataset.name_ref]['sample_factor'],
        #                        exp_time=self.code_options[dataset.name_ref]['exp_time'],
        #                        nthreads=self.code_options['nthreads'])

    def precompute(self, parameter_values, dataset):

        for key, key_val in parameter_values.items():
            if np.isnan(key_val):
                pass 

        if self.precomputed == True:
            pass


    def compute(self, parameter_values, dataset, x0_input=None):

        """
        :param parameter_values:
        :param dataset:
        :param x0_input:
        :return:
        """
        self.update_parameter_values(parameter_values)

        for key, key_val in parameter_values.items():
            if np.isnan(key_val):
                return 0.

        if x0_input is None:
            bjd = dataset.x0
            rv_rml =  dataset.x0 + np.random.normal(loc=0.0, scale=1, size=len(dataset.x0))

        else:
            bjd = x0_input
            rv_rml =  x0_input + np.random.normal(loc=0.0, scale=1, size=len(dataset.x0))


        print(' rml exomoon model: ', self.model_name, ' for dataset: ', dataset.name_ref)
        print(parameter_values)

        quit()
        return rv_rml
