from pyorbit.subroutines.common import np, OrderedSet, convert_ars_to_a, convert_RTaAU_to_insol
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import AbstractTransit

class Spiderman_Thermal(AbstractModel, AbstractTransit):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)  # this calls all constructors up to AbstractModel
        super(AbstractModel, self).__init__(*args, **kwargs)

        model_class = 'eclipse_phasecurve'
        unitary_model = True

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = OrderedSet([
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'omega',  # argument of pericenter (in radians)
            'R_Rs',  # planet radius (in units of stellar radii)
            'albedo',  # Bond Albedo
            'redist',  # Heat redistribution
            #'insol'
        ])

        self.use_semimajor_axis = False
        self.use_inclination = False
        self.use_time_inferior_conjunction = False
        self.use_stellar_radius = True
        self.use_stellar_temperature = True

        self.spiderman_params = None
        self.code_options = {}

    def initialize_model(self, mc, **kwargs):

        # Workarounf to allow spiderman import only when the code is axtually needed
        try:
            import spiderman
        except (ModuleNotFoundError,ImportError):
            print("ERROR: spiderman not installed, this will not work")
            quit()

        self._prepare_planetary_parameters(self, mc, **kwargs)
        self._prepare_star_parameters(self, mc, **kwargs)

        self.spiderman_params = spiderman.ModelParams(
            brightness_model='Louden', stellar_model="blackbody")

        print('   Warning on Spiderman Thermal model: null limb darkening parameters for the planet')
        print('   Warning on Spiderman Thermal model: null T_int')
        self.spiderman_params.p_u1 = 0.               # Planetary limb darkening parameter
        self.spiderman_params.p_u2 = 0.               # Planetary limb darkening parameter
        self.spiderman_params.T_int = 0.

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        self._prepare_dataset_options(self, mc, dataset, **kwargs)
        self.spiderman_params.thermal = True


    def compute(self, parameter_values, dataset, x0_input=None):
        """[summary]

        Args:
            parameter_values ([type]): [description]
            dataset ([type]): [description]
            x0_input ([type], optional): [description]. Defaults to None.
        """

        """ retrieve stellar temperature and radius either from the options or
            from the priors
        """

        stellar_radius = self.retrieve_radius(parameter_values)
        self.spiderman_params.T_s = self.retrieve_temperature(parameter_values)

        self.spiderman_params.per = parameter_values['P']  # orbital period
        # planet radius (in units of stellar radii)
        self.spiderman_params.rp = parameter_values['R_Rs']
        self.spiderman_params.ecc = parameter_values['e']  # eccentricity
        # argument of periastron (in degrees)
        self.spiderman_params.w = parameter_values['omega']

        self.spiderman_params.l1 = self.code_options[dataset.name_ref]['l1']
        self.spiderman_params.l2 = self.code_options[dataset.name_ref]['l2']

        # The absolute value of the semi-major axis [AU]
        self.spiderman_params.a_abs = convert_ars_to_a(self.spiderman_params.a,
                                                       stellar_radius)

        self.spiderman_params.insol = convert_RTaAU_to_insol(stellar_radius,
                                                            self.spiderman_params.T_s,
                                                            self.spiderman_params.a_abs)

        self.spiderman_params.albedo = parameter_values['albedo']
        self.spiderman_params.redist = parameter_values['redist']
        #self.spiderman_params.insol = parameter_values['insol']

        if x0_input is None:
            if self.spiderman_params.a < 1.0: return dataset.x0 * 0.

            if self.code_options[dataset.name_ref]['sample_factor'] > 1:

                time_array = np.linspace(-self.code_options[dataset.name_ref]['exp_time']/2.,
                                        self.code_options[dataset.name_ref]['exp_time']/2.,
                                        self.code_options[dataset.name_ref]['sample_factor'])

                lc_out = dataset.x0 * 0.
                for it, tt in enumerate(dataset.x0):
                    lc_out[it] = np.average(self.spiderman_params.lightcurve(tt + time_array))
                return lc_out - 1.
            else:
                return self.spiderman_params.lightcurve(dataset.x0) - 1.


        else:
            if self.spiderman_params.a < 1.0: return x0_input*0.

            if self.code_options[dataset.name_ref]['sample_factor'] > 1:

                time_array = np.linspace(-self.code_options[dataset.name_ref]['exp_time']/2.,
                                        self.code_options[dataset.name_ref]['exp_time']/2.,
                                        self.code_options[dataset.name_ref]['sample_factor'])

                lc_out = x0_input * 0.
                for it, tt in enumerate(x0_input):
                    lc_out[it] = np.average(self.spiderman_params.lightcurve(tt + time_array))
                return lc_out - 1.
            else:
                return self.spiderman_params.lightcurve(x0_input) - 1.
