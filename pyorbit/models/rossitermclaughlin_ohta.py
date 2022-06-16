from pyorbit.subroutines.common import np
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import AbstractTransit

try:
    from PyAstronomy import modelSuite as PyAstroModelSuite
except ImportError:
    pass

class RossiterMcLaughling_Ohta(AbstractModel, AbstractTransit):
    model_class = 'rossiter_mclaughlin'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)  # this calls all constructors up to AbstractModel
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            from PyAstronomy import modelSuite as PyAstroModelSuite
        except ImportError:
            print("ERROR: PyAstronomy not installed, this will not work")
            quit()

        self.unitary_model = False

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = {
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'omega',  # argument of pericenter (in radians)
            'lambda', # Sky-projected angle between stellar rotation axis and normal of orbit plane [deg]
            'R_Rs',  # planet radius (in units of stellar radii)
            'v_sini' # projected rotational velocity of the star
        }

        self.use_stellar_radius = True
        self.use_stellar_period = True
        self.use_stellar_inclination = False

        self.rm_ohta = None

    def initialize_model(self, mc, **kwargs):

        self._prepare_planetary_parameters(mc, **kwargs)
        self._prepare_star_parameters(mc, **kwargs)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

        """ Depending if the orbit is circular or not, a different function
            is selected
        """
        if mc.common_models[self.planet_ref].orbit == 'circular':
            self.orbit = 'circular'
            self.rm_ohta = PyAstroModelSuite.RmcL()
        else:
            self.orbit = 'keplerian'
            self.rm_ohta = PyAstroModelSuite.RmcLell()

        if len(self.ld_vars) > 1:
            print('WARNING on rossiter_mclaughlin ohta model:  ')
            print(' this model accepts only linear limb-darkening coefficients')
            print()

    def compute(self, variable_value, dataset, x0_input=None):
        """
        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """
        #t1_start = process_time()

        var_a, var_i = self.retrieve_ai(variable_value)
        var_tc = self.retrieve_t0(variable_value, dataset.Tref)
        var_Omega, var_Is = self.retrieve_Omega_Istar(variable_value)

        if self.orbit == 'circular':
            self.rm_ohta.assignValue({"a": var_a,
                            "lambda": variable_value['lambda']/180.*np.pi,
                            "epsilon": variable_value['ld_c1'],
                            "P": variable_value['P'],
                            "T0": var_tc,
                            "i": var_i/180.*np.pi,
                            "Is": var_Is/180.*np.pi,
                            "Omega": var_Omega,
                            "gamma": variable_value['R_Rs']})
        else:

            if self.use_time_of_transit:
                Tperi  = kepler_exo.kepler_Tc2Tperi_Tref(variable_value['P'],
                                                         var_tc,
                                                         variable_value['e'],
                                                         variable_value['omega'])
            else:
                Tperi  = kepler_exo.kepler_phase2Tperi_Tref(variable_value['P'],
                                                         variable_value['mean_long'],
                                                         variable_value['e'],
                                                         variable_value['omega'])

            self.rm_ohta.assignValue({"a": var_a,
                "lambda": variable_value['lambda']/180.*np.pi,
                "epsilon": variable_value['ld_c1'],
                "P": variable_value['P'],
                "tau": Tperi,
                "i": var_i/180.*np.pi,
                "w": variable_value['omega']/180.*np.pi-np.pi,
                "e":variable_value['e'],
                "Is": var_Is/180.*np.pi,
                "Omega": var_Omega,
                "gamma": variable_value['R_Rs']})

        if x0_input is None:
            return self.rm_ohta.evaluate(dataset.x0) * variable_value['radius'] * constants.Rsun * 1000.
        else:
            return self.rm_ohta.evaluate(x0_input) * variable_value['radius'] * constants.Rsun * 1000.

