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

        self.compute_Omega_rotation = True

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

    def compute(self, parameter_values, dataset, x0_input=None):
        """
        :param parameter_values:
        :param dataset:
        :param x0_input:
        :return:
        """
        #t1_start = process_time()

        self.update_parameter_values(parameter_values, dataset.Tref)

        if self.orbit == 'circular':
            self.rm_ohta.assignValue({"a": parameter_values['a_Rs'],
                            "lambda": parameter_values['lambda']/180.*np.pi,
                            "epsilon": parameter_values['ld_c1'],
                            "P": parameter_values['P'],
                            "T0": parameter_values['Tc'] - dataset.Tref,
                            "i": parameter_values['i'] *constants.deg2rad,
                            "Is": parameter_values['i_star']*constants.deg2rad,
                            "Omega": parameter_values['Omega_rotation'],
                            "gamma": parameter_values['R_Rs']})
        else:

            if self.use_time_of_transit:
                Tperi  = kepler_exo.kepler_Tc2Tperi_Tref(parameter_values['P'],
                                                         par_tc,
                                                         parameter_values['e'],
                                                         parameter_values['omega'])
            else:
                Tperi  = kepler_exo.kepler_phase2Tperi_Tref(parameter_values['P'],
                                                         parameter_values['mean_long'],
                                                         parameter_values['e'],
                                                         parameter_values['omega'])

            self.rm_ohta.assignValue({"a": par_a,
                "lambda": parameter_values['lambda']/180.*np.pi,
                "epsilon": parameter_values['ld_c1'],
                "P": parameter_values['P'],
                "tau": Tperi,
                "i": par_i/180.*np.pi,
                "w": parameter_values['omega']/180.*np.pi-np.pi,
                "e":parameter_values['e'],
                "Is": par_Is/180.*np.pi,
                "Omega": par_Omega,
                "gamma": parameter_values['R_Rs']})

        if x0_input is None:
            return self.rm_ohta.evaluate(dataset.x0) * parameter_values['radius'] * constants.Rsun * 1000.
        else:
            return self.rm_ohta.evaluate(x0_input) * parameter_values['radius'] * constants.Rsun * 1000.

