from pyorbit.subroutines.common import np, OrderedSet
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import AbstractTransit

try:
    from PyAstronomy import modelSuite as PyAstroModelSuite
except (ModuleNotFoundError,ImportError):
    pass

class RossiterMcLaughlin_Ohta(AbstractModel, AbstractTransit):
    model_class = 'rossiter_mclaughlin'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)  # this calls all constructors up to AbstractModel
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            from PyAstronomy import modelSuite as PyAstroModelSuite
        except (ModuleNotFoundError,ImportError):
            print("ERROR: PyAstronomy not installed, this will not work")
            quit()

        self.unitary_model = False

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = OrderedSet([
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'omega',  # argument of pericenter (in radians)
            'lambda', # Sky-projected angle between stellar rotation axis and normal of orbit plane [deg]
            'R_Rs',  # planet radius (in units of stellar radii)
        ])

        self.Omega_rotation_conversion = None

        self.rm_ohta = None

    def initialize_model(self, mc, **kwargs):

        mc.common_models[self.stellar_ref].use_stellar_inclination = True
        if not (mc.common_models[self.stellar_ref].use_projected_velocity
                or mc.common_models[self.stellar_ref].use_stellar_rotation_period):
            mc.common_models[self.stellar_ref].use_equatorial_velocity = True
            mc.common_models[self.stellar_ref].use_stellar_radius = True

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
            print('    this model accepts only linear limb-darkening coefficients')
            print()

        if mc.common_models[self.stellar_ref].use_equatorial_velocity \
            and mc.common_models[self.stellar_ref].use_stellar_radius:
                print("Angular rotation velocity from equatorial velocity and radius")
                self.Omega_rotation_conversion = self._omega_conversion_mod00
        elif self.compute_Omega_rotation and mc.common_models[self.stellar_ref].use_stellar_rotation:
                print("Angular rotation velocity from rotation period")
                self.Omega_rotation_conversion = self._omega_conversion_mod01
        else:
                print("Angular rotation velocity from vsini")
                self.Omega_rotation_conversion = self._omega_conversion_mod02



    def compute(self, parameter_values, dataset, x0_input=None):
        """
        :param parameter_values:
        :param dataset:
        :param x0_input:
        :return:
        """
        #t1_start = process_time()

        self.update_parameter_values(parameter_values, dataset.Tref)
        parameter_values['Omega_rotation'] = self.Omega_rotation_conversion(parameter_values)

        for key, key_val in parameter_values.items():
            if np.isnan(key_val):
                return 0.

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

            if self.use_time_inferior_conjunction:
                Tperi  = kepler_exo.kepler_Tc2Tperi_Tref(parameter_values['P'],
                                                         parameter_values['Tc'] - dataset.Tref,
                                                         parameter_values['e'],
                                                         parameter_values['omega'])
            else:
                Tperi  = kepler_exo.kepler_phase2Tperi_Tref(parameter_values['P'],
                                                         parameter_values['mean_long'],
                                                         parameter_values['e'],
                                                         parameter_values['omega'])

            self.rm_ohta.assignValue({"a": parameter_values['a_Rs'],
                "lambda": parameter_values['lambda']*constants.deg2rad,
                "epsilon": parameter_values['ld_c1'],
                "P": parameter_values['P'],
                "tau": Tperi,
                "i": parameter_values['i'] *constants.deg2rad,
                "w": parameter_values['omega']*constants.deg2rad,
                "e":parameter_values['e'],
                "Is": parameter_values['i_star']*constants.deg2rad,
                "Omega": parameter_values['Omega_rotation'],
                "gamma": parameter_values['R_Rs']})

        if x0_input is None:
            return self.rm_ohta.evaluate(dataset.x0) * parameter_values['radius'] * constants.Rsun * 1000.
        else:
            return self.rm_ohta.evaluate(x0_input) * parameter_values['radius'] * constants.Rsun * 1000.


    def _omega_conversion_mod00(self, parameter_values):
        return parameter_values['veq_star'] / (parameter_values['radius'] * constants.Rsun)

    def _omega_conversion_mod01(self, parameter_values):
        return 2* np.pi / ( parameter_values['rotation_period'] * constants.d2s)

    def _omega_conversion_mod02(self, parameter_values):
        return parameter_values['v_sini'] / (parameter_values['radius'] * constants.Rsun) / np.sin(parameter_values['i_star'] * constants.deg2rad)
