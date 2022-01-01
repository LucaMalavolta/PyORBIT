from pyorbit.classes.common import np, convert_rho_to_a, convert_b_to_i
import pyorbit.classes.constants as constants
import pyorbit.classes.kepler_exo as kepler_exo
from pyorbit.models.abstract_model import AbstractModel

class RossiterMcLaughling_Ohta(AbstractModel):
    model_class = 'rossiter_mclaughlin'
    unitary_model = False

    default_bounds = {}
    default_spaces = {}
    default_priors = {}

    recenter_pams_dataset = {}

    def __init__(self, *args, **kwargs):

        super(RossiterMcLaughling_Ohta, self).__init__(*args, **kwargs)

        try:
            from PyAstronomy import modelSuite as PyAstroModelSuite
        except ImportError:
            print("ERROR: PyAstronomy not installed, this will not work")
            quit()

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = {
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'o',  # argument of pericenter (in radians)
            'R',  # planet radius (in units of stellar radii)
            'o_star', # Sky-projected angle between stellar rotation axis and normal of orbit plane [deg]
            'i_star', # Inclination of the star
            'v_sini' # projected rotational velocity of the star
        }

        self.list_pams_dataset = {
            'ld_c1', #linear limb darkening
        }

        self.parametrization = 'Standard'
        self.orbit = 'circular'
        self.use_semimajor_axis = False
        self.use_inclination = False
        self.use_time_of_transit = False

        self.rm_ohta = None
        self.multivariate_mass_radius = False

    def initialize_model(self, mc, **kwargs):

        try:
            multivariate_vars = mc.common_models[self.stellar_ref].multivariate_vars
        except AttributeError:
            multivariate_vars = []

        if mc.common_models[self.planet_ref].use_semimajor_axis:
            """ a is the semi-major axis (in units of stellar radii) """
            self.list_pams_common.update({'a': None})
            self.use_semimajor_axis = True
        else:
            if 'mass' in multivariate_vars and 'radius' in multivariate_vars:
                self.list_pams_common.update({'mass': None, 'radius':None})
                self.multivariate_mass_radius = True
            else:
                """ rho is the density of the star (in solar units) """
                self.list_pams_common.update({'rho': None})
                self.list_pams_common.update({'radius': None})
                self.multivariate_mass_radius = False

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

        """ The appropriate function for variable conversion is stored internally
        """
        if self.use_semimajor_axis and self.use_inclination:
            self.retrieve_ai = self._internal_transformation_mod03
        elif self.use_semimajor_axis:
            if self.multivariate_mass_radius:
                self.retrieve_ai = self._internal_transformation_mod07
            else:
                self.retrieve_ai = self._internal_transformation_mod02
        elif self.use_inclination:
            self.retrieve_ai = self._internal_transformation_mod01
        else:
            if self.multivariate_mass_radius:
                self.retrieve_ai = self._internal_transformation_mod06
            else:
                self.retrieve_ai = self._internal_transformation_mod00

        if self.use_time_of_transit:
            self.retrieve_t0 = self._internal_transformation_mod04
        else:
            self.retrieve_t0 = self._internal_transformation_mod05

        """ Depending if the orbit is circular or not, a different function
            is selected
        """
        if mc.common_models[self.planet_ref].orbit == 'circular':
            self.orbit = 'circular'
            self.rm_ohta = PyAstroModelSuite.RmcL()
        else:
            self.orbit = 'keplerian'
            self.rm_ohta = PyAstroModelSuite.RmcLell()



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

        var_omega = variable_value['o'] * (180. / np.pi)


        Omega = variable_value['v_sini'] / (variable_value['radius'] * constants.Rsun) / np.cos(variable_value['i_star']/180.*np.pi)

        if self.orbit == 'circular':
            self.rm_ohta.assignValue({"a": var_a, 
                            "lambda": variable_value['o_star'], 
                            "epsilon": variable_value['ld_c1'],
                            "P": variable_value['P'],
                            "T0": var_tc, 
                            "i": var_i/180.*np.pi,
                            "Is": variable_value['i_star']/180.*np.pi, 
                            "Omega": Omega, 
                            "gamma": variable_value['R']})
        else:

            if self.use_time_of_transit:
                Tperi  = kepler_exo.kepler_Tc2Tperi_Tref(variable_value['P'],
                                                         var_tc,
                                                         variable_value['e'],
                                                         variable_value['o'])
            else:
                Tperi  = kepler_exo.kepler_phase2Tperi_Tref(variable_value['P'],
                                                         variable_value['f'],
                                                         variable_value['e'],
                                                         variable_value['o'])

            self.rm_ohta.assignValue({"a": var_a,
                "lambda": variable_value['o_star'],
                "epsilon": variable_value['ld_c1'],
                "P": variable_value['P'],
                "tau": Tperi,
                "i": variable_value['i']/180.*np.pi,
                "w": variable_value['o']-np.pi,
                "e":variable_value['e'],
                "Is": variable_value['i_star']/180.*np.pi,
                "Omega": Omega,
                "gamma": variable_value['R']})


        if x0_input is None:
            return self.rm_ohta(time) * variable_value['radius'] * constants.Rsun * 1000.
        else:
            return self.rm_ohta(x0_input) * variable_value['radius'] * constants.Rsun * 1000.

    """ function for internal transformation of variables, to avoid if calls
        copied & past from PyTransit_Transit class
    """
    def _internal_transformation_mod00(self, variable_value):
        """ this function transforms b and rho to i and a  """
        a = convert_rho_to_a(variable_value['P'], variable_value['rho'])
        i = convert_b_to_i(
            variable_value['b'], variable_value['e'], variable_value['o'], a)
        return a, i

    def _internal_transformation_mod01(self, variable_value):
        """ this function transforms b to i"""
        i = convert_b_to_i(
            variable_value['b'], variable_value['e'], variable_value['o'], variable_value['a'])
        return variable_value['a'], i

    def _internal_transformation_mod02(self, variable_value):
        """ this function transforms rho to a  """
        a = convert_rho_to_a(variable_value['P'], variable_value['rho'])
        return a, variable_value['i']

    def _internal_transformation_mod03(self, variable_value):
        """ no transformation needed  """
        return variable_value['a'], variable_value['i']

    def _internal_transformation_mod04(self, variable_value, Tref):
        """ this function transforms Tc into Tc- Tref t"""
        return variable_value['Tc'] - Tref

    def _internal_transformation_mod05(self, variable_value, Tref):
        """ this function transforms Tc into Tc- Tref t"""
        return kepler_exo.kepler_phase2Tc_Tref(variable_value['P'],
                                               variable_value['f'],
                                               variable_value['e'],
                                               variable_value['o'])

    def _internal_transformation_mod06(self, variable_value):
        """ this function transforms b, mass, radius to i and a 
            it replaces _internal_transformation_mod00 when mass & radius
            multivariate are used
        """
        rho = variable_value['mass']/variable_value['radius']**3
        a = convert_rho_to_a(variable_value['P'], rho)
        i = convert_b_to_i(
            variable_value['b'], variable_value['e'], variable_value['o'], a)
        return a, i

    def _internal_transformation_mod07(self, variable_value):
        """ this function transforms P,mass, radius to a
            it replaces _internal_transformation_mod02 when mass & radius
            multivariate are used
        """
        rho = variable_value['mass']/variable_value['radius']**3
        a = convert_rho_to_a(variable_value['P'], rho)
        return a, variable_value['i']