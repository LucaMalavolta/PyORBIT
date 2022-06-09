from pyorbit.subroutines.common import np
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import AbstractTransit
import astropy.units as u


class RossiterMcLaughling_Starry(AbstractModel, AbstractTransit):
    model_class = 'rossiter_mclaughlin'

    def __init__(self, *args, **kwargs):
        # this calls all constructors up to AbstractModel
        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            import starry
            starry.config.lazy = False
            starry.config.quiet = True
        except ImportError:
            print("ERROR: PyAstronomy not installed, this will not work")
            quit()

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = {
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'omega',  # argument of pericenter (in radians)
            'R_Rs',  # planet radius (in units of stellar radii)
            # Sky-projected angle between stellar rotation axis and normal of orbit plane [deg]
            'o_star',
            'i_star',  # Inclination of the star
            'v_sini'  # projected rotational velocity of the star
        }

        self.use_stellar_radius = True
        self.unitary_model = False

        self.rm_ohta = None

    def initialize_model(self, mc, **kwargs):

        self._prepare_planetary_parameters(mc, **kwargs)
        self._prepare_star_parameters(self, mc, **kwargs)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

        """ Depending if the orbit is circular or not, a different function
            is selected
        """
        #if mc.common_models[self.planet_ref].orbit == 'circular':
        #    self.orbit = 'circular'
        #    self.rm_ohta = PyAstroModelSuite.RmcL()
        #else:
        #    self.orbit = 'keplerian'
        #    self.rm_ohta = PyAstroModelSuite.RmcLell()

        #if len(self.ld_vars) > 1:
        #    print('WARNING on rossiter_mclaughlin ohta model:  ')
        #    print(' this model accepts only linear limb-darkening coefficients')
        #    print()

    def compute(self, variable_value, dataset, x0_input=None):
        """
        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """
        #t1_start = process_time()

        Omega = variable_value['v_sini'] / (variable_value['radius'] *
                                            constants.Rsun) / np.cos(variable_value['i_star']/180.*np.pi)

        veq = Omega * variable_value['radius'] * constants.Rsun
        veq = variable_value['v_sini'] / \
            np.cos(variable_value['i_star']/180.*np.pi) * 1000.  # (m/s)

        # a = 6.7                       # Semi major axis [stellar radii]
        # lambda_ang = 7.2/180.0*np.pi  # Sky-projected angle between stellar rotation axis and normal of orbit plane [rad]
        # epsilon = 0.5                 # linear limb dark
        # P = 1.74                      # Orbital period [d]
        # T0 = 0.2                      # Central transit time
        # i = 87.8/180.*np.pi           # Inclination of orbit [rad]
        # Is = 90.0/180.0*np.pi         # Is - Inclination of stellar rotation axis [rad]
        # Omega = 1.609e-5             # Omega - Angular rotation velocity (star) [rad/s]
        # gamma = 0.2                   # Rp/Rs (ratio of planetary and stellar radius)

        # obl = lambda_ang
        # inc = Is

        ''' 
        a = 6.039                       # Semi major axis [stellar radii]

        # Sky-projected angle between stellar rotation axis and normal of orbit plane [rad]
        lambda_ang = 7.2/180.0*np.pi
        epsilon = 0.5                 # linear limb dark
        P = 1.74                      # Orbital period [d]
        T0 = 0.2                      # Central transit time
        i = 87.8/180.*np.pi           # Inclination of orbit [rad]
        # Is - Inclination of stellar rotation axis [rad]
        Is = 90.0/180.0*np.pi
        # Omega - Angular rotation velocity (star) [rad/s]
        Omega = 1.609e-5
        # Rp/Rs (ratio of planetary and stellar radius)
        gamma = 0.2   r_ / Rs


        A = starry.Primary(
            starry.Map(udeg=2,
                       rv=True,
                       amp=1, #The overall amplitude of the map in arbitrary units.
                       veq=veq, #equatorial velocity, in m/s
                       alpha=0, #The rotational shear coefficient, a number in the range [0, 1].
                        obl=-lambda_ang/np.pi*180,
                        inc=Is/np.pi*180),
            r=radius_star,
            m=mass_star,
            length_unit=u.Rsun,
            angle_unit=u.deg,
            mass_unit=u.Msun,
        )
        A.map[1] = epsilon
        A.map[2] = 0.85

        # Define the planet
        b = starry.Secondary(
            starry.Map(rv=True, amp=0, veq=0),
            #a = 16.7,
            r=gamma*radius_star,
            porb=P,
            m=0.0,
            t0=T0,
            inc=i,
            ecc=0.35,
            w=-np.pi/2,
            length_unit=u.Rsun,
            mass_unit=u.Msun,
            angle_unit=u.deg,
            time_unit=u.day,
        )

        # Define the system
        sys = starry.System(A, b)
        # Compute the flux & RV signal
        time = np.linspace(-0.5, 2.5, 1000)
        flux = sys.flux(time)
        rv = sys.rv(time)

        var_a, var_i = self.retrieve_ai(variable_value)
        var_tc = self.retrieve_t0(variable_value, dataset.Tref)

        Omega = variable_value['v_sini'] / (variable_value['radius'] *
                                            constants.Rsun) / np.cos(variable_value['i_star']/180.*np.pi)

        if self.orbit == 'circular':
            self.rm_ohta.assignValue({"a": var_a,
                                      "lambda": variable_value['o_star'],
                                      "epsilon": variable_value['ld_c1'],
                                      "P": variable_value['P'],
                                      "T0": var_tc,
                                      "i": var_i/180.*np.pi,
                                      "Is": variable_value['i_star']/180.*np.pi,
                                      "Omega": Omega/180.*np.pi,
                                      "gamma": variable_value['R_Rs']})
        else:

            if self.use_time_of_transit:
                Tperi = kepler_exo.kepler_Tc2Tperi_Tref(variable_value['P'],
                                                        var_tc,
                                                        variable_value['e'],
                                                        variable_value['omega'])
            else:
                Tperi = kepler_exo.kepler_phase2Tperi_Tref(variable_value['P'],
                                                           variable_value['f'],
                                                           variable_value['e'],
                                                           variable_value['omega'])

            self.rm_ohta.assignValue({"a": var_a,
                                      "lambda": variable_value['o_star'],
                                      "epsilon": variable_value['ld_c1'],
                                      "P": variable_value['P'],
                                      "tau": Tperi,
                                      "i": variable_value['i']/180.*np.pi,
                                      "w": variable_value['omega']/180.*np.pi-np.pi,
                                      "e": variable_value['e'],
                                      "Is": variable_value['i_star']/180.*np.pi,
                                      "Omega": Omega/180.*np.pi,
                                      "gamma": variable_value['R_Rs']})

        if x0_input is None:
            return self.rm_ohta(dataset.x0) * variable_value['radius'] * constants.Rsun * 1000.
        else:
            return self.rm_ohta(x0_input) * variable_value['radius'] * constants.Rsun * 1000.

    '''