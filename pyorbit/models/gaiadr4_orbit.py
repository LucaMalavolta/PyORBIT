from os import system

from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_astrometry import AbstractAstrometry
from pyorbit.subroutines.kepler_exo import kepler_E
import pyorbit.subroutines.constants as constants

class GaiaDR4AstrometricOrbit(AbstractModel, AbstractAstrometry):
    model_class = 'gaia_astrometry_orbit'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

        ''' Orbital parameters to be used in the astrometric fit '''
        self.list_pams_common = OrderedSet([
            'P',     # Period in days
            #'M_Me', # planet mass in Earth masses
            'Omega', # longitude of ascending node
            'e',     # eccentricity, uniform prior - to be fixed
            'i',
            'omega', # argument of pericenter
            #'i',     # inclination in degrees
            'mass', #stellar mass
            'parallax', #stellar parallax
        ])

        self.list_pams_dataset = OrderedSet()
        self.accept_multiple_planets = False

        self.scan_angle_column = "scan_angle_deg"
        self.parallax_factor_column = "parallax_factor_al"
        self.scan_angle_sign = 1.0
        self.scan_angle_offset_deg = 0.0
        self.parallax_factor_sign = 1.0

        self.gaia_instrumental = {}

    def print_warning(self):
        print("{0:s} WARNING:".format(self.model_name))
        print('    Astrometric orbit modelling requires the use of the stellar mass')
        print('    This may cause a clash with models requiring stellar density and radius, e.g., RM modelling')
        print('    The use of a multivariate approach is strongly suggested')
        print('    You can control the behaviour of mass/radius/density with the specific keywords')
        print('    compute_mass, compute_radius, compute_density')
        print()


    def initialize_model(self, mc, **kwargs):
        self.Tref = mc.Tref
        """Read fixed geometry conventions from the model configuration."""
        self.scan_angle_column = kwargs.get(
            "scan_angle_column", self.scan_angle_column
        )
        self.parallax_factor_column = kwargs.get(
            "parallax_factor_column", self.parallax_factor_column
        )
        self.scan_angle_sign = float(
            kwargs.get("scan_angle_sign", self.scan_angle_sign)
        )
        self.scan_angle_offset_deg = float(
            kwargs.get("scan_angle_offset_deg", self.scan_angle_offset_deg)
        )
        self.parallax_factor_sign = float(
            kwargs.get("parallax_factor_sign", self.parallax_factor_sign)
        )

    def initialize_model_parameters(self, mc, **kwargs):

        self._prepare_astrometry_parameters(mc, **kwargs)


    def initialize_model_dataset(self, mc, dataset, **kwargs):

        self.gaia_instrumental[dataset.name_ref] = {}


        if dataset.ancillary is None:
            raise ValueError(
                "Gaia AL model requires ancillary columns scan_angle_deg and parallax_factor_al"
            )

        names = dataset.ancillary.dtype.names or ()
        required = (self.scan_angle_column, self.parallax_factor_column)
        missing = [name for name in required if name not in names]
        if missing:
            raise ValueError(
                "Missing Gaia AL ancillary column(s): " + ", ".join(missing)
            )

        for name in required:
            values = np.asarray(dataset.ancillary[name], dtype=float)
            if values.shape[0] != dataset.n:
                raise ValueError(
                    f"Ancillary column {name!r} has {values.shape[0]} rows; "
                    f"dataset has {dataset.n}."
                )
            if not np.all(np.isfinite(values)):
                raise ValueError(f"Ancillary column {name!r} contains non-finite values.")

        psi_deg = np.asarray(dataset.ancillary[self.scan_angle_column], dtype=np.double)
        psi = (self.scan_angle_sign * psi_deg + self.scan_angle_offset_deg) * constants.deg2rad

        self.gaia_instrumental[dataset.name_ref]['sin_psi'] = np.sin(psi)
        self.gaia_instrumental[dataset.name_ref]['cos_psi'] = np.cos(psi)
        self.gaia_instrumental[dataset.name_ref]['parallax_factor_al'] = np.asarray(dataset.ancillary[self.parallax_factor_column]*self.parallax_factor_sign, dtype=float)

    @staticmethod
    def _eccentric_orbit_coordinates(parameter_values, x0):
        """Return the dimensionless orbital-plane coordinates X and Y.

        X = cos(E) - e
        Y = sqrt(1-e^2) sin(E)

        The mean anomaly at the common PyORBIT reference epoch is
        M0 = lambda0 - omega - Omega.
        """

        mean_anomaly = (
            2.0 * np.pi * x0 / parameter_values["P"]
            + np.deg2rad(
                parameter_values["mean_long"] - parameter_values["omega"] - parameter_values["Omega"]
            )
        ) % (2.0 * np.pi)

        eccentricity = float(parameter_values["e"])
        eccentric_anomaly = kepler_E(mean_anomaly, eccentricity)

        X = np.cos(eccentric_anomaly) - eccentricity
        Y = np.sqrt(max(0.0, 1.0 - eccentricity * eccentricity)) * np.sin(
            eccentric_anomaly
        )
        return np.asarray(X, dtype=np.double), np.asarray(Y, dtype=np.double)

    @staticmethod
    def _planet_sky_coordinates(parameter_values, X, Y):
        """Project the planet's relative orbit onto east and north.

        Omega is a position angle measured from north through east. The z axis is
        directed toward the observer, so the node at Omega is the crossing toward
        the observer. Astrometry alone has the usual 180-degree node ambiguity;
        simultaneous RVs select the physically consistent node.
        """
        omega = np.deg2rad(float(parameter_values["omega"]))
        Omega = np.deg2rad(float(parameter_values["Omega"]))
        inclination = np.deg2rad(float(parameter_values["i"]))

        sin_Omega = np.sin(Omega)
        cos_Omega = np.cos(Omega)
        sin_omega = np.sin(omega)
        cos_omega = np.cos(omega)
        cos_i = np.cos(inclination)

        # The node direction is (east, north) = (sin Omega, cos Omega).
        east_X = cos_omega * sin_Omega + sin_omega * cos_i * cos_Omega
        east_Y = -sin_omega * sin_Omega + cos_omega * cos_i * cos_Omega
        north_X = cos_omega * cos_Omega - sin_omega * cos_i * sin_Omega
        north_Y = -sin_omega * cos_Omega - cos_omega * cos_i * sin_Omega

        east_planet = east_X * X + east_Y * Y
        north_planet = north_X * X + north_Y * Y
        return east_planet, north_planet

    def compute(self, parameter_values, dataset, x0_input=None):

        if x0_input is not None:
            return np.zeros_like(x0_input, dtype=float)

        self.update_parameter_values(parameter_values)

        delta_year = dataset.x0 * constants.day2year

        X, Y = self._eccentric_orbit_coordinates(parameter_values, dataset.x0)
        east_planet, north_planet = self._planet_sky_coordinates(parameter_values, X, Y)


        """Compute the angular semimajor axis of the host-star reflex orbit.

        Kepler's third law in AU, yr, and solar masses gives

            a_rel[AU]^3 = (Mstar + Mp) * P[yr]^2.

        The host-star barycentric semimajor axis is

            a_star = a_rel * Mp / (Mstar + Mp),

        and 1 AU subtends ``parallax_mas`` milliarcseconds at the source distance.
        """

        planet_mass =  parameter_values['M_Me'] * constants.Mears
        a_star_mas = (parameter_values['a_AU'] 
            * planet_mass / (parameter_values['mass'] + planet_mass)) * parameter_values['parallax']


        # PyORBIT omega belongs to the planet. The host star is on the opposite
        # side of the system barycentre, hence the explicit minus sign.
        east_star = -a_star_mas * east_planet
        north_star = -a_star_mas * north_planet

        return east_star * self.gaia_instrumental[dataset.name_ref]['sin_psi']  + north_star * self.gaia_instrumental[dataset.name_ref]['cos_psi']