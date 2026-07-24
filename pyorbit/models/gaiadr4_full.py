from os import system

from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_astrometry import AbstractAstrometry
from pyorbit.subroutines.kepler_exo import kepler_E
import pyorbit.subroutines.constants as constants

class GaiaDR4AstrometricFull(AbstractModel, AbstractAstrometry):
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
            "offset_ra",
            "offset_dec",
            "pm_ra",
            "pm_dec",
        ])

        self.list_pams_dataset = OrderedSet()

        self.accept_multiple_planets = True
        self.force_model_dataset_initialization = True

        self.scan_angle_column = "scan_angle_deg"
        self.parallax_factor_column = "parallax_factor_al"
        self.scan_angle_sign = 1.0
        self.scan_angle_offset_deg = 0.0
        self.parallax_factor_sign = 1.0

        self.gaia_instrumental = {}

    def print_warning(self):
        pass 

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

    def print_info(self):
        print("*** model {0:s}:".format(self.model_name))
        print("WARNING:")
        print('    Astrometric orbit modelling requires the use of the stellar mass')
        print('    This may cause a clash with models requiring stellar density and radius, e.g., RM modelling')
        print('    At least TWO out of mass, radius, density of the star must be provided as priors')
        print('    The use of a multivariate approach is strongly suggested')
        print('    You can control the behaviour of mass/radius/density with the specific keywords')
        print('    compute_mass, compute_radius, compute_density')
        print()

    @staticmethod
    def _orbit_offsets(x0, period, eccentricity, omega_deg, mean_long_deg,
                       i_deg, Omega_deg, a_star_mas):
        """Return east (dalpha*) and north (ddelta) photocentre offsets."""
        omega = omega_deg * constants.deg2rad
        Omega = Omega_deg * constants.deg2rad
        inc = i_deg * constants.deg2rad

        mean_anomaly = (
            2.0 * np.pi * x0 / period
            + np.deg2rad(
                mean_long_deg - omega_deg - Omega_deg
            )
        ) % (2.0 * np.pi)

        eccentricity = float(eccentricity)
        eccentric_anomaly = kepler_E(mean_anomaly, eccentricity)

        X = np.cos(eccentric_anomaly) - eccentricity
        Y = np.sqrt(max(0.0, 1.0 - eccentricity * eccentricity)) * np.sin(
            eccentric_anomaly
        )

        sin_Omega = np.sin(Omega)
        cos_Omega = np.cos(Omega)
        sin_omega = np.sin(omega)
        cos_omega = np.cos(omega)
        cos_i = np.cos(inc)

        # Omega is measured from north through east. These expressions return
        # east=dalpha* and north=ddelta.

        # The node direction is (east, north) = (sin Omega, cos Omega).
        east_X = cos_omega * sin_Omega + sin_omega * cos_i * cos_Omega
        east_Y = -sin_omega * sin_Omega + cos_omega * cos_i * cos_Omega
        north_X = cos_omega * cos_Omega - sin_omega * cos_i * sin_Omega
        north_Y = -sin_omega * cos_Omega - cos_omega * cos_i * sin_Omega

        # PyORBIT omega belongs to the planet. The host star is on the opposite
        # side of the system barycentre, hence the explicit minus sign.
        east_planet = -a_star_mas*(east_X * X + east_Y * Y)
        north_planet = -a_star_mas*(north_X * X + north_Y * Y)
        return east_planet, north_planet

    def compute(self, not_used, dataset, x0_input=None):


        parameter_values = self.parameter_values
        if x0_input is not None:
            return np.zeros_like(x0_input, dtype=float)

        delta_year = dataset.x0 * constants.day2year

        east = (
            parameter_values["offset_ra"]
            + parameter_values["pm_ra"] * delta_year
        )
        north = (
            parameter_values["offset_dec"]
            + parameter_values["pm_dec"] * delta_year
        )

        for planet_name in self.planet_list:
            prepend = planet_name + '__'

            planet_mass =  parameter_values[prepend+'M_Me'] * constants.Mears
            a_star_mas = (self.parameter_values[prepend+'a_AU'] 
                * planet_mass / (self.parameter_values['mass'] + planet_mass)) * self.parameter_values['parallax']

            east_orbit, north_orbit = self._orbit_offsets(
                x0=dataset.x0,
                period=self.parameter_values[prepend+"P"],
                eccentricity=self.parameter_values[prepend+"e"],
                omega_deg=self.parameter_values[prepend+"omega"],
                mean_long_deg=self.parameter_values[prepend+"mean_long"],
                i_deg=self.parameter_values[prepend+"i"],
                Omega_deg=self.parameter_values[prepend+"Omega"],
                a_star_mas=a_star_mas,
            )

            east += east_orbit
            north += north_orbit

        return (east * self.gaia_instrumental[dataset.name_ref]['sin_psi'] 
                + north * self.gaia_instrumental[dataset.name_ref]['cos_psi'] 
                + parameter_values["parallax"] * self.gaia_instrumental[dataset.name_ref]['parallax_factor_al']*self.parallax_factor_sign)
