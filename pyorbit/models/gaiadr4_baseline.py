from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.subroutines.kepler_exo import kepler_E
import pyorbit.subroutines.constants as constants


class GaiaDR4AstrometricBaseline(AbstractModel):
    model_class = 'gaia_astrometry_baseline'
    time_independent_model = True

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        ''' Orbital parameters to be used in the astrometric fit '''
        self.list_pams_common = OrderedSet([
            'parallax', #stellar parallax
            "offset_ra",
            "offset_dec",
            "pm_ra",
            "pm_dec",
        ])

        self.list_pams_dataset = OrderedSet()

        self.scan_angle_column = "scan_angle_deg"
        self.parallax_factor_column = "parallax_factor_al"
        self.scan_angle_sign = 1.0
        self.scan_angle_offset_deg = 0.0
        self.parallax_factor_sign = 1.0

        self.gaia_instrumental = {}

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


    def compute(self, parameter_values, dataset, x0_input=None):
        """Return position + proper motion + annual parallax in AL, in mas."""
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

        return east * self.gaia_instrumental[dataset.name_ref]['sin_psi'] \
                + north * self.gaia_instrumental[dataset.name_ref]['cos_psi']\
                + parameter_values["parallax"] * self.gaia_instrumental[dataset.name_ref]['parallax_factor_al']*self.parallax_factor_sign

