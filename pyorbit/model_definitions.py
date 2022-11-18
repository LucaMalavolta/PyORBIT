#from pyorbit.common.dataset import Dataset
from pyorbit.common.planets import CommonPlanets
from pyorbit.common.activity import CommonActivity

from pyorbit.common.polynomial_trend import CommonPolynomialTrend
from pyorbit.common.common_offset import CommonOffset
from pyorbit.common.common_jitter import CommonJitter

from pyorbit.common.harmonics import CommonHarmonics

from pyorbit.common.dilution_factor import CommonDilutionFactor
from pyorbit.common.normalization_factor import CommonNormalizationFactor
from pyorbit.common.star_parameters import CommonStarParameters


from pyorbit.common.limb_darkening import LimbDarkening_Linear, \
    LimbDarkening_Quadratic, \
    LimbDarkening_SquareRoot, LimbDarkening_Logarithmic, \
    LimbDarkening_Exponential, LimbDarkening_Power2, \
    LimbDarkening_NonLinear

from pyorbit.common.cheops_modelling import CommonCheopsModelling
from pyorbit.common.correlation import CommonCorrelation

from pyorbit.common.lightcurve_detrending import CommonLightcurveDetrending

from pyorbit.models.radial_velocities import \
    RVkeplerian, RVdynamical, \
    TransitTimeKeplerian, TransitTimeDynamical, DynamicalIntegrator

from pyorbit.models.harmonics import Harmonics
from pyorbit.models.pytransit_transit import PyTransit_Transit
from pyorbit.models.pytransit_transit_with_ttv import PyTransit_Transit_With_TTV

from pyorbit.models.batman_transit import Batman_Transit
from pyorbit.models.batman_transit_with_ttv import Batman_Transit_With_TTV
from pyorbit.models.batman_transit_ttv_subset import Batman_Transit_TTV_Subset
from pyorbit.models.batman_transit_secondary_phasecurve import Batman_Transit_Eclipse_PhaseCurve
from pyorbit.models.spiderman_thermal import Spiderman_Thermal

from pyorbit.models.gp_framework_quasiperiodic_activity import \
    GP_Framework_QuasiPeriodicActivity

from pyorbit.models.gp_multidimensional_quasiperiodic_activity import \
    GP_Multidimensional_QuasiPeriodicActivity

from pyorbit.models.gp_pyaneti_quasiperiodic_activity import \
    GP_Pyaneti_QuasiPeriodicActivity

from pyorbit.models.gp_quasiperiodic_activity import \
    GaussianProcess_QuasiPeriodicActivity
from pyorbit.models.gp_quasiperiodic_activity_common import \
    GaussianProcess_QuasiPeriodicActivity_Common
from pyorbit.models.gp_quasiperiodic_activity_alternative import \
    GaussianProcess_QuasiPeriodicActivity_Alternative
from pyorbit.models.gp_quasiperiodic_cosine_activity import \
    GaussianProcess_QuasiPeriodicCosineActivity

from pyorbit.models.gp_quasiperiodic_activity_derivative import \
    GaussianProcess_QuasiPeriodicActivity_Derivative

from pyorbit.models.celerite2_granulation_rotation import Celerite2_Granulation_Rotation
from pyorbit.models.celerite2_rotation import Celerite2_Rotation
from pyorbit.models.celerite2_matern32 import Celerite2_Matern32

from pyorbit.models.celerite_rotation import Celerite_Rotation
from pyorbit.models.celerite_matern32 import Celerite_Matern32
from pyorbit.models.celerite_matern32_common import Celerite_Matern32_Common

from pyorbit.models.correlation import LocalCorrelation
from pyorbit.models.correlated_jitter import LocalCorrelatedJitter
from pyorbit.models.common_offset import Offset
from pyorbit.models.common_jitter import Jitter
from pyorbit.models.sinusoid_common_period import SinusoidCommonPeriod
from pyorbit.models.harmonics import Harmonics

from pyorbit.models.dilution_factor import DilutionFactor, LocalDilutionFactor
from pyorbit.models.normalization_factor import NormalizationFactor, LocalNormalizationFactor, SubsetNormalizationFactor

from pyorbit.models.rossitermclaughlin_ohta import RossiterMcLaughling_Ohta

from pyorbit.models.polynomial_trend import PolynomialTrend, LocalPolynomialTrend, SubsetPolynomialTrend

from pyorbit.models.spectral_rotation import SpectralRotation, SubsetSpectralRotation, SubsetSpectralRotationPolynomial

from pyorbit.models.cheops_detrending import CheopsDetrending
from pyorbit.models.cheops_factormodel import CheopsFactorModel

from pyorbit.models.lightcurve_linear_detrending import LightcurveLinearDetrending, LocalLightcurveLinearDetrending
from pyorbit.models.lightcurve_detrending import LightcurveDetrending, LocalLightcurveDetrending
from pyorbit.models.lightcurve_poly_detrending import LightcurvePolyDetrending, LocalLightcurvePolyDetrending

from pyorbit.models.celerite2_granulation_oscillation_rotation import Celerite2_Granulation_Oscillation_Rotation
from pyorbit.models.tinygp_quasiperiodic_activity import TinyGaussianProcess_QuasiPeriodicActivity
from pyorbit.models.tinygp_multidimensional_quasiperiodic_activity import TinyGP_Multidimensional_QuasiPeriodicActivity

"""
 model_requires_planets: all those models that requires AT LEAST one of the planets in the system must be listed here
    this is the case for dataset that contains the signature of multiple planets, e.g., RVs or transit light curve
 single_planet_model: the model is associated to a specific planet, e.g., time of transits
"""

model_requires_planets = ['radial_velocities' 'rv_planets',
                           'batman_transit', 'pytransit_transit',
                           'batman_transit_with_ttv', 'pytransit_transit_with_ttv',
                           'subset_batman_transit_ttv', 'batman_transit_ttv_subset',
                          'rossitermclaughlin_ohta',
                          'spiderman_thermal', 'batman_transit_eclipse_phasecurve']
single_planet_model = ['Tc_planets', 'transit_times']
transit_time_model = ['Tc_planets', 'transit_times']
model_requires_limbdarkening = ['transit',
                                'transit_eclipse_phasecurve',
                                'spectral_rotation',
                                'subset_spectral_rotation',
                                'subset_spectral_rotation_polynomial']

define_common_type_to_class = {
    'planets': CommonPlanets,
    'activity': CommonActivity,
    'polynomial_trend': CommonPolynomialTrend,
    'common_offset': CommonOffset,
    'common_jitter': CommonJitter,
    'harmonics': CommonHarmonics,
    'ld_linear': LimbDarkening_Linear,
    'ld_quadratic': LimbDarkening_Quadratic,
    'ld_square-root': LimbDarkening_SquareRoot,
    'ld_logarithmic': LimbDarkening_Logarithmic,
    'ld_exponential': LimbDarkening_Exponential,
    'ld_power2': LimbDarkening_Power2,
    'ld_nonlinear': LimbDarkening_NonLinear,
    'dilution_factor': CommonDilutionFactor,
    'normalization_factor': CommonNormalizationFactor,
    'star_parameters': CommonStarParameters,
    'cheops_modelling': CommonCheopsModelling,
    'correlation': CommonCorrelation,
    'lightcurve_detrending': CommonLightcurveDetrending,
}

define_type_to_class = {
    'radial_velocities': {'circular': RVkeplerian,
                          'keplerian': RVkeplerian,
                          'dynamical': RVdynamical},
    'rv_planets': {'circular': RVkeplerian,
                   'keplerian': RVkeplerian,
                   'dynamical': RVdynamical},
    'Tc_planets': {'circular': TransitTimeKeplerian,
                   'keplerian': TransitTimeKeplerian,
                   'dynamical': TransitTimeDynamical},
    'transit_times': {'circular': TransitTimeKeplerian,
                      'keplerian': TransitTimeKeplerian,
                      'dynamical': TransitTimeDynamical},
    'pytransit_transit': PyTransit_Transit,
    'pytransit_transit_with_ttv': PyTransit_Transit_With_TTV,
    'batman_transit': Batman_Transit,
    'batman_transit_ttv_subset': Batman_Transit_TTV_Subset,
    'subset_batman_transit_ttv': Batman_Transit_TTV_Subset,
    'batman_transit_with_ttv': Batman_Transit_With_TTV,
    'batman_transit_eclipse_phasecurve': Batman_Transit_Eclipse_PhaseCurve,
    'spiderman_thermal': Spiderman_Thermal,
    'gp_framework_quasiperiodic': GP_Framework_QuasiPeriodicActivity,
    'gp_pyaneti_quasiperiodic': GP_Pyaneti_QuasiPeriodicActivity,
    'gp_multidimensional_quasiperiodic': GP_Multidimensional_QuasiPeriodicActivity,
    'gp_quasiperiodic': GaussianProcess_QuasiPeriodicActivity,
    'gp_quasiperiodic_common': GaussianProcess_QuasiPeriodicActivity_Common,
    'gp_quasiperiodic_alternative': GaussianProcess_QuasiPeriodicActivity_Alternative,
    'gp_quasiperiodic_derivative': GaussianProcess_QuasiPeriodicActivity_Derivative,
    'gp_quasiperiodic_cosine': GaussianProcess_QuasiPeriodicCosineActivity,
    'celerite2_granulation_rotation': Celerite2_Granulation_Rotation,
    'celerite2_granulation_oscillation_rotation': Celerite2_Granulation_Oscillation_Rotation,
    'celerite2_rotation': Celerite2_Rotation,
    'celerite2_matern32': Celerite2_Matern32,
    'celerite_matern32': Celerite_Matern32,
    'celerite_matern32_common': Celerite_Matern32_Common,
    'celerite_rotation': Celerite_Rotation,
    'local_correlation': LocalCorrelation,
    'correlation': LocalCorrelation,
    'polynomial_trend': PolynomialTrend,
    'local_polynomial_trend': LocalPolynomialTrend,
    'subset_polynomial_trend': SubsetPolynomialTrend,
    'polynomial_trend_subset': SubsetPolynomialTrend,
    'common_offset': Offset,
    'common_jitter': Jitter,
    'harmonics': Harmonics,
    'sinusoid_common_period': SinusoidCommonPeriod,
    'dilution_factor': DilutionFactor,
    'local_dilution_factor': LocalDilutionFactor,
    'normalization_factor': NormalizationFactor,
    'local_normalization_factor': LocalNormalizationFactor,
    'subset_normalization_factor': SubsetNormalizationFactor,
    'local_correlated_jitter': LocalCorrelatedJitter,
    'rossitermclaughlin_ohta': RossiterMcLaughling_Ohta,
    'spectral_rotation': SpectralRotation,
    'subset_spectral_rotation': SubsetSpectralRotation,
    'subset_spectral_rotation_polynomial': SubsetSpectralRotationPolynomial,
    'cheops_detrending': CheopsDetrending,
    'cheops_factormodel': CheopsFactorModel,
    'lightcurve_linear_detrending': LightcurveLinearDetrending,
    'local_lightcurve_detrending': LocalLightcurveLinearDetrending,
    'lightcurve_detrending': LightcurveDetrending,
    'local_lightcurve_detrending': LocalLightcurveDetrending,
    'lightcurve_poly_detrending': LightcurvePolyDetrending,
    'local_lightcurve_poly_detrending': LocalLightcurvePolyDetrending,
    'tinygp_quasiperiodic': TinyGaussianProcess_QuasiPeriodicActivity,
    'tinygp_multidimensional_quasiperiodic': TinyGP_Multidimensional_QuasiPeriodicActivity
}

accepted_extensions = ['.yaml', '.yml', '.conf', '.config', '.input', ]

star_properties_list = ['limb_darkening', 'dilution_factor']

# Trying to guess all the possible mistakes....
datatype_definition = {
    'RV': ['RV', 'RVs', 'rv', 'rvs'],
    'Tcent': ['Tcent', 'TCent', 'Tc', 'TC', 'T0', 'TT'],
    'astrometry': ['astrometry', 'Astrometry', 'AstroMetry', 'Astro', 'astro', 'AM', 'Gaia', 'gaia'],
    'H-alpha': ['H', 'HA', 'h', 'ha', 'Halpha', 'H-alpha', 'halpha', 'h-alpha'],
    'Phot': ['P', 'Ph', 'p', 'ph', 'PHOT', 'Phot', 'phot', 'Photometry', 'photometry'],
    'FWHM': ['FWHM', 'fwhm'],
    'BIS': ['BIS', 'bis'],
    'Ca_HK': ['Ca', 'Cahk', 'CaHK', 'Ca_hk', 'Ca_HK',
              'logR', 'logRhk', 'logRHK', 'logR_hk', 'logR_HK',
              'log(R)', 'log(Rhk)', 'log(RHK)', 'log(R_hk)', 'log(R_HK)'],
    'S_index': ['S', 'S_index', 'Shk', 'SHK', 'S_HK', 'S_hk'],
    'CCF': ['CCF', 'CCFs', 'ccf', 'ccfs']
}