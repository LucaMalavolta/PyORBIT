(photometry)=

# Photometry

This section collects light-curve models for planetary transits, transit-time
variations, secondary eclipses and orbital phase curves, together with the
baseline and dilution models used in photometric fits.

The basic distinction is whether a planet is described by one linear ephemeris,
or whether each observed transit is assigned an independent mid-transit time.
The secondary-eclipse and phase-curve model still uses the same orbital
ephemeris as the primary transit, but adds the occultation depth and the
brightness modulation over the orbit.

```{toctree}
:maxdepth: 1
photometry/transit_light_curve_models.md
photometry/transit_models_for_ttv_measurements.md
photometry/secondary_eclipse_and_phase_curve.md
photometry/normalization_factor.md
photometry/dilution_factor.md
photometry/polynomial_normalization.md
```
