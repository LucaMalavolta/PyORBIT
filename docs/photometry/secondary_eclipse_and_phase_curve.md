# Secondary eclipse and phase curve

The `batman_transit_eclipse_phasecurve` model combines the primary transit, the
secondary eclipse and a sinusoidal orbital phase curve.  It is a photometric
model for systems where the same linear ephemeris must describe both the
transit and the occultation.

The alias `batman_transit_secondary_phasecurve` is also accepted.

## Model definition and requirements

- model name: `batman_transit_eclipse_phasecurve`
- alias: `batman_transit_secondary_phasecurve`

Use this model when the light curve contains the primary transit and enough
out-of-transit orbital coverage to constrain the planet-to-star flux ratio.  The
model uses the same orbital parameters as the standard transit model, then adds
dataset-level parameters for the eclipse depth and the phase-curve amplitude.

For a pure transit model, use [Transit light-curve models](transit_light_curve_models.md).
For independent mid-transit times, use
[Transit models for TTV measurements](transit_models_for_ttv_measurements.md).

## Model parameters

| Parameter | Scope | Meaning |
| --- | --- | --- |
| `P` | planet common | Orbital period. |
| `Tc` | planet common | Reference time of inferior conjunction. |
| `R_Rs` | planet common | Planet-to-star radius ratio. |
| `e`, `omega` | planet common | Eccentricity and argument of periastron, or the selected parametrization. |
| `phase_off` | planet common | Phase offset in degrees, included when `phase_offset: True`. |
| `phase_amp` | dataset | Phase-curve semi-amplitude. |
| `delta_occ` | dataset | Secondary-eclipse depth, included when `nightside_emission: True`. |
| limb-darkening coefficients | limb-darkening common | Coefficients used by the transit part of the `batman` model. |

The model also uses the transit geometry prepared by the abstract transit
machinery: `b` or `i`, and either `a_Rs` or the stellar density, depending on
the common-object configuration.

## Keywords

Model-wide keywords, with default values in boldface where applicable.

**planets**
* accepted values: list of planet common objects (required)
* selects the planets included in the model.

**limb_darkening**
* accepted values: limb-darkening common-object name (required)
* selects the limb-darkening common object.

**supersample_factor**
* accepted values: positive integer | **`1`**
* sets the number of sub-exposures used to integrate finite exposure times.

**exposure_time**
* accepted values: exposure time in seconds | **`0.01`**
* sets the exposure time used with `supersample_factor`; use the dataset cadence when it represents the exposure duration.

**nthreads**
* accepted values: positive integer | **`1`**
* sets the number of threads passed to `batman`.

**nightside_emission**
* accepted values: **`True`** | `False`
* if `False`, `delta_occ` is removed and the occulted flux is tied to `phase_amp`.

**phase_offset**
* accepted values: **`True`** | `False`
* if `False`, `phase_off` is removed and the phase maximum is fixed to the default phase.

## Examples

```yaml
planet_b:
  common: planets
  orbit: circular
  parametrization: Eastman2013
  use_time_inferior_conjunction: True

phasecurve_model:
  model: batman_transit_eclipse_phasecurve
  planets: [b]
  limb_darkening: ld_quadratic
  supersample_factor: 5
  exposure_time: 0.02043365
  nightside_emission: True
  phase_offset: True

input:
  LCdata:
    file: full_orbit_lightcurve.dat
    kind: Phot
    models:
      - phasecurve_model
```

When `nightside_emission` is enabled, `delta_occ` represents the flux blocked
during the secondary eclipse.  The phase curve is normalized internally with
`phase_amp` and `delta_occ`, so the priors on these two parameters should be
chosen consistently with the flux units of the photometric dataset.
