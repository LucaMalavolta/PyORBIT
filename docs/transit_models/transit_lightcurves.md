# Transit light-curve models

PyORBIT provides two main families of transit models.  The models described in
this page use one orbital ephemeris for each planet: the period `P` and one
reference time of inferior conjunction `Tc` define all the predicted transits as
`Tc + N * P`.

Use these models when the transit times are assumed to follow a linear
ephemeris, or when the timing deviations are not part of the fit.  If each
transit must have its own fitted mid-time, use the TTV models documented in
[Transit models for TTV measurements](transit_ttv.md).

## Planet setup

For a standard transit fit, define the planet as a transiting planet by enabling
the time of inferior conjunction:

```yaml
planet_b:
  common: planets
  orbit: circular
  parametrization: Eastman2013
  use_time_inferior_conjunction: True
```

With this option the planet common object exposes `Tc` as the epoch parameter.
If it is not enabled, PyORBIT uses the orbital longitude parametrization and
derives the inferior-conjunction time internally.

The transit shape is controlled by the usual planet, stellar and
limb-darkening parameters:

| Parameter | Scope | Meaning |
| --- | --- | --- |
| `P` | planet common | Orbital period. |
| `Tc` | planet common | Reference time of inferior conjunction. |
| `R_Rs` | planet common | Planet-to-star radius ratio. |
| `b` or `i` | planet common | Impact parameter or inclination, depending on the planet setup. |
| `a_Rs` or stellar `density` | planet/star common | Scaled semi-major axis, or the stellar density used to compute it. |
| `e`, `omega` | planet common | Eccentricity and argument of periastron, or the selected eccentricity parametrization. |
| limb-darkening coefficients | limb-darkening common | Coefficients used by the selected limb-darkening law. |

## Available models

| Model name | Backend | Use case |
| --- | --- | --- |
| `batman_transit` | `batman` | Standard transit model with one linear ephemeris per planet. |
| `pytransit_transit` | `PyTransit` | Standard transit model using the PyTransit backend. It uses the RoadRunner model by default when available. |
| `batman_transit_rprs_subset` | `batman` | Same linear ephemeris, but with a different `R_Rs` value for each dataset subset. |
| `pytransit_dynamical` | `PyTransit` | Dynamical transit model. Transit times are predicted by the dynamical model, not fitted as independent TTV parameters. |

The alias `subset_batman_transit_rprs` is accepted for the subset radius-ratio
model.  The transit, secondary-eclipse and phase-curve model is documented in
[Secondary eclipse and phase curve](secondary_eclipse_phasecurve.md).

## Shared model keywords

| Keyword | Models | Meaning |
| --- | --- | --- |
| `planets` | all | List of planet common objects included in the light-curve model. |
| `limb_darkening` | transit models | Limb-darkening common object used by the transit backend. |
| `supersample_factor` | all | Number of sub-exposures used to integrate long-cadence observations. |
| `exposure_time` | all | Exposure time used together with `supersample_factor`. |
| `nthreads` | `batman` models | Number of threads passed to the `batman` backend. |
| `use_roadrunner` | `pytransit_transit` | Use the PyTransit RoadRunner implementation when possible. |

## Minimal examples

A standard `batman` transit model:

```yaml
lc_model:
  model: batman_transit
  planets: [b]
  limb_darkening: ld_quadratic
  supersample_factor: 5
  exposure_time: 0.02043365
```

The equivalent PyTransit setup:

```yaml
lc_model:
  model: pytransit_transit
  planets: [b]
  limb_darkening: ld_quadratic
  use_roadrunner: True
```

The dataset using the model can then be declared in the usual way:

```yaml
input:
  LCdata:
    file: lightcurve.dat
    kind: Phot
    models:
      - lc_model
```

## Radius-ratio subsets

Use `batman_transit_rprs_subset` when all subsets share the same ephemeris and
orbital shape, but each subset needs its own radius ratio.  The input dataset
must include a `subset` column, and PyORBIT creates parameters such as `R_Rs_0`,
`R_Rs_1`, and so on for the active subset identifiers.

```yaml
lc_model:
  model: batman_transit_rprs_subset
  planets: [b]
  limb_darkening: ld_quadratic
```

This is useful for multi-band or multi-instrument light curves where the
transit depth can change, but the timing is still described by a single
`P` and `Tc`.

Choose the TTV models instead when the goal is to measure an independent
mid-transit time for each observed transit.
