(star)=

# Star: stellar parameters

The `star` container holds stellar common objects. Its `star_parameters` entry
stores quantities shared by transit, radial-velocity, astrometric, and
stellar-activity models. A model uses only the parameters it requires.

```{note}
Define only the parameters whose defaults you need to override. Parameters
omitted from the YAML configuration retain the bounds, priors, spaces, and
fixed values declared in `CommonStarParameters`.
```

## Model definition and requirements

- container name: `star`
- common object name: `star_parameters`
- source class: `CommonStarParameters`
- requirements: place the object under `common: star:`; selected models use
  only the stellar quantities they need. Models that share stellar-activity
  hyperparameters can refer to `star_parameters` explicitly in their `common`
  list.

In a configuration file the object is defined as:

```yaml
common:
  star:
    star_parameters:
      priors:
        mass: ['Gaussian', 0.806, 0.048]
```

## Model parameters

| Name | Parameter | Unit |
| :--- | :-------- | :--- |
| `radius` | Stellar radius | Solar radii |
| `mass` | Stellar mass | Solar masses |
| `density` | Mean stellar density | Solar mean densities |
| `i_star` | Stellar spin inclination | degrees |
| `cosi_star` | Cosine of the stellar spin inclination | unitless |
| `v_sini` | Projected stellar rotational velocity | km/s |
| `rotation_period` | Stellar rotation period | days |
| `activity_decay` | Decay timescale of active regions | days |
| `temperature` | Effective temperature of the photosphere | K |
| `natural_contrast` | Intrinsic stellar-line contrast | relative depth |
| `natural_broadening` | Intrinsic stellar-line broadening | km/s |
| `rv_center` | Stellar-line velocity centroid | km/s |
| `veq_star` | Equatorial stellar rotational velocity | km/s |
| `alpha_rotation` | Differential-rotation coefficient | unitless |
| `convective_c1` | First convective-polynomial coefficient | unitless |
| `convective_c2` | Second convective-polynomial coefficient | unitless |
| `convective_c3` | Third convective-polynomial coefficient | unitless |
| `offset_ra` | Right-ascension position offset | mas |
| `offset_dec` | Declination position offset | mas |
| `pm_ra` | Right-ascension proper motion | mas/yr |
| `pm_dec` | Declination proper motion | mas/yr |
| `parallax` | Stellar parallax | mas |
| `macroturbulence` | Macroturbulent velocity | km/s |

Depending on the parametrization, `PyORBIT` derives `i_star` from
`cosi_star`; `veq_star` from `rotation_period` and `radius`; `v_sini` from
equatorial velocity and inclination; or `rotation_period` from `veq_star`
and `radius`. It also derives one of `mass`, `radius`, and `density` from the
other two. The stellar-line parameters `natural_contrast` and
`natural_broadening` belong to this common object; instrument-specific line
parameters such as `line_contrast` and `line_fwhm` belong to the relevant
model.

## Keywords

The default value is highlighted in boldface.

**use_stellar_rotation_period**
* accepted values: `True` | **`False`**
* if `True`, samples `rotation_period`, `radius`, and stellar inclination so
  that `veq_star` and `v_sini` can be derived. Do not also force
  `use_equatorial_velocity` for this parametrization.

**use_equatorial_velocity**
* accepted values: `True` | **`False`**
* includes `veq_star` as a sampled parameter.

**use_stellar_inclination**
* accepted values: `True` | **`False`**
* includes `i_star` among the sampled parameters.

**use_cosine_stellar_inclination**
* accepted values: `True` | **`False`**
* samples `cosi_star` instead of `i_star`; `i_star` is derived from it.

**use_projected_velocity**
* accepted values: **`True`** | `False`
* includes `v_sini` as a direct parameter when a selected model requires it.

**use_differential_rotation**
* accepted values: `True` | **`False`**
* enables `alpha_rotation`. Without the rotation-period parametrization,
  `veq_star` and stellar inclination are also required.

**use_stellar_radius**
* accepted values: `True` | **`False`**
* includes `radius` as a direct parameter when the chosen stellar-rotation
  parametrization requires it.

**compute_mass**
* accepted values: **`True`** | `False`
* derives `mass` from `radius` and `density`. If both `mass` and `radius`
  have priors but `density` does not, `PyORBIT` selects `compute_density`
  automatically unless explicitly overridden.

**compute_radius**
* accepted values: `True` | **`False`**
* derives `radius` from `mass` and `density`.

**compute_density**
* accepted values: `True` | **`False`**
* derives `density` from `mass` and `radius`.

```{warning}
Select at most one of `compute_mass`, `compute_radius`, and `compute_density`.
If all three are set to `False`, `PyORBIT` falls back to
`compute_mass: True`.
```

**convective_order**
* accepted values: **`0`** | `1` | `2` | `3`
* includes `convective_c1`, `convective_c2`, and `convective_c3` up to the
  requested order in models that support convective-polynomial terms.

## Examples

The stellar mass, radius, and density priors below come from
[HD189733_example05_TESSandRV.yaml](https://github.com/LucaMalavolta/PyORBIT_examples/blob/main/quickstart/HD189733_example05_TESSandRV.yaml).
The selected transit and radial-velocity models use the stellar quantities
they need:

```yaml
common:
  star:
    star_parameters:
      priors:
        mass: ['Gaussian', 0.806, 0.048]
        radius: ['Gaussian', 0.756, 0.018]
        density: ['Gaussian', 1.864, 0.175]
```

The two-season activity fit in
[RV_GPtrained_2seasons.yaml](https://github.com/LucaMalavolta/PyORBIT_examples/blob/main/gp_multiple_seasons/RV_GPtrained_2seasons.yaml) supplies a stellar rotation
prior to two Gaussian-process models. Here `use_stellar_rotation_period` is a
keyword of each *activity model*: it makes that model use the
`star_parameters.rotation_period` value rather than the activity object's
`Prot`. The same-named keyword under `common: star: star_parameters:` instead
controls how the stellar rotational velocities are parametrized.

```yaml
common:
  activity_s01:
    model: activity
    boundaries:
      Pdec: [10.0, 100.0]
      Oamp: [0.001, 1.0]
  star:
    star_parameters:
      boundaries:
        rotation_period: [8.0, 10.0]
      priors:
        rotation_period: ['Gaussian', 8.8, 0.1]
models:
  gp_quasiperiodic_s01:
    model: gp_quasiperiodic
    common: [activity_s01, star_parameters]
    use_stellar_rotation_period: True
```

For an astrometric fit,
[HD5388_test03.yaml](https://github.com/LucaMalavolta/PyORBIT_examples/blob/main/astrometry/HD5388_test03.yaml) also provides a parallax prior:

```yaml
common:
  star:
    star_parameters:
      priors:
        mass: ['Gaussian', 1.21, 0.05]
        radius: ['Gaussian', 1.91, 0.05]
        parallax: ['Gaussian', 30.56, 0.09]
```
