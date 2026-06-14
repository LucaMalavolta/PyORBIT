(star)=

# Star

The `star` container hosts star-related common objects. Its main entry is
`star_parameters`, which stores stellar quantities shared across models such as
transits, Rossiter-McLaughlin effects, and stellar-activity regressions.

```{note}
Only the parameters required by the selected models need to be defined in the
YAML file. `PyORBIT` will use the default bounds, priors, spaces, and fixed
values declared in the source code for any parameter that is not explicitly
overridden.
```

## Model definition

- container name: `star`
- common object name: `star_parameters`
- source class: `CommonStarParameters`

In a configuration file the object is defined as:

```yaml
common:
  star:
    star_parameters:
      ...
```

## Parameters

| Name | Parameter | Unit |
| :--- | :-------- | :--- |
| `radius` | Stellar radius | Solar radii |
| `mass` | Stellar mass | Solar masses |
| `density` | Stellar density | Solar density |
| `i_star` | Stellar inclination | degrees |
| `cosi_star` | Cosine of the stellar inclination | unitless |
| `v_sini` | Projected stellar rotational velocity | km/s |
| `rotation_period` | Stellar rotation period | days |
| `activity_decay` | Decay timescale of active regions | days |
| `temperature` | Effective temperature of the photosphere | K |
| `line_contrast` | CCF line contrast | percent |
| `line_fwhm` | CCF line full width at half maximum | km/s |
| `rv_center` | CCF line centroid | km/s |
| `veq_star` | Equatorial stellar rotational velocity | km/s |
| `alpha_rotation` | Differential-rotation coefficient | unitless |
| `convective_c1` | First convective-polynomial coefficient | unitless |
| `convective_c2` | Second convective-polynomial coefficient | unitless |
| `convective_c3` | Third convective-polynomial coefficient | unitless |

## Keywords

The default keyword is highlighted in boldface.

**use_stellar_rotation_period**
* accepted values: `True` | **`False`**
* if `True`, the stellar rotation is parametrized through `rotation_period`,
  `radius`, and stellar inclination rather than through `veq_star` or `v_sini`.
  In this setup `veq_star` and `v_sini` are derived quantities.

**use_equatorial_velocity**
* accepted values: `True` | **`False`**
* forces the use of `veq_star` as a sampled parameter.

**use_stellar_inclination**
* accepted values: `True` | **`False`**
* includes `i_star` among the sampled parameters.

**use_cosine_stellar_inclination**
* accepted values: `True` | **`False`**
* samples `cosi_star` instead of `i_star`; `i_star` is then derived from
  `cosi_star`.

**use_projected_velocity**
* accepted values: **`True`** | `False`
* includes `v_sini` as a direct parameter. This is the default behaviour when
  no rotation-period-based parametrization is requested.

**use_differential_rotation**
* accepted values: `True` | **`False`**
* activates the differential-rotation parametrization through
  `alpha_rotation`. When enabled without `use_stellar_rotation_period`, the
  code also requires `veq_star` and stellar inclination.

**compute_mass**
* accepted values: **`True`** | `False`
* when `mass` and `radius` are sampled, `density` is derived from them by
  default.

**compute_radius**
* accepted values: `True` | **`False`**
* makes `radius` a derived quantity from `mass` and `density`.

**compute_density**
* accepted values: `True` | **`False`**
* makes `density` a derived quantity from `mass` and `radius`.

```{warning}
Only one of `compute_mass`, `compute_radius`, or `compute_density` should be
active at a time. If all three are set to `False`, `PyORBIT` falls back to
`compute_mass: True`.
```

**convective_order**
* accepted values: integer, usually `0` to `3`
* enables the convective polynomial terms `convective_c1`, `convective_c2`, and
  `convective_c3` up to the requested order in Rossiter-McLaughlin-like models.

## Derived quantities

Depending on the selected keywords, `PyORBIT` can derive:

- `i_star` from `cosi_star`
- `veq_star` from `rotation_period` and `radius`
- `v_sini` from `veq_star` plus `i_star` or `cosi_star`
- `rotation_period` from `veq_star` and `radius`
- one among `mass`, `radius`, and `density` from the other two

This lets you choose the parametrization that is most natural for the dataset
being modelled while keeping the physically linked stellar quantities
consistent.

## Examples

The most common use is to provide informative priors on the stellar bulk
properties:

```yaml
common:
  star:
    star_parameters:
      priors:
        mass: ['Gaussian', 0.806, 0.048]
        radius: ['Gaussian', 0.756, 0.018]
        density: ['Gaussian', 1.864, 0.175]
```

When the stellar rotation period is known and you want `PyORBIT` to derive the
projected and equatorial velocities consistently, you can switch to the
rotation-based parametrization:

```yaml
common:
  star:
    star_parameters:
      use_stellar_rotation_period: True
      use_cosine_stellar_inclination: True
      boundaries:
        rotation_period: [10.0, 20.0]
        radius: [0.60, 0.90]
        cosi_star: [0.0, 1.0]
      priors:
        rotation_period: ['Gaussian', 14.0, 0.5]
        radius: ['Gaussian', 0.68, 0.02]
```

For Rossiter-McLaughlin analyses that need differential rotation and a simple
convective polynomial:

```yaml
common:
  star:
    star_parameters:
      use_equatorial_velocity: True
      use_stellar_inclination: True
      use_differential_rotation: True
      convective_order: 2
      boundaries:
        veq_star: [1.0, 20.0]
        i_star: [0.0, 180.0]
        alpha_rotation: [0.0, 1.0]
        convective_c1: [0.0, 2.0]
        convective_c2: [-2.0, 0.0]
```
