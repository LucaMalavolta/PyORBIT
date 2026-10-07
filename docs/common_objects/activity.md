(activity)=

# Stellar activity

The `activity` common object stores the hyperparameters used by stellar-activity
models. A model chooses which parameters it needs and whether each parameter is
shared across datasets or fitted separately for each dataset. For example,
quasi-periodic models usually share the rotation period and activity decay
timescale, while fitting an amplitude for each dataset.

```{note}
Only the parameters required by the selected models need to be defined in the
YAML file. `PyORBIT` uses the default bounds, priors, spaces, and fixed values
declared in the source code for parameters that are not explicitly overridden.
```

## Model definition and requirements

- common object name: `activity`
- source class: `CommonActivity`
- requirement: reference the common object with `common: activity` in an
  activity model under `models`. The model determines the parameters it uses.

In a configuration file the object is defined as:

```yaml
common:
  activity:
    boundaries:
      Prot: [10.0, 20.0]
      Pdec: [20.0, 1000.0]
      Oamp: [0.001, 1.0]
models:
  gp_quasiperiodic:
    model: gp_quasiperiodic
    common: activity
```

For several independent activity objects, give each one a distinct name and
declare `model: activity` inside it. When a model takes `rotation_period` or
`activity_decay` from `star_parameters`, also list `star_parameters` among that
model's common objects. The [quasi-periodic GP documentation](../gaussian_process/quasiperiodic_kernel.md)
describes the model-specific requirements.

## Model parameters

The table lists parameters available from `CommonActivity`. The selected model
uses only a subset; consult its model page to determine whether a parameter is
common or dataset-specific. “As input” follows the unit declaration in the
source code for model-dependent scales and amplitudes.

| Name | Parameter | Unit |
| :--- | :-------- | :--- |
| `Prot` | Stellar rotation period | days |
| `Pdec` | Decay timescale of active regions | days |
| `Pcyc` | Long-timescale activity cycle or squared-exponential timescale | days |
| `Oamp` | Coherence scale of the periodic component | as input |
| `Hamp` | Covariance amplitude | as input |
| `Camp` | Secondary covariance amplitude, used for a cosine, derivative, or cycle component | as input |
| `rot_sigma` | Amplitude of a celerite2 rotation term | as input |
| `rot_Q0` | Base quality factor of a celerite2 rotation term | as input |
| `rot_deltaQ` | Difference between the rotation-mode quality factors | as input |
| `rot_fmix` | Fractional amplitude of the secondary rotation mode | as input |
| `grn_period` | Granulation SHO timescale | days |
| `grn_sigma` | Granulation SHO amplitude | as input |
| `grn_k*_period` | Timescale of the `k`-th granulation SHO term | days |
| `grn_k*_sigma` | Amplitude of the `k`-th granulation SHO term | as input |
| `osc_k*_period` | Timescale of the `k`-th oscillation SHO term | days |
| `osc_k*_sigma` | Amplitude of the `k`-th oscillation SHO term | as input |
| `osc_k*_Q0` | Quality factor of the `k`-th oscillation SHO term | as input |
| `sho_scale` | SHO undamped period | days |
| `sho_decay` | SHO damping timescale | as input |
| `sho_sigma` | SHO process amplitude | as input |
| `matern32_scale` | Matérn-3/2 scale | as input |
| `matern32_sigma` | Matérn-3/2 covariance amplitude | as input |
| `rot_amp` | Coefficient of the derivative component of a latent GP | as input |
| `con_amp` | Coefficient of the latent GP itself | as input |
| `cos_amp` | Coefficient of the cosine latent component | as input |
| `cos_der` | Coefficient of the derivative of the cosine component | as input |
| `cyc_amp` | Coefficient of the squared-exponential cycle component | as input |
| `cyc_der` | Coefficient of the derivative of the cycle component | as input |
| `matern32_multigp_sigma` | Coefficient of the Matérn-3/2 latent GP | as input |
| `matern32_multigp_sigma_deriv` | Coefficient of the derivative of the Matérn-3/2 latent GP | as input |
| `Vc`, `Vr`, `Lc`, `Bc`, `Br` | Rajpaul-framework coefficients | as input, except `Vc` (days) |
| `sin_P` | Period of a sinusoidal activity term | days |
| `sin_K` | Semi-amplitude of a sinusoidal activity term | as input |
| `sin_f` | Phase of a sinusoidal activity term | degrees |

In the indexed granulation and oscillation names, `*` stands for an integer
from `0` to `9`, for example `grn_k0_period`.

## Keywords

The default keyword is highlighted in boldface. These switches are useful
when the selected activity model supports sharing stellar parameters with
[`star_parameters`](star.md).

**use_stellar_rotation_period**
* accepted values: `True` | **`False`**
* if `True`, replaces `Prot` with `rotation_period` from `star_parameters`.
  This allows activity models for different datasets or seasons to share one
  stellar rotation period.

**use_stellar_activity_decay**
* accepted values: `True` | **`False`**
* if `True`, replaces `Pdec` with `activity_decay` from `star_parameters`.

The flags can be declared in the `activity` common object or in a model that
supports them. Include `star_parameters` in that model's `common` list when
using either flag. See also the [simultaneous photometry and spectroscopy GP
example](../advanced_use/gp_rv_lc_fit.md).

## Examples

The [K2-141 RV and activity-indicator example](https://github.com/LucaMalavolta/PyORBIT_examples/blob/main/gaussian_processes/RV_GPtrained_2p.yaml)
uses one `activity` object for a quasi-periodic GP fitted to RV, BIS, and
S-index data. Its relevant configuration is:

```yaml
common:
  activity:
    boundaries:
      Prot: [10.0, 20.0]
      Pdec: [20.0, 1000.0]
      Oamp: [0.001, 1.0]
models:
  gp_quasiperiodic:
    model: gp_quasiperiodic
    common: activity
    hyperparameters_condition: True
    rotation_decay_condition: True
    boundaries:
      Hamp: [0.0, 100.0]
```

The [TOI-1807 two-season example](https://github.com/LucaMalavolta/PyORBIT_examples/blob/main/gp_multiple_seasons/RV_GPtrained_2seasons.yaml)
defines a separate activity object for each season, while both GP models use
the same stellar rotation period:

```yaml
common:
  activity_s01:
    model: activity
    boundaries:
      Pdec: [10.0, 100.0]
      Oamp: [0.001, 1.0]
  activity_s02:
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
    common:
      - activity_s01
      - star_parameters
    use_stellar_rotation_period: True
  gp_quasiperiodic_s02:
    model: gp_quasiperiodic
    common:
      - activity_s02
      - star_parameters
    use_stellar_rotation_period: True
```

The [two-component quasi-periodic example](https://github.com/LucaMalavolta/PyORBIT_examples/blob/main/gp_newkernels/model05_RV_tinygp_trainedQPSE.yaml)
also sets `Pcyc` in `activity` and gives `Hamp` and `Camp` separate bounds for
each dataset.
