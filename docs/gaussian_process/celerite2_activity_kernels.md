(celerite2_activity_kernels)=

# Granulation, Oscillation, Rotation kernels

This model provides multi-component Gaussian Process regression to account for stellar oscillation, granulation, annd rotational modulation, as proposed by [Barros et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...634A..75B/abstract), with some modifications as detailed in Section 5.1 of [Nardiello et al. 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...664A.163N/abstract). The former also provides a practical example of its application.

The kernel comprises:




If you use this model, please cite the `celerite2` [Foreman-Mackey et al. 2017](https://ui.adsabs.harvard.edu/abs/2017AJ....154..220F/abstract) and [Foreman-Mackey 2018](https://ui.adsabs.harvard.edu/abs/2018RNAAS...2...31F/abstract)


## Model definition and requirements

**model name**: `celerite2_rotation`
- required common object: `activity`
- rotation term only

**model name**: `celerite2_granulation_rotation`
- required common object: `activity`
- one granulation SHO term plus one rotation term

**model name**: `celerite2_granulation_oscillation_rotation`
- required common object: `activity`
- configurable number of granulation and oscillation SHO terms plus one rotation term


## Model parameters

### `celerite2_rotation`

| Name | Parameter | Common? | Definition | Notes |
| :--- | :-------- | :------ | :--------- | :---- |
| `Prot` | Primary period of variability | common | `activity` | Replaced by `rotation_period` when `use_stellar_rotation_period: True` |
| `rot_Q0` | Base quality factor of the secondary oscillation | common | `activity` | |
| `rot_deltaQ` | Difference between the quality factors of the two modes | common | `activity` | |
| `rot_fmix` | Fractional amplitude of the secondary mode | common | `activity` | |
| `rot_sigma` | Standard deviation of the process | dataset | `activity` | Common when `use_shared_hyperparameters: True` |

### `celerite2_granulation_rotation`

| Name | Parameter | Common? | Definition | Notes |
| :--- | :-------- | :------ | :--------- | :---- |
| `grn_period` | Granulation SHO timescale | common | `activity` | |
| `Prot` | Primary period of variability | common | `activity` | Replaced by `rotation_period` when `use_stellar_rotation_period: True` |
| `rot_Q0` | Base quality factor of the secondary oscillation | common | `activity` | |
| `rot_deltaQ` | Difference between the quality factors of the two modes | common | `activity` | |
| `rot_fmix` | Fractional amplitude of the secondary mode | common | `activity` | |
| `grn_sigma` | Granulation amplitude | dataset | `activity` | |
| `rot_sigma` | Rotation amplitude | dataset | `activity` | |

### `celerite2_granulation_oscillation_rotation`

| Name | Parameter | Common? | Definition | Notes |
| :--- | :-------- | :------ | :--------- | :---- |
| `Prot` | Primary period of variability | common | `activity` | Replaced by `rotation_period` when `use_stellar_rotation_period: True` |
| `rot_Q0` | Base quality factor of the secondary oscillation | common | `activity` | |
| `rot_deltaQ` | Difference between the quality factors of the two modes | common | `activity` | |
| `rot_fmix` | Fractional amplitude of the secondary mode | common | `activity` | |
| `rot_sigma` | Rotation amplitude | dataset | `activity` | Common when `common_amplitudes: True` |
| `grn_k*_period` | Timescale of each granulation SHO term | common | `activity` | `*` is the kernel index, starting at 0 |
| `grn_k*_sigma` | Amplitude of each granulation SHO term | dataset | `activity` | Common when `common_amplitudes: True` |
| `osc_k*_period` | Timescale of each oscillation SHO term | common | `activity` | `*` is the kernel index, starting at 0 |
| `osc_k*_Q0` | Quality factor of each oscillation SHO term | common | `activity` | |
| `osc_k*_sigma` | Amplitude of each oscillation SHO term | dataset | `activity` | Common when `common_amplitudes: True` |



## Keywords

Model-wide keywords, with the default value in boldface.

**use_stellar_rotation_period**
* accepted values: `True` | **`False`**
* replaces `Prot` with `rotation_period` from `star_parameters`.

**use_shared_hyperparameters**
* accepted values: `True` | **`False`**
* for `celerite2_rotation`, also makes `rot_sigma` a common parameter.

**rotation_kernels**
* accepted values: integer | **`1`**
* used by `celerite2_granulation_oscillation_rotation`. Only one rotation kernel is currently supported.

**granulation_kernels**
* accepted values: integer | **`2`**
* used by `celerite2_granulation_oscillation_rotation`. Adds parameters named `grn_k0_period`, `grn_k0_sigma`, `grn_k1_period`, `grn_k1_sigma`, and so on.

**oscillation_kernels**
* accepted values: integer | **`1`**
* used by `celerite2_granulation_oscillation_rotation`. Adds parameters named `osc_k0_period`, `osc_k0_sigma`, `osc_k0_Q0`, and so on.

**common_amplitudes**
* accepted values: `True` | **`False`**
* used by `celerite2_granulation_oscillation_rotation`. If activated, the `*_sigma` amplitudes are common parameters instead of dataset-specific parameters.

## Examples

Rotation-only model:

```yaml
inputs:
  LCdata:
    file: datasets/lightcurve_PyORBIT.dat
    kind: Phot
    models:
      - celerite2_rotation
common:
  activity:
    boundaries:
      Prot: [10.0, 20.0]
      rot_Q0: [0.1, 100.0]
      rot_deltaQ: [0.1, 100.0]
      rot_fmix: [0.01, 1.0]
models:
  celerite2_rotation:
    model: celerite2_rotation
    common: activity
    boundaries:
      rot_sigma: [0.0, 1.0]
```

Granulation, oscillation, and rotation model:

```yaml
models:
  celerite2_activity:
    model: celerite2_granulation_oscillation_rotation
    common: activity
    granulation_kernels: 2
    oscillation_kernels: 1
    common_amplitudes: False
    boundaries:
      rot_sigma: [0.0, 1.0]
      grn_k0_sigma: [0.0, 1.0]
      grn_k1_sigma: [0.0, 1.0]
      osc_k0_sigma: [0.0, 1.0]
```
