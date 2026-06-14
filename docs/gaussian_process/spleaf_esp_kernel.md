(spleaf_esp_kernel)=

# Exponential-sine periodic kernel (S+LEAF)

The `S+LEAF` exponential-sine periodic (ESP) kernel is a fast approximation of the quasi-periodic kernel. It is useful when one wants an independent covariance matrix for each dataset while sharing the stellar-activity hyperparameters.

The hyperparameters use the same names as the quasi-periodic models:

- `Prot`: stellar rotation period
- `Pdec`: decay timescale of active regions
- `Oamp`: coherence scale of the periodic component
- `Hamp`: dataset-specific covariance amplitude

```{note}
If you use this model, please cite `S+LEAF`: [Delisle et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...638A..95D/abstract) and [Delisle et al. 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...659A.182D/abstract).
```

## Model definition and requirements

**model name**: `spleaf_esp`
- required common object: `activity`
- implemented with the `S+LEAF` package
- default implementation for the one-dimensional ESP kernel

**model aliases**
- `spleaf_exponentialsineperiodic`
- `spleaf_esp_slow`
- `spleaf_exponentialsineperiodic_slow`

The `slow` implementation builds the `S+LEAF` covariance directly at each likelihood call. The default `spleaf_esp` implementation keeps a reusable covariance object and is normally preferred.

## Keywords

Model-wide keywords, with the default value in boldface.

**n_harmonics**
* accepted values: integer | **`4`**
* number of harmonics used by the ESP approximation.

**hyperparameters_condition**
* accepted values: `True` | **`False`**
* activates the quasi-periodic hyperparameter condition described in the [quasi-periodic kernel](quasiperiodic_kernel).

**rotation_decay_condition**
* accepted values: `True` | **`False`**
* if activated, requires `Pdec > 2 Prot`.

**use_stellar_rotation_period**
* accepted values: `True` | **`False`**
* replaces `Prot` with `rotation_period` from `star_parameters`.

**use_stellar_activity_decay**
* accepted values: `True` | **`False`**
* replaces `Pdec` with `activity_decay` from `star_parameters`.

## Example

```yaml
inputs:
  LC_K2:
    file: datasets/K2_lightcurve_PyORBIT.dat
    kind: Phot
    models:
      - spleaf_esp
common:
  activity:
    boundaries:
      Prot: [10.0, 20.0]
      Pdec: [20.0, 1000.0]
      Oamp: [0.001, 1.0]
models:
  spleaf_esp:
    model: spleaf_esp
    common: activity
    n_harmonics: 4
    hyperparameters_condition: True
    rotation_decay_condition: True
    boundaries:
      Hamp: [0.0, 1.0]
```

## Model parameters

The following parameters will be inherited from the common model (column *Common?: common*) or a different value will be assigned for each dataset (column *Common?: dataset*).

| Name | Parameter | Common? | Definition | Notes |
| :--- | :-------- | :------ | :--------- | :---- |
| `Prot` | Rotational period of the star | common | `activity` | Replaced by `rotation_period` when `use_stellar_rotation_period: True` |
| `Pdec` | Decay timescale of active regions | common | `activity` | Replaced by `activity_decay` when `use_stellar_activity_decay: True` |
| `Oamp` | Coherence scale | common | `activity` | |
| `Hamp` | Amplitude of the kernel | dataset | `activity` | |
