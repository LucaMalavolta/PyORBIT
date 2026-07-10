(sho_matern32_kernels)=

# SHO and Matern-3/2 kernels


The Stochastically-driven Harmonic Oscillator (SHO) and Matern-3/2 kernels are generic GP kernels commonly used to model correlated noise regardless of the specific origin. By default, their parameters are dataset-specific, but it is possible to change this behaviour using the available keywords.

```{tip}
For their general purpose, these kernels will not inherit the hyperparameters from the `activity` common model. However, it is possible to relate the hyperparameters to the star's rotation period and activity-decay timescale. 
```

The kernels have been implemented using [`celerite2`](https://celerite2.readthedocs.io/en/latest/), [`tinygp`](https://tinygp.readthedocs.io/en/stable/), and [`S+LEAF`](https://gitlab.unige.ch/delisle/spleaf). 
- If you use `tinygp`, cite the [Zenodo repository](https://zenodo.org/records/19035246) and check the [Caveats on the use of `tinyGP`](../running_pyorbit/tinygp_caveats)
- If you use `celerite2`, please cite [Foreman-Mackey et al. 2017](https://ui.adsabs.harvard.edu/abs/2017AJ....154..220F/abstract) and [Foreman-Mackey 2018](https://ui.adsabs.harvard.edu/abs/2018RNAAS...2...31F/abstract)
- If you use `S+LEAF` kernels, please cite [Delisle et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...638A..95D/abstract) and [Delisle et al. 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...659A.182D/abstract).


### SHO Kernel

This kernel has been introduced by [Foreman-Mackey et al. 2017](https://ui.adsabs.harvard.edu/abs/2017AJ....154..220F/abstract) with the following power spectral density:


```{math}
:label: quasiperiodic_squaredexponential_pyorbit

        S(\omega) = \sqrt{\frac{2}{\pi}} \frac{S_0\,\omega_0^4}
        {(\omega^2-{\omega_0}^2)^2 + {\omega_0}^2\,\omega^2/Q^2}
```

Where the hyperparameters represent:

- $\omega_0$: The undamped angular frequency
- $Q$: the quality factor
- $S0$: scale factor so that the power of the PSD is $S(0) =S_0\sqrt{2/\pi}$ 

In `PyORBIT`, the alternative parametrisation implemented in `celerite2` has been extended to all the `SHO` models:

- $\rho$: the undamped period of the oscillator, defined as $\rho = 2\,\pi / \omega_0$,
- $\tau$; the damping timescale of the process, defined as $\tau = 2\,Q / \omega_0$
- $sigma$_ the standard deviation of the process, defined as $\sigma = \sqrt{S_0\,\omega_0\,Q} $

If you use the SHO kernel, please cite [Foreman-Mackey et al. 2017](https://ui.adsabs.harvard.edu/abs/2017AJ....154..220F/abstract) 


## Model definition and requirements

**model name**: `tinygp_sho`
- required common object: `activity`
- implemented with `tinygp`

**model name**: `celerite2_sho`
- required common object: `activity`
- implemented with `celerite2`

**model name**: `spleaf_sho`
- required common object: `activity`
- implemented with `S+LEAF`

**model name**: `tinygp_matern32`
- required common object: `activity`
- implemented with `tinygp`

**model name**: `celerite2_matern32`
- required common object: `activity`
- implemented with `celerite2`

## Model parameters

### SHO models

| Name | Parameter | Common? | Definition | Notes |
| :--- | :-------- | :------ | :--------- | :---- |
| `sho_scale` | SHO undamped period, also known as `rho` | dataset | `activity` | Common when `use_shared_scale: True` or `use_shared_hyperparameters: True` |
| `sho_decay` | SHO damping timescale, also known as  `tau` | dataset | `activity` | Common when `use_shared_decay: True` or `use_shared_hyperparameters: True` |
| `sho_sigma` | Standard deviation of the process | dataset | `activity` | Common when `use_shared_hyperparameters: True` |

### Matern-3/2 models

| Name | Parameter | Common? | Definition | Notes |
| :--- | :-------- | :------ | :--------- | :---- |
| `matern32_scale` | Matern-3/2 scale, also known as  `rho` | dataset | `activity` | Common when `use_shared_scale: True` or `use_shared_hyperparameters: True` |
| `matern32_sigma` | Standard deviation of the process | dataset | `activity` | Common when `use_shared_hyperparameters: True` |

```{tip}
Some older result files may contain `matern32_rho`. The current models map it internally to `matern32_scale`.
```


## Keywords

Model-wide keywords, with the default value in boldface.

**use_shared_hyperparameters**
* accepted values: `True` | **`False`**
* moves all kernel parameters from dataset-specific to common parameters. For SHO models this includes `sho_scale`, `sho_decay`, and `sho_sigma`; for Matern-3/2 models this includes `matern32_scale` and `matern32_sigma`.

**use_shared_scale**
* accepted values: `True` | **`False`**
* makes only the kernel timescale common. For SHO models this is `sho_scale`; for Matern-3/2 models this is `matern32_scale`.

**use_shared_decay**
* accepted values: `True` | **`False`**
* makes only the SHO damping timescale `sho_decay` common.

**use_activity_Prot**
* accepted values: `True` | **`False`**
* replaces the kernel scale with `Prot` from the `activity` common object.

**use_activity_Pdec**
* accepted values: `True` | **`False`**
* replaces the kernel decay or scale with `Pdec` from the `activity` common object.

**use_stellar_rotation_period**
* accepted values: `True` | **`False`**
* replaces the kernel scale with `rotation_period` from `star_parameters`.

**use_stellar_activity_decay**
* accepted values: `True` | **`False`**
* replaces the kernel decay or scale with `activity_decay` from `star_parameters`.

For Matern-3/2 models, the replacement flags are mutually exclusive because there is only one timescale.

The SHO models also support some of the hyperparameter-condition keywords involving `Prot` and `Pdec` described in the [quasi-periodic kernel](quasiperiodic_kernel.md). In that case, the SHO scale is interpreted as `Prot`, and the SHO decay as `Pdec`.

## Examples

Dataset-specific SHO kernel:

```yaml
models:
  celerite2_sho:
    model: celerite2_sho
    boundaries:
      sho_scale: [0.1, 100.0]
      sho_decay: [0.1, 1000.0]
      sho_sigma: [0.0, 1.0]
```

Matern-3/2 kernel with a shared timescale:

```yaml
common:
  activity:
    boundaries:
      matern32_scale: [0.1, 100.0]
models:
  tinygp_matern32:
    model: tinygp_matern32
    common: activity
    use_shared_scale: True
    boundaries:
      matern32_sigma: [0.0, 1.0]
```

SHO kernel using stellar rotation and activity decay from `star_parameters`:

```yaml
common:
  activity:
    model: activity
  star:
    star_parameters:
      boundaries:
        rotation_period: [10.0, 20.0]
        activity_decay: [20.0, 1000.0]
models:
  spleaf_sho:
    model: spleaf_sho
    common:
      - activity
      - star_parameters
    use_stellar_rotation_period: True
    use_stellar_activity_decay: True
    boundaries:
      sho_sigma: [0.0, 1.0]
```


