(activity)=

# Stellar activity

The `activity` common object stores hyperparameters that are shared by stellar-activity models. A model can still decide that some of these quantities are dataset-specific; for example, the rotation period and decay timescale are usually shared, while the covariance amplitude is often fitted independently for each dataset.

```{note}
The `activity` common object only defines the available parameters, boundaries, priors, spaces, and fixed values. The model listed in the `models` section decides which parameters are used and whether they are common to all datasets or local to one dataset.
```

## Model definition

- common object name: `activity`
- default common object for most Gaussian-process models of stellar activity

## Shared stellar parameters

Several activity models can take the stellar rotation period or the activity decay timescale from the `star_parameters` common object instead of from `activity`.

**use_stellar_rotation_period**
* accepted values: `True` | **`False`**
* replaces `Prot` with `rotation_period` from the `star_parameters` common object.

**use_stellar_activity_decay**
* accepted values: `True` | **`False`**
* replaces `Pdec` with `activity_decay` from the `star_parameters` common object.

The same flags can be activated either in the `activity` common object or in the individual model definition. See the [simultaneous photometry and spectroscopy GP example](../advanced_use/gp_rv_lc_fit) for a complete configuration.

## Example

```yaml
common:
  activity:
    boundaries:
      Prot: [10.0, 20.0]
      Pdec: [20.0, 1000.0]
      Oamp: [0.001, 1.0]
    priors:
      Prot: ['Gaussian', 14.0, 0.5]
  star:
    star_parameters:
      boundaries:
        rotation_period: [10.0, 20.0]
        activity_decay: [20.0, 1000.0]
```

## Parameters

Only the parameters required by the selected models have to be defined in the YAML file.

| Name | Parameter | Typical use |
| :--- | :-------- | :---------- |
| `Prot` | Stellar rotation period | Quasi-periodic, rotation, and multidimensional activity kernels |
| `Pdec` | Decay timescale of active regions | Quasi-periodic and exponential-sine periodic kernels |
| `Pcyc` | Long-timescale cycle or squared-exponential timescale | Quasi-periodic plus squared-exponential kernels |
| `Oamp` | Coherence scale of the periodic component | Quasi-periodic and exponential-sine periodic kernels |
| `Hamp` | Covariance amplitude | Dataset-specific amplitude in trained GP models |
| `Camp` | Secondary covariance amplitude | Cosine, derivative, or cycle component depending on the model |
| `rot_sigma` | Amplitude of a celerite2 rotation term | celerite2 rotation kernels |
| `rot_Q0` | Base quality factor of the celerite2 rotation term | celerite2 rotation kernels |
| `rot_deltaQ` | Difference between the rotation-mode quality factors | celerite2 rotation kernels |
| `rot_fmix` | Fractional amplitude of the secondary rotation mode | celerite2 rotation kernels |
| `grn_period` | Granulation SHO timescale | celerite2 granulation plus rotation |
| `grn_sigma` | Granulation SHO amplitude | celerite2 granulation plus rotation |
| `grn_k*_period` | Timescale of the `k`-th granulation SHO term | configurable celerite2 granulation/oscillation/rotation |
| `grn_k*_sigma` | Amplitude of the `k`-th granulation SHO term | configurable celerite2 granulation/oscillation/rotation |
| `osc_k*_period` | Timescale of the `k`-th oscillation SHO term | configurable celerite2 granulation/oscillation/rotation |
| `osc_k*_sigma` | Amplitude of the `k`-th oscillation SHO term | configurable celerite2 granulation/oscillation/rotation |
| `osc_k*_Q0` | Quality factor of the `k`-th oscillation SHO term | configurable celerite2 granulation/oscillation/rotation |
| `sho_scale` | SHO undamped period, or `rho` | SHO kernels |
| `sho_decay` | SHO damping timescale, or `tau` | SHO kernels |
| `sho_sigma` | SHO process amplitude | SHO kernels |
| `matern32_scale` | Matern-3/2 scale, also accepted as `matern32_rho` in some models | Matern-3/2 kernels |
| `matern32_sigma` | Matern-3/2 covariance amplitude | Single-output Matern-3/2 kernels |
| `rot_amp` | Coefficient of the derivative component of a latent GP | Multidimensional GP models |
| `con_amp` | Coefficient of the latent GP itself | Multidimensional GP models |
| `cos_amp` | Coefficient of the cosine latent component | Multidimensional quasi-periodic plus cosine kernel |
| `cos_der` | Coefficient of the derivative of the cosine component | Multidimensional quasi-periodic plus cosine kernel |
| `cyc_amp` | Coefficient of the squared-exponential cycle component | Multidimensional quasi-periodic plus squared-exponential kernels |
| `cyc_der` | Coefficient of the derivative of the squared-exponential cycle component | Multidimensional quasi-periodic plus squared-exponential kernel |
| `matern32_multigp_sigma` | Coefficient of the Matern-3/2 latent GP | Multidimensional Matern-3/2 kernel |
| `matern32_multigp_sigma_deriv` | Coefficient of the derivative of the Matern-3/2 latent GP | Multidimensional Matern-3/2 kernel |
| `Vc`, `Vr`, `Lc`, `Bc`, `Br` | Original Rajpaul-framework coefficients | Legacy `gp_framework_quasiperiodic` model |
