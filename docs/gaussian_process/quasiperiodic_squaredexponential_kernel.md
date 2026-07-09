(quasiperiodic_squaredexponential_kernel)=

# Quasi-periodic plus squared-exponential kernel


The quasi-periodic plus squared-exponential kernel combines the usual stellar-rotation term with an additional  squared-exponential component to model long-term activity cycles. This kernel has been introduced in [Basilicata et al. 2024](https://ui.adsabs.harvard.edu/abs/2024A%26A...686A.127B/abstract)

The full kernel $\gamma$ is given by the combination of a quasi-period kernel $\gamma_{\rm QP}$ and a squared-exponential (also knownas exponential squared) kernel $\gamma{\rm SE}$ 

```{math}
:label: quasiperiodic_squaredexponential_pyorbit

\gamma(t_i,t_j) =
H_\mathrm{amp}^2 \gamma_\mathrm{QP}(t_i,t_j; P_\mathrm{rot}, P_\mathrm{dec}, O_\mathrm{amp}) +
C_\mathrm{amp}^2 \gamma_\mathrm{SE}(t_i,t_j; P_\mathrm{cyc})
```

where the  squared-exponential kernel is equal to:

```{math}
:label: squaredexponential_noamp_pyorbit

\gamma_\mathrm{SE} (t_i, t_j) = \exp{ - \frac{(t_i-t_j)^2}{2 P_\mathrm{cyc}^2} \right \}
```

with $ P_\mathrm{cyc}$  being the correlation decay timescale of the activity cycle.

## Model definition and requirements

**model name**: `tinygp_quasiperiodicsquaredexponential`
- required common object: `activity`
- implemented with `tinygp`
- read [Caveats on the use of `tinyGP`](../running_pyorbit/tinygp_caveats) carefully

**model alias**
- `tinygp_quasiperiodic_squaredexponential`

There is no direct implementation of this kernel available.

## Keywords

Model-wide keywords, with the default value in boldface.

**hyperparameters_condition**
* accepted values: `True` | **`False`**
* activates the quasi-periodic hyperparameter condition described in the [quasi-periodic kernel](quasiperiodic_kernel).

**rotation_decay_condition**
* accepted values: `True` | **`False`**
* if activated, requires `Pdec > 2 Prot`.

**halfrotation_decay_condition**
* accepted values: `True` | **`False`**
* if activated, requires `Pdec > 0.5 Prot`.

**decay_rotation_factor** or **rotation_decay_factor**
* accepted values: float | **not used**
* if provided, requires `Pdec` to be larger than the specified factor times `Prot`.

**use_stellar_rotation_period**
* accepted values: `True` | **`False`**
* replaces `Prot` with `rotation_period` from `star_parameters`.

**use_stellar_activity_decay**
* accepted values: `True` | **`False`**
* replaces `Pdec` with `activity_decay` from `star_parameters`.

## Example

```yaml
inputs:
  LCdata:
    file: datasets/lightcurve_PyORBIT.dat
    kind: Phot
    models:
      - gp_qp_cycle
common:
  activity:
    boundaries:
      Prot: [10.0, 20.0]
      Pdec: [20.0, 1000.0]
      Pcyc: [100.0, 5000.0]
      Oamp: [0.001, 1.0]
models:
  gp_qp_cycle:
    model: tinygp_quasiperiodicsquaredexponential
    common: activity
    hyperparameters_condition: True
    rotation_decay_condition: True
    boundaries:
      Hamp: [0.0, 1.0]
      Camp: [0.0, 1.0]
```

## Model parameters

| Name | Parameter | Common? | Definition | Notes |
| :--- | :-------- | :------ | :--------- | :---- |
| `Prot` | Rotational period of the star | common | `activity` | Replaced by `rotation_period` when `use_stellar_rotation_period: True` |
| `Pdec` | Decay timescale of active regions | common | `activity` | Replaced by `activity_decay` when `use_stellar_activity_decay: True` |
| `Pcyc` | Timescale of the squared-exponential component | common | `activity` | |
| `Oamp` | Coherence scale | common | `activity` | |
| `Hamp` | Amplitude of the quasi-periodic component | dataset | `activity` | |
| `Camp` | Amplitude of the squared-exponential component | dataset | `activity` | |
