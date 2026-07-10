(quasiperiodic_derivative_kernel)=

# Quasi-periodic plus derivative kernel


The kernel employed in the `gp_quasiperiodic_derivative` model is a combination of the quasi-periodic kernel and its first derivative, not different from using the [multivariate approach](../multidimensional_gps.md)on a single dataset. If we consider the correlated noise of a given dataset $\mathrm{D}$ as the combination of a Gaussian process and its first derivative:

```{math}
:label: gp_framework_original_onedataset

\Delta \mathrm{D} & = H_\mathrm{amp} G(t) + C_\mathrm{amp} G^\prime (t) \\

```

Then the final covariance between observations of G and its derivative at times $t_i$ and $t_j$ is given by:

```{math}
:label: quasiperiodic_derivative_pyorbit

\gamma (t_i, t_j)_{\rm full}  = H_\mathrm{amp}^2 \gamma_{\rm GP}^{(G,G)} (t_i, t_j)  + C_\mathrm{amp}^2 \gamma_{\rm GP}^{(dG,dG)} (t_i, t_j) (t)
```

Where $\gamma_{\rm GP}^{(G,G)} (t_i, t_j)$ 
is the [quasi-periodic kernel](quasiperiodic_kernel) and  $ \gamma^{(dG,dG)}(t_i, t_j) = \left.\left.\frac{\partial}{\partial t'} \frac{\partial}{\partial t} \gamma^{(G,G)}(t, t') \right|_{t=t_i} \right|_{t'=t_j} $. 
Whit a single dataset, the terms  $ \gamma^{(G,dG)}(t_i, t_j)$ and  $ \gamma^{(dG,G)}(t_i, t_j)$ cancel each other. See Section 3.3 of [Rajpaul et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.452.2269R/abstract) for more details.

This model can be useful when a single dataset needs a more flexible stellar-activity covariance than the standard quasi-periodic kernel, without switching to a multidimensional GP where the first derivative is usually employed.


## Model definition and requirements

The fastest implementation relies on `tinyGP`, but it requires a few extra tricks in the configuration file and execution (see [Caveats on the use of `tinyGP`](../running_pyorbit/tinygp_caveats)). If you use this model, cite the [Zenodo repository](https://zenodo.org/records/19035246).

An independent implementation that relies only on basic packages is maintained for legacy reasons; however, it is much slower, and I don't recommend using it.  

**model name**: `tinygp_quasiperiodic_derivative` 
- **available since version 11.2.6**
- required common object: `activity`
- implemented using  `tinygp` (version 0.3.0, [link to documentation](https://tinygp.readthedocs.io/en/stable/))
- GPU acceleration supported (instruction incoming)
- Read [Caveats on the use of `tinyGP`](../running_pyorbit/tinygp_caveats) carefully

**model name**: `gp_quasiperiodic_derivative`
- required common object: `activity`
- *direct* implementation relying only on `numpy` and `scipy`
- independent covariance matrix for each dataset


## Model parameters

| Name | Parameter | Common? | Definition | Notes |
| :--- | :-------- | :------ | :--------- | :---- |
| `Prot` | Rotational period of the star | common | `activity` | Replaced by `rotation_period` when `use_stellar_rotation_period: True` |
| `Pdec` | Decay timescale of active regions | common | `activity` | Replaced by `activity_decay` when `use_stellar_activity_decay: True` |
| `Oamp` | Coherence scale | common | `activity` | |
| `Hamp` | Amplitude of the quasi-periodic component | dataset | `activity` | |
| `Camp` | Amplitude of the derivative component | dataset | `activity` | |

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
  RVdata:
    file: datasets/star_RV_PyORBIT.dat
    kind: RV
    models:
      - radial_velocities
      - gp_qp_derivative
common:
  activity:
    boundaries:
      Prot: [10.0, 20.0]
      Pdec: [20.0, 1000.0]
      Oamp: [0.001, 1.0]
models:
  gp_qp_derivative:
    model: gp_quasiperiodic_derivative
    common: activity
    hyperparameters_condition: True
    rotation_decay_condition: True
    boundaries:
      Hamp: [0.0, 100.0]
      Camp: [0.0, 100.0]
```


