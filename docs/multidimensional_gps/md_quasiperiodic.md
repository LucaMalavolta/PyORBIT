(md_quasiperiodic)=

# Quasi-periodic kernel

The quasi-periodic kernel is the preferred choice for the multidimensional GP. We follow the expression given by [Rajpaul et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.452.2269R/abstract):

```{math}
:label: quasiperiodic_grunblatt

\gamma ^{(G,G)}_{i,j} = \exp{ \left \{-\frac{\sin^2{[\pi(t_i - t_j)/ P_\mathrm{rot}]}}{2 O_\mathrm{amp} ^2} - \frac{(t_i-t_j)^2}{2 P_\mathrm{dec}^2} \right \} }
```

Where $P_\mathrm{rot}$ is equivalent to the rotation period of the star, $O_\mathrm{amp}$ is the coherence scale, and $P_\mathrm{dec}$ is usually associated with the decay time scale of the active regions.

```{note}
Check the documentation page on the [quasi-periodic kernel](../gaussian_process/quasiperiodic_kernel) for additional information on this kernel
```



### Improvements

In `PyORBIT` the framework implementation has been expanded to support an unlimited number of datasets. Additionally, it is not required for the datasets to have the same dimensions, thus allowing the use of several RV datasets even if activity indicators are not available for all of them.

The coefficients $X_c$ and $X_c$ are thus renamed in *con_amp* and *rot_amp*, with the reference dataset easily identifiable from the terminal output.

## Model definition and requirements

**model name**: ``tinygp_multidimensional_quasiperiodic``
- required common objects: ``activity``
- GPU acceleration supported (instruction incoming)
- Read [Caveats on the use of `tinyGP`](../running_pyorbit/tinygp_caveats) carefully

**model name**: ``gp_multidimensional_quasiperiodic``
- required common objects: ``activity``
- *direct* implementation using `numpy` and `scipy` packages.

**model name**: ``gp_multidimensional_quasiperiodic_numba``
- required common objects: ``activity``
- *direct* implementation using `numpy` and `scipy` packages with `numba` acceleration.


## Keywords

Model-wide keywords, with the default value in bold face.

**hyperparameters_condition**
* accepted values: `True` | **`False`**
* activate the conditions $ P_\mathrm{rot}  ^ 2 > \frac{3}{2 \pi} P_\mathrm{rot} ^2 O_\mathrm{amp} ^ 2 $ from [Rajpaiul 2017](https://ui.adsabs.harvard.edu/abs/2017PhDT.......229R/abstract) and [Rajpaul et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.1847R/abstract), to ensure that the QP function has at least one non-trivial turning point.

**rotation_decay_condition**
* accepted values: `True` | **`False`**
* if activated, it ensures that the decay time scale of the activity regions $\lambda$ is at least twice the rotational period of the star $\theta$

**derivative**
* accepted values: list of datasets using the model
* if not provided, default is **`True`** for all datasets except `H-alpha` `S_index` `Ca_HK`, and `FWHM`. If provided, default is **`False`** for all the datasets not explicitly mentioned.
* If **`True`**, the first derivative of the underlying Gaussian process is included, otherwise, `rot_amp`F is fixed to zero.

## Examples

In the following example, one RV dataset comes together with two activity indicators, the `BIS` and the `FWHM` of the CCF, the latter as a replacement of  $\log{R^{\prime}_\mathrm{HK}}$.

Here `gp_multidimensional` is the label that we assign to the `tinygp_multidimensional_quasiperiodic` model. The label is assigned to each dataset, including the RV dataset for which also the RV model is included.

The `common:planets` section includes three planets, of which one is transiting - with Gaussian priors on the period and time of inferior conjunction -  and two non-transiting planets with a prior on the eccentricity.
The `common:activity` section provides the hyperparameters for the GP shared among all the datasets. The example shows how to assign boundaries and priors to the parameters.
The model keywords and the boundaries for the dataset-specific parameters are listed in `models:gp_multidimensional`

```yaml
inputs:
  RVdata:
    file: RVS_PyORBIT.dat
    kind: RV
    models:
      - radial_velocities
      - gp_multidimensional
  BISdata:
    file: BIS_PyORBIT.dat
    kind: BIS
    models:
      - gp_multidimensional
  FWHMdata:
    file: FWHM_PyORBIT.dat
    kind: FWHM
    models:
      - gp_multidimensional
common:
  planets:
    b:
      orbit: circular
      use_time_inferior_conjunction: True
      boundaries:
        P: [2.21000, 2.240000]
        K: [0.001, 20.0]
        Tc: [59144.60, 59144.63]
      priors:
        P: ['Gaussian', 2.2241951, 0.00000030]
        Tc: ['Gaussian', 59144.616171, 0.000284]
    c:
      orbit: keplerian
      boundaries:
        P: [2.0, 100.0]
        K: [0.001, 30.0]
        e: [0.00, 0.70]
      priors:
        e: ['Gaussian', 0.00, 0.098]
    d:
      orbit: keplerian
      boundaries:
        P: [2.0, 100.0]
        K: [0.001, 30.0]
        e: [0.00, 0.70]
      priors:
        e: ['Gaussian', 0.00, 0.098]
  activity:
    boundaries:
      Prot: [20.0, 30.0]
      Pdec: [30.0, 1000.0]
      Oamp: [0.01, 1.0]
    priors:
      Prot: ['Gaussian', 14.00, 0.50]
      Oamp: ['Gaussian', 0.35, 0.035]
  star:
    star_parameters:
      priors:
        mass: ['Gaussian', 0.65, 0.05]
        radius: ['Gaussian', 0.624, 0.005]
        density: ['Gaussian', 2.65, 0.08]
models:
  radial_velocities:
    planets:
      - b
      - c
      - d
  gp_multidimensional:
    model: tinygp_multidimensional_quasiperiodic
    common: activity
    hyperparameters_condition: True
    rotation_decay_condition: True
    RVdata:
      boundaries:
        rot_amp: [0.0, 20.0] #at least one must be positive definite
        con_amp: [-20.0, 20.0]
      derivative: True
    BISdata:
      boundaries:
        rot_amp: [-20.0, 20.0]
        con_amp: [-20.0, 20.0]
      derivative: True
    FWHMdata:
      boundaries:
        con_amp: [-50., 50.]
      derivative: False
parameters:
  Tref: 59200.00
  safe_reload: True
  low_ram_plot: True
  plot_split_threshold: 1000
  cpu_threads: 16
solver:
  pyde:
    ngen: 50000
    npop_mult: 4
  emcee:
    npop_mult: 4
    nsteps: 100000
    nburn: 25000
    nsave: 25000
    thin: 100
    #use_threading_pool: False
  nested_sampling:
    nlive: 1000
    sampling_efficiency: 0.30
  recenter_bounds: True

```

```{tip}
Since the coefficients are always coupled, there is a degeneracy between a given set and the opposite one (i.e., all the coefficient with opposite sign). This degeneracy can be solved by force one coefficient to be positive, e.g., `rot_amp` for the `RV_data` in the example.
```

A simpler way to prepare the configuration file if you are not specifying the boundaries is to list the datasets with derivatives under the same keyword:

```yaml
...
  gp_multidimensional:
    model: gp_multidimensional_quasiperiodic
    common: activity
    hyperparameters_condition: True
    rotation_decay_condition: True
    derivative:
      RVdata: True
      BISdata: True
      FWHMdata: False
parameters:
...

```


## Model parameters

The following parameters will be inherited from the common model (column *Common?: common*) or a different value will be assigned for each dataset (column *Common?: dataset*)

| Name        | Parameter | Common?  | Definition  | Notes |
| :---        | :-------- | :-------------  | :-----  | :---- |
| Prot      | rotational period of the star $\theta$ | common | ``activity``     | |
| Pdec      | Decay time scale of active regions $\lambda$ | common | ``activity``     | |
| Oamp | Coherence scale $w$ | common | ``activity`` |   |
| con_amp    | coefficient of $G(t)$ component | dataset | ``activity``     | |
| rot_amp    | coefficient of $G^\prime (t)$ component | dataset | ``activity``     | |
