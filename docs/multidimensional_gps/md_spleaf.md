(md_spleaf)=

# Exponential-sine periodic kernel (S+LEAF)

The Exponential-sine periodic (ESP) kernel is an approximation of the quasi-periodic (QP) kernel implemented in the `S+LEAF` software.
In the [`S+LEAF documentation`](https://obswww.unige.ch/~delisle/spleaf/doc/_autosummary/spleaf.term.ESPKernel.html#spleaf.term.ESPKernel), the QP kernel is referred to as the squared-exponential periodic (SEP) kernel. Despite the different name, the QP periodic kernel definition is the same as the one implemented in `PyORBIT` in either the [trained](../gaussian_process/quasiperiodic_kernel) or the [multidimensional](./md_quasiperiodic) approach. 




The quasi-periodic kernel is the preferred choice for the multidimensional GP. We follow the expression given by [Rajpaul et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.452.2269R/abstract):

```{math}
:label: quasiperiodic_grunblatt

\gamma ^{(G,G)}_{i,j} = \exp{ \left \{-\frac{\sin^2{[\pi(t_i - t_j)/ P_\mathrm{rot}]}}{2 O_\mathrm{amp} ^2} - \frac{(t_i-t_j)^2}{2 P_\mathrm{dec}^2} \right \} }


```{note}
If you use the multidimensional GP through `S+LEAF`, please don't forget to cite [Delisle et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...638A..95D/abstract) and [Delisle et al. 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...659A.182D/abstract)
```


### Implementation of the `S+LEAF` package

As `PyORBIT` naturally supports an unlimited number of datasets with heterogeneous cadences, the inclusion of the `S+LEAF` package has been quite straightforward. 
The GP hyperparameters preserve the same name as in the other kernels: the rotation period of the star $P_\mathrm{rot}$, the coherence scale $O_\mathrm{amp}$, and the decay time scale of the active regions $P_\mathrm{dec}$, corresponds to $P$, $\eta$, $\rho$ in `S+LEAF` documentation.

The coefficients $\alpha$ and $\beta$ introduced in the [`S+LEAF` multiGP [example](https://obswww.unige.ch/~delisle/spleaf/doc/multi.html)  are consequently renamed in *con_amp* and *rot_amp*, with the reference dataset easily identifiable from the terminal output.

## Model definition and requirements

**model name**: ``spleaf_multidimensional_esp``
- required common objects: ``activity``

## Keywords

Model-wide keywords, with the default value in bold face.

**n_harmonics**
* accepted values: integer | **4**
* Number of harmonics to include in the ESP approximation of the QP kernel

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

This example is nearly identical to the one presented for the [multidimensional QP kernel](md_quasiperiodic). The three main differences are:
- The model `spleaf_multidimensional_esp` is replacing `tinygp_multidimensional_quasiperiodic`
- the additional keyword `n_harmonics` is included in the example
- the `safe_reload` keyword in the `parameters` section is not longer required

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
    model: spleaf_multidimensional_esp
    common: activity
    n_harmonics: 4
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
