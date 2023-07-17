(multidimensional_gaussian_processes)=

# Multidimensional GPs

In the original Gaussian process framework ([Rajpaul et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.452.2269R/abstract), [Rajpaul et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.1847R/abstract)) the radial velocity datasets ($\Delta \mathrm{RV}$, after removing the deterministic part) and two activity indicators (in this example, $\mathrm{BIS}$ and  $\log{R^{\prime}_\mathrm{HK}}$) are modelled as a liner combination of an underlying Gaussian proccess $G(t)$  and its first derivative $G^\prime (t)$.

```{math}
:label: gp_framework_original

\Delta \mathrm{RV} & = V_c G(t) + V_r G^\prime (t) \\
\mathrm{BIS} & = B_c G(t) + B_r G^\prime (t) \\
\log{R^{\prime}_\mathrm{HK}} & = L_c G(t) \\
```

`PyORBIT` implementation produces results that are perfectly consistent with the GP framework by [Rajpaul et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.452.2269R/abstract)  and the multidimensional GP by [Barragán et al. 2022](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509..866B/abstract), as shown in [Nardiello et al. 2022, appendix D](https://ui.adsabs.harvard.edu/abs/2022A%26A...664A.163N/abstract). Note that the signs of the coefficients may vary according to the employed defintion of  $\Delta t = (t_i-t_j)$. 

```{admonition} Give credits to the authors
If you use the multidimensional GP, remeber to cite [Rajpaul et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.452.2269R/abstract)  and [Barragán et al. 2022](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509..866B/abstract).
For more details on the `PyORBIT` implementation of  multidimensional GP, please refer to [Nardiello et al. 2022, appendix D](https://ui.adsabs.harvard.edu/abs/2022A%26A...664A.163N/abstract)

```


### Kernels

The only implemented kernel for the underlying GP is the quasi-periodic one, following the expression given by [Grunblatt et al. 2015](https://ui.adsabs.harvard.edu/abs/2015ApJ...808..127G/abstract):

```{math}
:label: quasiperiodic_grunblatt

\gamma ^{(G,G)}_{i,j} = \exp{ \left \{-\frac{\sin^2{[\pi(t_i - t_j)/\theta]}}{2 w ^2} - \left (\frac{t_i-t_j}{\lambda} \right )^2 \right \} }
```

where $\theta$ is equivalent to the rotation period of the star, $w$ is the coherence scale, and $\lambda$ is usually associated to the decay time scale fo the active regions

```{important}
It is common to have a factor 2 in the denominator of the aperiodic variation (i.e., $2 \lambda$ rather than $\lambda$) in {eq}`quasiperiodic_grunblatt`. In such a case, it is sufficient to multiply the value of $\lambda$ of `PyORBIT` by a factor $\sqrt(2)$ - keep it in mind when assigning priors!
```


### Improvements

In `PyORBIT` the framework implementation has been expanded to support an unlimited number of datasets. Additionally, it is not required for the datasets to have the same dimensions, thus allowing the use of several RV datasets even if activity indicators are not available for all of them.

The coefficients $X_c$ and $X_c$ are thus renamed in *con_amp* and *rot_amp*, with the reference dataset easily identifiable from the terminal output.

## Model definition and requirements

- model name: ``gp_multidimensional_quasiperiodic``
- required common objects : ``activity``

## Keywords

Model-wide keywords, with the default value in bold face.

**hyperparameters_condition**
* accepted values: `True` | **`False`**
* activate the conditions $ \lambda ^ 2 > (3/4 \pi) \theta ^2 w ^ 2 $ (adapted from [Rajpaiul 2017](https://ui.adsabs.harvard.edu/abs/2017PhDT.......229R/abstract) and [Rajpaul et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.1847R/abstract) to take into account the factor 2 in the denominator of the aperiodic variation) to ensure that the QP function has at least one non-trivial turning point.

**rotation_decay_condition**
* accepted values: `True` | **`False`**
* if activated, it ensures that the decay time scale of the activity regions $\lambda$ is at least twice the rotational period of the star $\theta$

**derivative**
* accepted values: list of dataset using the model
* if not provided, default is **`True`** for all datasets except `H-alpha` `S_index` `Ca_HK`, and `FWHM`. If provided, default is **`False`** for all the dataset not explicitly mentioned.
* If **`True`**, the first derivative of the ubderlying Gaussian process is included, otherwise `rot_amp`F is fixed to zero.

## Examples

In the following example, ore RV dataset comes together with two activity indicators, the `BIS` and the `FWHM` of the CCF, the latter as replacement of  $\log{R^{\prime}_\mathrm{HK}}$.

Here `gp_multidimensional` is the label that we assign to the `gp_multidimensional_quasiperiodic` model. The label is assigned to each dataset, including the RV dataset for which also the RV model is included.

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
    model: gp_multidimensional_quasiperiodic
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
Since the coefficients are always copuled, there is a degenracy between a given set and the opposite one (i.e., all the coefficient with opposite sign). This degeneracy can be solved by force one coefficient to be positive, e.g., `rot_amp` for the `RV_data` in the example.
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
| rot_amp    | coefficient of $G(t)$ component | dataset | ``activity``     | |
| con_amp    | coefficient of $G^\prime (t)$ component | dataset | ``activity``     | |
