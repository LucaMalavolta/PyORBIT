(quasi_periodic_kernel)=

# Quasi-periodic kernel

The most common and most reliable kernel for stellar activity is the quasi-periodic one, following the expression given by [Rajpaul et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.452.2269R/abstract)

```{math}
:label: quasiperiodic_pyorbit

\gamma (t_i, t_j) = H_\mathrm{amp}^2  \exp{ \left \{-\frac{\sin^2{[\pi(t_i - t_j)/ P_\mathrm{rot}]}}{2 O_\mathrm{amp} ^2} - \frac{(t_i-t_j)^2}{2 P_\mathrm{dec}^2} \right \} }
```

where $P_\mathrm{rot}$ is equivalent to the rotation period of the star, $O_\mathrm{amp}$ is the coherence scale, and $P_\mathrm{dec}$ is usually associated with the decay time scale of the active regions.

```{important}
It is common to have a factor 2 in the denominator of the aperiodic variation or in the numerator of the sinusoidal term. In such a case, it is sufficient to multiply/divide the value value of `PyORBIT` by a factor $\sqrt(2)$ - keep it in mind when assigning priors!

A comparison between different formulations of the quasi-periodic kernel is provided in the section [Differences among various parametrization](differences-among-various-parametrizations) of this page.

```

## Model definition and requirements

The fastest implementation relies on `tinyGP`, but it requires a few extra tricks in the configuration file and execution (see [Caveats on the use of `tinyGP`](../running_pyorbit/tinygp_caveats) )
The original implementation based on `george` is still available. An independent implementation relying only on basic packages is available, however it is much slower.

**model name**: `tinygp_quasiperiodic`
- required common object: `activity`
- implemented using  `tinygp` (version 0.3.0, [link to documentation](https://tinygp.readthedocs.io/en/stable/))
- GPU acceleration supported (instruction incoming)
- Read [Caveats on the use of `tinyGP`](../running_pyorbit/tinygp_caveats) carefully

**model name**: `gp_quasiperiodic`
- required common object: `activity`
- implemented using  `george` (version 0.4.0, [Ambikasaram et al. 2015](https://ui.adsabs.harvard.edu/abs/2015ITPAM..38..252A), [link to documentation](https://george.readthedocs.io/))

'''warning
Starting with `PyORBIT 10.3`, only version >=0.3 of `tinygp` is supported. 
'''

**model name**: `gp_quasiperiodic_alternative`
- required common object: `activity`
- *direct* implementation relying only on `numpy` and `scipy`


## Model parameters

The following parameters will be inherited from the common model (column *Common?: common*) or a different value will be assigned for each dataset (column *Common?: dataset*)

| Name        | Parameter | Common?  | Definition  | Notes |
| :---        | :-------- | :-------------  | :-----  | :---- |
| Prot      | Rotational period of the star | common | ``activity``     | |
| Pdec      | Decay time scale of active regions | common | ``activity``     | |
| Oamp | Coherence scale | common | ``activity`` |   |
| Hamp  | Amplitude of the kernel | dataset | ``activity``     | |


## Keywords

Model-wide keywords, with the default value in boldface.

**hyperparameters_condition**
* accepted values: `True` | **`False`**
* activate the conditions $ P_\mathrm{rot}  ^ 2 > \frac{3}{2 \pi} P_\mathrm{rot} ^2 O_\mathrm{amp} ^ 2 $ from [Rajpaiul 2017](https://ui.adsabs.harvard.edu/abs/2017PhDT.......229R/abstract) and [Rajpaul et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.1847R/abstract), to ensure that the QP function has at least one non-trivial turning point.

**rotation_decay_condition**
* accepted values: `True` | **`False`**
* if activated, it ensures that the decay time scale of the activity regions $\lambda$ is at least twice the rotational period of the star $\theta$

**use_stellar_rotation_period**
* accepted values: `True` | **`False`**
* if activated, the parameter `Prot` from the `activity` *common model* will be replaced by the parameter `rotation_period` from the `star_parameters` *common model*. In this way, a unique parameter can be used by different models, e.g., stellar activity and Rossiter-McLaughlin modeling. It can also be useful if you want to use independent GP hyperparameters over several observational seasons while using a single parameter for the rotational period of the star.

## Examples

In the following example, a Radial Velocity (`RV`) dataset comes together with two activity indicators, the Bisector Inverse Span (`BIS`)  of the CCF and the Mount Wilson S index (`S_index`).

The `common:activity` section provides the hyperparameters for the GP shared among all the datasets. The example shows how to assign boundaries and priors to the parameters.
The model keywords and the boundaries for the dataset-specific parameters are listed in `models:gp_quasiperiodic`

```{code-block} yaml
:lineno-start: 1

inputs:
  RVdata:
    file: datasets/K2-141_RV_PyORBIT.dat
    kind: RV
    models:
      - radial_velocities
      - gp_quasiperiodic
  BISdata:
    file: datasets/K2-141_BIS_PyORBIT.dat
    kind: BIS
    models:
      - gp_quasiperiodic
  Sdata:
    file: datasets/K2-141_Sindex_PyORBIT.dat
    kind: S_index
    models:
      - gp_quasiperiodic
common:
  planets:
    b:
      orbit: circular
      use_time_inferior_conjunction: True
      boundaries:
        P: [0.2750, 0.2850]
        K: [0.001, 20.0]
        Tc: [57744.00, 57744.10]
      priors:
        P: ['Gaussian', 0.280324956, 0.000000067]
        Tc: ['Gaussian', 57744.071508, 0.000103]
      spaces:
        P: Linear
        K: Linear
    c:
      orbit: keplerian
      parametrization: Eastman2013
      use_time_inferior_conjunction: True
      boundaries:
        P: [7.70, 7.80]
        K: [0.001, 20.0]
        Tc: [58371.00, 58371.10]
        e: [0.00, 0.70]
      priors:
        e: ['Gaussian', 0.00, 0.098]
        P: ['Gaussian', 7.7489943, 0.0000149]
        Tc: ['Gaussian', 58371.07415, 0.000652]
      spaces:
        P: Linear
        K: Linear
  activity:
    boundaries:
      Prot: [10.0, 20.0]
      Pdec: [20.0, 1000.0]
      Oamp: [0.001, 1.0]
    #priors:
    #  Oamp: ['Gaussian', 0.35, 0.035]
  star:
    star_parameters:
      priors:
        mass: ['Gaussian', 0.708, 0.028]
        radius: ['Gaussian', 0.681, 0.018]
        density: ['Gaussian', 2.65, 0.08]
models:
  radial_velocities:
    planets:
      - b
      - c
  gp_quasiperiodic:
    model: gp_quasiperiodic
    common: activity
    hyperparameters_condition: True  # Condition from Rajpaul 2017, Rajpaul+2021
    rotation_decay_condition: True # It forces the decay timescale to be at least twice the rotational period
    boundaries:
      Hamp: [0.0, 100.0] # same range for all datasets
parameters:
  Tref: 59200.00
solver:
  pyde:
    ngen: 50000
    npop_mult: 4
  emcee:
    npop_mult: 4
    nsteps: 50000
    nburn: 20000
    nsave: 25000
    thin: 100
  recenter_bounds: True

```


In the example below, one dataset is divided into several observing seasons. The decay time scale of active regions, the coherence scale,  and the amplitude of the kernel are independent for each season, while the rotational period of the star is shared across all the datasets. \
Note the use of two `activity` *common models*  marked with the labels `activity_s01` and `activity_s02`. Correspondingly, two  `gp_quasiperiodic` *models* are employed, marked as `gp_quasiperiodic_s01` and `gp_quasiperiodic_s02`. The prior on the `rotation_period` is stored in the `star_aparameters` *common_model*. As a consequence, both the `activity` and `star_parameters` *common models* must be recalled in the `gp_quasiperiodic` *model*.

The use of a `common_offset` ensures the use of a single offset parameter for each type of dataset, avoiding degeneracies in the presence of long-period signals. \
The use of a `common_jitter` model may be motivated as well.


```{code-block} yaml
:lineno-start: 1
inputs:
  RVdata_s01:
    file: datasets/TOI1807_RV_s01_PyORBIT.dat
    kind: RV
    models:
      - radial_velocities
      - gp_quasiperiodic_s01
      - common_offset_RV
  BISdata_s01:
    file: datasets/TOI1807_BIS_s01_PyORBIT.dat
    kind: BIS
    models:
      - gp_quasiperiodic_s01
      - common_offset_BIS
  ...
  RVdata_s02:
    file: datasets/TOI1807_RV_s02_PyORBIT.dat
    kind: RV
    models:
      - radial_velocities
      - gp_quasiperiodic_s02
      - common_offset_RV
  BISdata_s02:
    file: datasets/TOI1807_BIS_s02_PyORBIT.dat
    kind: BIS
    models:
      - gp_quasiperiodic_s02
      - common_offset_BIS
    ...
common:
  planets:
    b:
      orbit: circular
      use_time_inferior_conjunction: True
      boundaries:
        P: [0.54, 0.56]
        K: [0.001, 20.0]
        Tc: [8899.30, 8899.40]
      priors:
        P: ['Gaussian', 0.549374, 0.000013]
        Tc: ['Gaussian', 8899.3449, 0.0008]
      spaces:
        P: Linear
        K: Linear
  activity_s01:
    model: activity
    boundaries:
      Pdec: [10.0, 100.0]
      Oamp: [0.001, 1.0]
    #priors:
    #  Oamp: ['Gaussian', 0.35, 0.035]
  activity_s02:
    model: activity
    boundaries:
      Pdec: [10.0, 100.0]
      Oamp: [0.001, 1.0]
    #priors:
    #  Oamp: ['Gaussian', 0.35, 0.035]
  star:
    star_parameters:
      boundaries:
        rotation_period: [8.0, 10.0]
      priors:
        mass: ['Gaussian', 0.76, 0.03]
        radius: ['Gaussian', 0.690, 0.036]
        density: ['Gaussian', 2.3, 0.4]
        rotation_period: ['Gaussian', 8.8, 0.1]
  common_offset_RV:
    model: common_offset
  common_offset_BIS:
    model: common_offset
  common_offset_logRHK:
    model: common_offset
models:
  radial_velocities:
    planets:
      - b
  gp_quasiperiodic_s01:
    model: gp_quasiperiodic
    common:
      - activity_s01
      - star_parameters
    use_stellar_rotation_period: True
    hyperparameters_condition: True  # Condition from Rajpaul 2017, Rajpaul+2021
    rotation_decay_condition: True # It forces the decay timescale to be at least twice the rotational period
    boundaries:
      Hamp: [0.0, 100.0] # same range for all datasets
  gp_quasiperiodic_s02:
    model: gp_quasiperiodic
    common:
      - activity_s02
      - star_parameters
    use_stellar_rotation_period: True
    hyperparameters_condition: True  # Condition from Rajpaul 2017, Rajpaul+2021
    rotation_decay_condition: True # It forces the decay timescale to be at least twice the rotational period
    boundaries:
      Hamp: [0.0, 100.0] # same range for all datasets
```


## Differences among various parametrizations

`PyORBIT`, following [Rajpaul et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.452.2269R/abstract)

```{math}
:label: quasiperiodic_pyorbit

G(t_i, t_j) = H_\mathrm{amp}^2  \exp{ \left \{-\frac{\sin^2{[\pi(t_i - t_j)/ P_\mathrm{rot}]}}{\mathbf{2} O_\mathrm{amp} ^2} - \frac{(t_i-t_j)^2}{\mathbf{2} P_\mathrm{dec}^2} \right \} }
```


- [Rajpaul et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.452.2269R/abstract)

```{math}
:label: quasiperiodic_rajpaul

G(t_i, t_j) = \eta_1^2  \exp{ \left \{-\frac{\sin^2{[\pi(t_i - t_j)/P]}}{\mathbf{2} \lambda_p^2} - \frac{(t_i-t_j)^2}{\mathbf{2} \lambda_e^2} \right \} }
```


- [Haywood et al. 2014](https://ui.adsabs.harvard.edu/abs/2014MNRAS.443.2517H/abstract)

```{math}
:label: quasiperiodic_haywood

G(t_i, t_j) = \eta_1^2  \exp{ \left \{-\frac{ \mathbf{2} \sin^2{[\pi(t_i - t_j)/\eta_3]}}{\eta_4^2} - \left ( \frac{t_i-t_j}{\mathbf{2} \eta_2} \right )^2 \right \} }
```


- [Grunblatt et al. 2015](https://ui.adsabs.harvard.edu/abs/2015ApJ...808..127G/abstract):

```{math}
:label: quasiperiodic_grunblatt

G(t_i, t_j) = h^2  \exp{ \left \{-\frac{\sin^2{[\pi(t_i - t_j)/\theta]}}{\mathbf{2} w ^2} - \left ( \frac{t_i-t_j}{\lambda} \right )^2 \right \} }
```

- [Lopez-Morales et al. 2016](https://ui.adsabs.harvard.edu/abs/2016AJ....152..204L/abstract)

```{math}
:label: quasiperiodic_lopezmorales

G(t_i, t_j) = \eta_1^2  \exp{ \left \{-\frac{ \sin^2{[\pi(t_i - t_j)/\eta_3]}}{\eta_4^2} - \left ( \frac{t_i-t_j}{\eta_2} \right )^2 \right \} }
```