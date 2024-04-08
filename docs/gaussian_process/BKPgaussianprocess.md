(gaussian_process_regressionOLD)=

# Gaussian process regression

Gaussian process regression (GPR) is a nonparametric, Bayesian approach to regression which has been very succesful for the analysis of radial velocity datasets in the presence of stellar activity, e.g., [Haywood at al. 2014](https://ui.adsabs.harvard.edu/abs/2014MNRAS.443.2517H/abstract),  [Grunblatt et al. 2015](https://ui.adsabs.harvard.edu/abs/2015ApJ...808..127G/abstract). The first stable and tested implementation of GPR in `PyORBIT` dates back to 2018. 

In this section, we will discuss only models encompassing independent covariance matrix among datasets, with some hyperparameters in common. I will refer to these approach as *classic* or *standard* Gaussian processes. 



### Kernels

#### Quasi-periodic

The most common and most reliable kernel for stellar activity  is the quasi-periodic one, following the expression given by [Grunblatt et al. 2015](https://ui.adsabs.harvard.edu/abs/2015ApJ...808..127G/abstract):

```{math}
:label: quasiperiodic_grunblatt

G(t_i, t_j) = h^2  \exp{ \left \{-\frac{\sin^2{[\pi(t_i - t_j)/\theta]}}{2 w ^2} - \left ( \frac{t_i-t_j}{\lambda} \right )^2 \right \} }
```

where $\theta$ is equivalent to the rotation period of the star, $w$ is the coherence scale, and $\lambda$ is usually associated to the decay time scale fo the active regions

```{important}
It is common to have a factor 2 in the denominator of the aperiodic variation (i.e., $2 \lambda$ rather than $\lambda$) in {eq}`quasiperiodic_grunblatt`. In such a case, it is sufficient to multiply the value of $\lambda$ of `PyORBIT` by a factor $\sqrt(2)$ - keep it in mind when assigning priors!
```

This kernel relies on the `george` package. An independent implementation not relying on any package is available, however it is much slower. 
A new implementation using  `tinyGP` is now available, but it requires a few extra tricks in the configuration file and execution.


#### Quasi-periodic with cosine 

**work in progress** 

This kernel has been introduced by [Perger et al. 2021](https://ui.adsabs.harvard.edu/abs/2021A%26A...645A..58P/abstract). The kernell has been implemented in `PyORBIT` without relying on any other package.

If we define $\tau = t_i-t_j$ : 

```{math}
:label: quasiperiodic_cosine

G(\tau ) = \exp{\frac{-2 \tau^2}{\lambda^2}} *  \left [ h_1^2 \exp{-\frac {
  \sin^2{( \pi \tau / \theta )}}{2 w ^2}} + h_2^2  \cos \frac{ 4\pi \tau}{\theta} \right ]

```

Within `PyORBIT`, $h_1$ and $h_2$ have been labelled as `Hamp` and `Camp` respectively.


#### Matern 3/2

**work in progress** 

#### Damped harmonic oscillator

**work in progress** 

#### Rotation

**work in progress** 

#### Granulation and Rotation

**work in progress** 

#### Granulation, Oscillations, and Rotation

**work in progress** 


### Improvements

In `PyORBIT` unlimited number of additional datasets can be included for the simultaneous training of the hyperparameters. All the parameters will be shared except the amplitude of covariance matrix, which is dataset-dependent. Each dataset will be characterized by its own covariance matrix.

## Model definition and requirements

- model name: ``gp_quasiperiodic``
- required common objects : ``activity``

## Keywords

Model-wide keywords, with the default value in bold face.

**hyperparameters_condition**
* accepted values: `True` | **`False`**
* activate the conditions $ \lambda ^ 2 > (3/4 \pi) \theta ^2 w ^ 2 $ (adapted from [Rajpaiul 2017](https://ui.adsabs.harvard.edu/abs/2017PhDT.......229R/abstract) and [Rajpaul et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.1847R/abstract) to take into account the factor 2 in the denominator of the aperiodic variation) to ensure that the QP function has at least one non-trivial turning point.

**rotation_decay_condition**
* accepted values: `True` | **`False`**
* if activated, it ensures that the decay time scale of the activity regions $\lambda$ is at least twice the rotational period of the star $\theta$


## Examples

In the following example, ore RV dataset comes together with two activity indicators, the `BIS`  of the CCF and the `S index`

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



## Model parameters

The following parameters will be inherited from the common model (column *Common?: common*) or a different value will be assigned for each dataset (column *Common?: dataset*)

| Name        | Parameter | Common?  | Definition  | Notes |
| :---        | :-------- | :-------------  | :-----  | :---- |
| Prot      | rotational period of the star $\theta$ | common | ``activity``     | |
| Pdec      | Decay time scale of active regions $\lambda$ | common | ``activity``     | |
| Oamp | Coherence scale $w$ | common | ``activity`` |   |
| Hamp  | Amplitude of the kernel | dataset | ``activity``     | |
| Camp  | Amplitude of the cosine part of the kernel | dataset | ``activity``     | |

