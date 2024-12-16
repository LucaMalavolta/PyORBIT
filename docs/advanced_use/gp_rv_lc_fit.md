(gp_rv_lc_fit)=

# Simultaneous GPs on photometry and spectroscopy

Imagine we have both photometric and spectroscopic time series of a rather active star. Photometry will provide essential information on the rotational period of the star and the decay time scale of activity regions, while spectroscopic activity indexes are essential to disentangle planetary and stellar signal in radial velocities.

In this example, we will analyze the planetary system around K2-141 [Malavolta et al. 2018](https://ui.adsabs.harvard.edu/abs/2018AJ....155..107M/abstract). The spectroscopic dataset is taken from [Bonomo et al. 2023](https://ui.adsabs.harvard.edu/abs/2023A%26A...677A..33B/abstract) and includes Radial velocity data, Bisector Inverse Span, and S-index, all gathered with HARPS-N.

Photometric data comes from K2 campaign 12, corrected fro systematics as described in [Vanderburg and Johnson 2014](https://ui.adsabs.harvard.edu/abs/2014PASP..126..948V/abstract) and available at the [MAST K2SFF website](https://archive.stsci.edu/hlsp/k2sff). Data from TESS [sector 42](https://tess.mit.edu/observations/sector-42/) and [sector 70](https://tess.mit.edu/observations/sector-70/) are available, but we are not going to use them in this example.

```{figure} plots/K2-141_lightcurves.png
:alt: lightcurve of K2-141
:width: 100 %
Photometric data for K2-141. The stellar activity modulation is clearly identifiable.
```

We want to model stellar activity in photometry simultaneously with the spectroscopic time series.

- For photometric data, we use the [exponential-sine periodic kernel] (ESP) implemented in 'S+LEAF'. Each dataset will have its independent covariance matrix.
- For spectroscopic data, we use the [multidimensional quasi-periodic kernel](../multidimensional_gps/md_quasiperiodic.md)
 through `tinygp`.

For the photometric data, we also include the transit model for the two planets, and a normalization factor, independent for each dataset (if you have more than one light curve). Check the [](#lightcurve_fit) section in [](#quickstart) for additional information.

For the spectroscopic data, we only have to include the `radial_velocities` model due to the planets in the radial velocity dataset.

The input section takes this form:

```{code-block} yaml
:lineno-start: 1
inputs:
  LC_K2_c12:
    file: datasets/K2-141_K2_KSFF_PyORBIT.dat
    kind: Phot
    models:
      - spleaf_esp
      - lc_model
      - normalization_factor
  RVdata:
    file: datasets/K2-141_RV_PyORBIT.dat
    kind: RV
    models:
      - radial_velocities
      - gp_multidimensional
  BISdata:
    file: datasets/K2-141_BIS_PyORBIT.dat
    kind: BIS
    models:
      - gp_multidimensional
  Sdata:
    file: datasets/K2-141_Sindex_PyORBIT.dat
    kind: S_index
    models:
      - gp_multidimensional
```

In the `common` section, we specify:
- the two transiting **planets**:  `planets:b` and `planet:c`. Note that we do not specify any prior on the period or transit time, as they will be fitted directly in photometric data;
- **activity**: the rotational period of the star `Prot`, the decay timescale of activity regions `Pdec`, and the coherence scale `Oamp`, will be shared among the trained ESP and the multidimensional QP.
- **stellar parameters**: in addition to he priors on mass, radius, and density, we also include priors on the limb darkening coefficients for the two instruments involved in our analysis;


```{code-block} yaml
:lineno-start: 25
common:
  planets:
    b:
      orbit: circular
      use_time_inferior_conjunction: True
      boundaries:
        P: [0.2750, 0.2850]
        K: [0.001, 20.0]
        Tc: [57744.00, 57744.10]
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
      spaces:
        P: Linear
        K: Linear
  activity:
    boundaries:
      Prot: [10.0, 20.0]
      Pdec: [20.0, 1000.0]
      Oamp: [0.001, 1.0]
  star:
    star_parameters:
      priors:
        mass: ['Gaussian', 0.708, 0.028]
        radius: ['Gaussian', 0.681, 0.018]
        density: ['Gaussian', 2.65, 0.08]
    limb_darkening:
      type: ld_quadratic
      priors:
        ld_c1: ['Gaussian', 0.68, 0.10]
        ld_c2: ['Gaussian', 0.05, 0.10]
```

In the `models` section we include two models for Gaussian Process regression: the `spleaf_esp` will perform GP regression on each light curve independently (if you have more than one light curve), i.e., each light curve will have its own covariance matrix and its own amplitude of the covariance, similarly when performing the *trained* GP regression.

The `gp_multidimensional` will perform multidimensional GP regression over the datasets that are employing these model, namely for this specific example, the RV dataset, the BIs dataset, and the S index.

You may notice that both `spleaf_esp` and `gp_multidimensional` are recalling the same common model `common:activity`. The GP hyperparameters will be shared, in this sense we may say that photometry is training the hyperparameters of the multidimensional GP.

```{code-block} yaml
:lineno-start: 67
models:
  radial_velocities:
    planets:
      - b
      - c
  lc_model:
    kind: pytransit_transit
    limb_darkening: limb_darkening
    planets:
      - b
      - c
  normalization_factor:
    kind: local_normalization_factor
    boundaries:
      n_factor: [0.9, 1.1]
  spleaf_esp:
    model: spleaf_esp
    common: activity
    n_harmonics: 4
    hyperparameters_condition: True
    rotation_decay_condition: True
    LC_K2_c12:
      boundaries:
        Hamp: [0.000, 1.00]
  gp_multidimensional:
    model: tinygp_multidimensional_quasiperiodic
    common: activity
    hyperparameters_condition: True
    rotation_decay_condition: True
    RVdata:
      boundaries:
        rot_amp: [-100.0, 100.0] #at least one must be positive definite
        con_amp: [0.0, 100.0]
      derivative: True
    BISdata:
      boundaries:
        rot_amp: [-100.0, 100.0]
        con_amp: [-100.0, 100.0]
      derivative: True
    Sdata:
      boundaries:
        con_amp: [-1.0, 1.0]
      derivative: False
```

The rest of the configuration file looks as usual.

```{code-block} yaml
:lineno-start: 129
parameters:
  Tref: 59200.00
  safe_reload: True
  use_tex: False
  low_ram_plot: True
  plot_split_threshold: 1000
  cpu_threads: 32
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

## Sharing only some of the hyperparameters

:::{admonition} Improved interface starting from version 10.9
The following example applies to analysis performed with `PyORBIT` version >= 10.9
Previous version salready implemented the following analysis, with a less user-friendly interface. 
:::

The previous configuration assumes that you want to share the rotational period of the star `Prot`, the decay time scale of active regions `Pdec`, and the coherence scale `Oamp` among all the datasets. There are reasons to believe that the `Pdec` may change with time, and `Oamp` may change from datasets to dataset, depending on the specific aspects of stellar activity captured by each observational method.
To share only some of the hyperparameters involved in the GP regression, we associated the photometric and the spectroscopic models to two separate common activity models. In this way, the hyperparameters will be independent.
Since it would not make sense to have *all* the hyperparameters being independent (as it wold be equivalent to run two independent fits), we specify that the rotational period of the star `Prot` will be extracted from the stellar properties rather then from the activity model by setting the `use_stellar_rotation_period` flag to `True`:

```{code-block} yaml
  activity_photometry:
    model: activity
    use_stellar_rotation_period: True
    boundaries:
      #Prot: [10.0, 20.0]
      Pdec: [20.0, 1000.0]
      Oamp: [0.001, 1.0]
  activity_spectroscopy:
    model: activity
    use_stellar_rotation_period: True
    boundaries:
      #Prot: [10.0, 20.0]
      Pdec: [20.0, 1000.0]
      Oamp: [0.001, 1.0]
```

In this way, there will be independent `Pdec` and `Oamp` depending on which *common* activity model each model will recall.

Note how we have to specify two different names for the models, `activity_photometry` and `activity_spectroscopy`, while declaring that they are both `model:activity`. Since the rotational period hyperparameter is not taken from the activity model anymore, there is no need to specify its boundaries in the activity model. Instead , you can specifiy boundaries, prior, and space foir the `rotation_period` in the `star_parameters` section. 

```{code-block} yaml
  star:
    star_parameters:
      boundaries:
        rotation_period: [10.0, 20.0]
      priors:
        mass: ['Gaussian', 0.708, 0.028]
        radius: ['Gaussian', 0.681, 0.018]
        density: ['Gaussian', 2.65, 0.08]
    limb_darkening:
      type: ld_quadratic
      priors:
        ld_c1: ['Gaussian', 0.68, 0.10]
        ld_c2: ['Gaussian', 0.05, 0.10]
   #use_stellar_activity_decay
   #activity_decay
```

Finally, remember to associate the corresponding *common* activity to the model used for each dataset, as in the example below (with some omessions to highlight the differences with respect to the previous code block). 

```{code-block} yaml
models:
  spleaf_esp:
    model: spleaf_esp
    common: activity_photometry 
    n_harmonics: 4
    hyperparameters_condition: True
    ...
  gp_multidimensional:
    model: tinygp_multidimensional_quasiperiodic
    common: activity_spectroscopy
    hyperparameters_condition: True
    ...
```

The `rotation_period` parameter will be taken automatically from the stellar porperties. If there is only one stellar model, the association will be automatic without requiring the explicit declaration in the `models:*model_name*:common` section.

If you want to share also the decay time scale `Pdec`, you can follow the same procedure as the rotational period `Prot`, noting that the flag to be activated is `use_stellar_activity_decay` rather than `use_stellar_rotation_period`. In the example below, both the `Prot` and `Pdec` are shared among the models.

```{code-block} yaml
  activity_photometry:
    model: activity
    use_stellar_rotation_period: True
    use_stellar_activity_decay: True
    boundaries:
      #Prot: [10.0, 20.0]
      #Pdec: [20.0, 1000.0]
      Oamp: [0.001, 1.0]
  activity_spectroscopy:
    model: activity
    use_stellar_rotation_period: True
    use_stellar_activity_decay: True
    boundaries:
      #Prot: [10.0, 20.0]
      #Pdec: [20.0, 1000.0]
      Oamp: [0.001, 1.0]
```

In the stellar properties, the decay time scale takes the name of `activity_decay`:

```{code-block} yaml
  star:
    star_parameters:
      boundaries:
        rotation_period: [10.0, 20.0]
        activity_decay: [20.0, 1000.0]
      priors:
        mass: ['Gaussian', 0.708, 0.028]
        ...
```

The `models` ection would be the same as in the case of the shared rotational period. 

Finally, it is possible to create mixed cases where a given hyperparameters is shared among some datasets but not with others. In the following (fictional) examples, the `Prot` is shared among all the parameters, while `Pdec' is shared only among the photometric models.

```{code-block} yaml
  activity_ground_photometry:
    model: activity
    use_stellar_rotation_period: True
    use_stellar_activity_decay: True
    boundaries:
      #Prot: [10.0, 20.0]
      #Pdec: [20.0, 1000.0]
      Oamp: [0.001, 1.0]
  activity_space_photometry:
    model: activity
    use_stellar_rotation_period: True
    use_stellar_activity_decay: True
    boundaries:
      #Prot: [10.0, 20.0]
      #Pdec: [20.0, 1000.0]
      Oamp: [0.001, 1.0]
activity_spectroscopy:
    model: activity
    use_stellar_rotation_period: True
    boundaries:
      #Prot: [10.0, 20.0]
      Pdec: [20.0, 1000.0]
      Oamp: [0.001, 1.0]
```

Note that omitting the declaration for `use_stellar_activity_decay` is equivalent to declare `use_stellar_activity_decay: False`.

:::{tip}
For now, there is no flag to share the coherence scale `Oamp`, as this value is the last one you want to share among models. In other words, if you inclided to share `Oamp`, then likely you wantto share `Prot` and `Pdec` as well, following back to the case highlighted at the beginning of this page.
::: 