(TTV_and_RV_fit)=

# TTV and RV fit 

In this example we will model the planetary system around TOI-1130 ([Huang et al. 2020](https://ui.adsabs.harvard.edu/abs/2020ApJ...892L...7H/abstract)) , using the Transit Timing Variations (TTVs) from [Borsato et al. 2024] (https://ui.adsabs.harvard.edu/abs/2024A%26A...689A..52B/abstract) and the Radial Velocity (RV)measurements from [Korth et al. 2023](https://ui.adsabs.harvard.edu/abs/2023A%26A...675A.115K/abstract)

```{note}
Both `PyORBIT` and `TRADES` measure and predict the *times of inferior cinjunction $T_c$* , which may differ from the *central times of transit* for eccentric orbits. For short, we will call the former **transit times**
```

## Preparing the datasets 

The dataset must be prepared following the instructions in <project:../prepare_datasets.md#times-of-inferior-conjunction-transit-time> for the TTV files, and <project:../prepare_datasets.md#how-to-properly-format-a-dataset> for the RV files.

Keep in mind the following facts when preparing the dataset files:
- The code will assume that the transit time measurements included in a file belong to the same planet, so you cannot mix $T_c$ measurements of different planets within the same file.
- You cannot include in the same file two $T_c$ measurements referring to the same transits, e.g. from two different instruments.

Each dataset will have its corresponding jitter parameter. To use the same jitter parameter across different datasets, use the `common_jitter` model.

```{eval-rst}
.. code-block:: yaml
  :linenos:

  inputs:
    TCdata_b:
      file: ../dataset_pyorbit/borsato24_t0_planetb_MOD.dat
      kind: transit_time
      models: ttv_b
    TCdata_c:
      file: ../dataset_pyorbit/borsato24_t0_planetc_MOD.dat
      kind: transit_time
      models: ttv_c
    HARPS_RVdata:
      file: ../dataset_pyorbit/korth23_harps_rvs.dat
      kind: RV
      models:
        - radial_velocities
    PFS_RVdata:
      file: ../dataset_pyorbit/korth23_pfs_rvs.dat
      kind: RV
      models:
        - radial_velocities
```

## The planet model

When modelling the radial velocities or the transit of non-interacting planets, the relative geometry of the orbits (how the orbital planes are positioned with respect to each other) is not relevant, and that allows us to fix several orbital parameters without affecting the final outcome. For example, it is customary to fix the longitude of the ascending note $\Omega$ (either to $0^{\circ}$ or $180^{\circ}$) and set the eccentricity to zero (circular orbit). For the modelling of radial velocity, we don't even need to know the orbital inclination or the mass of the planets.

Dynamical modelling works differently. First of all, we need to include the true masses of the planets (and therefore also the mass of the star) to properly compute how they interact with each other. We must also know how the orbits are oriented with each other, using one of the planet as reference. These are the orbital parameter we have to include for each planet:

| Name        | Symbol | Parameter | Unit     |
| :---        | :-- | :----     | :---     |
| P      | $P$ | Orbital period of the planet                      | days     |
| mean_long | $L_0$ | Mean longitude of the orbit at $T_{\rm ref}$  | degrees |
| e      | $e$ | eccentricity of the orbit | adimensional |
| omega  | $\omega_p$ | argument of periastron of the *planet*  | degrees |
| M_Me | $M_{\rm p}$ | planet mass in Earth masses | $M_\oplus$ |
| i | $i$ | orbital inclination with respect to the plane of the sky | degrees |
| Omega  | $\Omega$ | longitude of the ascending note | degrees |
%| R_Rs | $k$ | planet radius in stellar radii | $R_\star$ |

To force the use of the mass and the orbital inclination, we have to activate the corresponding flags in the  `planet` model:

```{code-block} yaml
    planets:
      b:
        use_mass: True
        use_inclination: True
```

In transit fitting, the maximum value allowed for the inclination is set to $90^{\circ}$, as two inclinations equidistant from $90^{\circ}$ would produce the same signal. For the same reason, the impact parameter $b$ is forced to be positive.

In dynamical modelling, mutual inclinations play a role: two planets $i_b=87^{\circ}$ and $i_c=93^{\circ}$ would produce the same transit signal of two planets with $i_b=i_c=87^{\circ}$, but different dynamical signature.

After fixing the quadrant of the reference planet by forcing its inclination below $90^{\circ}$, we must allow for the other planets to expler both quadrants. If we are using the impact parameter $b$, we have to allow for negative values for all the planets except the first one.

For similar reasons, we fix the longitude of the ascending node the first planet $\Omega$ to $180^{\circ}$, and leave it free for the other planets.

In the example below, planet $b$ acts as the reference planet:

```{code-block} yaml
  :lineno-start: 20

  common:
    planets:
      b:
        orbit: dynamical
        parametrization: Eastman2013
        use_mass: True
        use_inclination: True
        boundaries:
          P: [4.00, 4.10]
          M_Me: [10., 30.0]
          i: [85., 90.]
          e: [0.00, 0.20]
        priors:
          i: ['Gaussian', 87.5, 0.1]
          P: ['Gaussian', 4.074554, 0.001]
        spaces:
          P: Linear
          M_Me: Linear
        fixed:
          Omega: [180.000, 0.001]
      c:
        orbit: dynamical
        parametrization: Eastman2013
        use_mass: True
        use_inclination: True
        boundaries:
          P: [8.30, 8.40]
          M_Me: [300.0, 350.0]
          i: [85., 95.]
          e: [0.00, 0.20]
        priors:
          i: ['Gaussian', 87.6, 0.1]
          P: ['Gaussian', 8.3501898, 0.0001]
        spaces:
          P: Linear
          M_Me: Linear
        fixed:
          Omega: [180.000, 0.001]
```