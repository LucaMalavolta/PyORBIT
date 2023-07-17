(combined_fit)=

# Combined radial velocity and light curve fit 

To perform a combined radial velocity and photometry fit we just need to merge the two configuration files:

```{code-block} yaml
:lineno-start: 1

inputs:
  LCdata_TESS:
    file: ./HD189733_TESS_PyORBIT.dat
    kind: Phot
    models:
      - lc_model
  Bouchy2005_RV01_noRML:
    file: Bouchy2005_RV01_noRML_PyORBIT.dat
    kind: RV
    models:
      - radial_velocities
common:
  planets:
    b:
      orbit: circular
      use_time_inferior_conjunction: True
      boundaries:
        P: [2.2185600, 2.2185800]
        Tc: [2459770.4100, 2459770.4110]
        K: [0.01, 300.0]
      spaces:
        P: Linear
  star:
    star_parameters:
      priors:
        mass: ['Gaussian', 0.806, 0.048]
        radius: ['Gaussian', 0.756, 0.018]
        density: ['Gaussian', 1.864, 0.175] #in Solar unit!!!!!!!
    limb_darkening:
      model: ld_quadratic
      parametrization: Kipping
models:
  radial_velocities:
    planets:
      - b
  lc_model:
    model: batman_transit
    limb_darkening: limb_darkening
    planets:
      - b
parameters:
  Tref: 2459750.00
solver:
  pyde:
    ngen: 50000
    npop_mult: 4
  emcee:
    npop_mult: 4
    nsteps: 100000
    nburn: 20000
    nsave: 10000
    thin: 100
  nested_sampling:
    nlive: 1000
  recenter_bounds: True
```

In other words:

- We include both datasets in the `input` sections
- We specify the boundaries for the planetary parameters that are involved in either the radial velocity modelling or the light curve modelling
- We specify both models in `models`.

The code will compute the likelihood using all the available data automatically, delivering a single set of orbital parameters able to reproduce the observed characteristics in all datasets.

Finally, you can add as many datasets and model as you wish:

```{code-block} yaml
:lineno-start: 1

inputs:
  Bouchy2005_RV01_noRML:
    file: Bouchy2005_RV01_noRML_PyORBIT.dat
    kind: RV
    models:
      - radial_velocities
  Winn2006_RV02_noRML:
    file: Winn2006_RV02_noRML_PyORBIT.dat
    kind: RV
    models:
      - radial_velocities
  Boisse2009_RV03_noRML:
    file: Boisse2009_RV03_noRML_PyORBIT.dat
    kind: RV
    models:
      - radial_velocities
  LCdata_TESS:
    file: HD189733_TESS_PyORBIT.dat
    kind: Phot
    models:
      - lc_model_TESS
  LCdata_Bakos2006_LC06:
    file: Bakos2006_LC06_PyORBIT.dat
    kind: Phot
    models:
      - lc_model_Ic_Cousins
  LCdata_Bakos2006_LC07:
    file: Bakos2006_LC07_PyORBIT.dat
    kind: Phot
    models:
      - lc_model_Ic_Cousins
  LCdata_Bakos2006_LC08:
    file: Bakos2006_LC08_PyORBIT.dat
    kind: Phot
    models:
      - lc_model_Ic_Cousins
  LCdata_Bakos2006_LC09:
    file: Bakos2006_LC09_PyORBIT.dat
    kind: Phot
    models:
      - lc_model_Ic_Cousins
common:
  planets:
    b:
      orbit: keplerian
      parametrization: Eastman2013
      use_time_inferior_conjunction: True
      boundaries:
        P: [2.2185600, 2.2185800]
        Tc: [2459770.4100, 2459770.4110]
        K: [0.01, 300.0]
        e: [0.00, 0.95]
      spaces:
        P: Linear
  star:
    star_parameters:
      priors:
        mass: ['Gaussian', 0.806, 0.048]
        radius: ['Gaussian', 0.756, 0.018]
        density: ['Gaussian', 1.864, 0.175] #in Solar unit!!!!!!!
    limb_darkening_TESS:
      model: ld_quadratic
      parametrization: Kipping
    limb_darkening_Ic_Cousins:
      type: ld_quadratic
      #parametrization: Kipping
      priors:
        ld_c1: ['Gaussian', 0.45, 0.05]
        ld_c2: ['Gaussian', 0.13, 0.05]
models:
  radial_velocities:
    planets:
      - b
  lc_model_TESS:
    model: batman_transit
    limb_darkening: limb_darkening_TESS
    planets:
      - b
  lc_model_Ic_Cousins:
    kind: batman_transit
    limb_darkening: limb_darkening_Ic_Cousins
    planets:
      - b
parameters:
  Tref: 2459750.00
solver:
  pyde:
    ngen: 50000
    npop_mult: 4
  emcee:
    npop_mult: 4
    nsteps: 100000
    nburn: 20000
    nsave: 10000
    thin: 100
  nested_sampling:
    nlive: 1000
  recenter_bounds: True
```