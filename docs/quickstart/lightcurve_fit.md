(lightcurve_fit)=

# Light curve fit

We will start this example by considering the [TESS](https://www.nasa.gov/tess-transiting-exoplanet-survey-satellite) light curve of HD189733, also known as [TOI 4470](https://exofop.ipac.caltech.edu/tess/target.php?id=256364928). 

For the purposes of this example, the light curve has been already *flattened* using `wōtan` [Hippke et al. 2019](https://ui.adsabs.harvard.edu/abs/2019AJ....158..143H/abstract), although one could perform light curve flattening within `PyORBIT` as well.

## Dataset 

Light curve datasets differ from radial velocity or activity datasets for a fundamental aspect: the effect of planets or stellar activity is a *multiplicative* effect, rather than *additive*. For example, a Jupiter-like transiting a Sun-like star will cause a drop in light of about 1% regardless of the absolute value of the flux recorded on the ground. Thus, removing an additive offset from the light curve would inevitably affect the observed transit depth. For this reason, a light curve dataset *must always have* the offset flag deactivated:

```
#epoch value error jitter_flag offset_flag
  2459422.028697         0.999750         0.000205 0 -1
  2459422.030086         1.000049         0.000205 0 -1
  2459422.031475         1.000001         0.000205 0 -1
...
```

This is equivalent in writing:
```
#epoch value error jitter_flag
  2459422.028697         0.999750         0.000205 0
  2459422.030086         1.000049         0.000205 0
  2459422.031475         1.000001         0.000205 0
...
```

Keep in mind that **0** activates a flag, **-1** deactives it, folllowing the Python notation, while the absence of a column automatically deactivate a flag.


## Configuration file 

```{code-cell} yaml
:lineno-start: 1

inputs:
  LCdata_TESS:
    file: ./HD189733_TESS_PyORBIT.dat
    kind: Phot
    models:
      - lc_model
common:
  planets:
    b:
      orbit: circular
      parametrization: Eastman2013_Tcent
      boundaries:
        P: [2.2185600, 2.2185800]
        Tc: [2459770.4100, 2459770.4110]
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

There is a lot to process:
