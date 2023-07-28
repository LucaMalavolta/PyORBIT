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

Keep in mind that **0** activates a flag, **-1** deactivates it, following the Python notation, while the absence of a column automatically deactivates a flag.

## Configuration file

In this example, the full configuration file to fit a TESS lightcurve is reported

```{code-block} yaml
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
      use_time_inferior_conjunction: True
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

**Fit of the time of inferior conjuction $T_c$** (equivalent to the central time of transit in the case of a circular orbit) is the way to go, as we now we have a good guess of the orbital period and time of transit. Although the $P$ and $T_C$ boundaries are quite tight, they can still be considered uninformative priors as the final posteriors will have much narrow distributions.
To enable the use of $T_C$ instead of the mean longitude `mean_long`, you need to activate the flag `use_time_inferior_conjunction`. 

```{code-block} yaml
:lineno-start: 7

common:
    planets:
        b:
            orbit: circular
            use_time_inferior_conjunction: True
            boundaries:
                P: [0.50, 5.0]
                K: [0.01, 300.0]
                e: [0.00, 0.95]
```

**Stellar density** is expressed in Solar units, so a star with one Solar mass and one Solar radius will have a density equal to one. If you know mass anda radius of a star in Solar units, then the density will be simply $\rho_\star = M_\star / R_\star^3$, without multiplicative constants.

**Limb darkening** coefficients are included under the `star` common model for conceptual reasons, although you must remeber that lightcurves obtained with different filters will require specific limb darkening paramters. When using a *quadratic* limb darkening law, you can use the parametrization introduced by [Kipping](https://ui.adsabs.harvard.edu/abs/2013MNRAS.435.2152K/abstract) cby simply activating the flag as in the example

**Light curve modellling** can be performed either with `batman` ([Kreidberg 2015](https://ui.adsabs.harvard.edu/abs/2015PASP..127.1161K/abstract)) or `PyTransit` ([Pairviainen 2015](https://arxiv.org/abs/1504.07433)). You can choose the model by specifying `batman_transit` or `pytransit_transit` respectively.

**Linear space for Period** will avoid the logarithmic transformation of this parameter, which is not required given the small range of the period.

## Multiband photometry 

Multiband photometry with specific limb darkening parameters is easy to accomodate.
In this example, we will add four lightcurves in Ic(Cousins) band from [Bakos et al. 2006](https://ui.adsabs.harvard.edu/abs/2006ApJ...641L..57B/abstract) to the TESS photometry

```{code-block} yaml
:lineno-start: 1

inputs:
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
    file: /Bakos2006_LC08_PyORBIT.dat
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

Note the main differences:

- There are now two limb darkening models under the `star` common model. I decided to rename the TESS limb darkening model to clearly distinguish from the Ic (Cousins) one, but the renaming is not strictly necessary.
- I have specified a set of prior for the limb darkening coefficient of the Ic (Cousins) filter. In doing so, there is no advantage in using the Kipping parametrization and it can be omitted
- There are now two light curve models under the `models` section, one for each set of limb darkening paramters. The number of models does not depend of the number of datasets, but on the number of photometric bands, as they influence the shape of a transit.
- In the `input` section, the correct model must be associated to each dataset.

```{attention} In this example we assumed that the transit depth is independent from wavelength.
```

It is indeed possible to obtain independent radius estimates for different photometric bands by using some of the advanced features of `PyORBIT`. An example will be provided in a dedicated tutorial. 

```{warning}
`batman` may get stuck in a loop when eccentricity is higher than 0.95. 

```