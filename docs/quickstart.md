(quickstart)=

# Quickstart 

Although it started as a simple program, `PyORBIT` got increasingly complex with time. 
Before you start delving into the many available models, let’s work together on a simple test case to better understand the philosophy behind `PyORBIT`. In this first example, we will fit just the radial velocity measurements of HD189733 from
[Bouchy et al. 2006](https://ui.adsabs.harvard.edu/abs/2005A%26A...444L..15B/abstract). 

To get started, we need two ingredients: the **dataset**, and a **configuration file**. 

## Dataset

The radial velocity dataset is just a text file with five columns, as in this example: 

```
#epoch value error jitter_flag offset_flag
  2453611.419560   -2229.0   5.1   0   0
  2453613.430197   -2156.0   5.1   0   0
  2453614.434525   -2575.0   5.1   0   0
```

- 1<sup>st</sup> column: epoch of the observation, espresed in BJD_TDB
- 2<sup>nd</sup> column: radial velocity measurements at the given epoch
- 3<sup>rd</sup> column: error associated to the measurement
- 4<sup>th</sup> column: flag to activate the jitter parameter(s)
- 5<sup>th</sup> column: flag to activate the offset parameter(s)

The **jitter** is a value added in quadrature to the error estimates to take into account underestimated error bars.

The **offset** is constant value added to the model - in this case it corresponds to the systemic velocity of the star.

Folllowing the Python notation,  **0** activates a flag, **-1** deactives it. 

Other flags and ancillary data can be included, as detailed in the *Prepare a dataset file* oage 

## Configuration file

The confoguration file is again just a text file that you can modfiy with any editor, this time written in [YAML](https://yaml.org/) language.
Given the simple model we are going to fit - it's just a sinusoid after all - the configuration file may seem rather complex! 
This structure reflects the logical organization of PyORBIT, and it will become extremely useful with increasingly more complex models. 



```{eval-rst}
.. code-block:: yaml
   :linenos:

    inputs:
        Bouchy2005_RV01_noRML:
            file: Bouchy2005_RV01_noRML_PyORBIT.dat
            kind: RV
            models:
            - radial_velocities
    common:
        planets:
            b:
            orbit: circular
            boundaries:
                P: [0.50, 5.0]
                K: [0.01, 300.0]
        star:
            star_parameters:
            priors:
                mass: ['Gaussian', 0.806, 0.048]
    models:
        radial_velocities:
            planets:
            - b
    parameters:
        Tref: 2459500.00
    solver:
        pyde:
            ngen: 50000
            npop_mult: 4
        emcee:
            npop_mult: 4
            nsteps: 20000
            nburn: 10000
            thin: 100
        nested_sampling:
            nlive: 1000
        recenter_bounds: True
```

