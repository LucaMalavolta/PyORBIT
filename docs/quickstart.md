(quickstart)=

# Quickstart 

Although it started as a simple program, `PyORBIT` got increasingly complex with time. 
Before you start delving into the many available models, let’s work together on a simple test case to better understand the philosophy behind `PyORBIT`. In this first example, we will fit just the radial velocity measurements of HD189733 from
[Bouchy et al. 2006](insert_ads_link). 

To get started, we need two ingredients: the **dataset**, and a **configuration file**. 

The dataset is just a text file with six columns 



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

