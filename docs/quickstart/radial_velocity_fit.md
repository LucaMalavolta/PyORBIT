(radial_velocity_fit)=

# Radial Velocities

In this first example, we will fit just the radial velocity measurements of HD189733 from
[Bouchy et al. 2006](https://ui.adsabs.harvard.edu/abs/2005A%26A...444L..15B/abstract). Later on we will include addiitonal datasets from the literature.

To get started, we need two ingredients: the **dataset**, and a **configuration file**.

(quickstart_rv_dataset)=
## Dataset

The radial velocity dataset is just a text file with five columns, as in this example:

```
#epoch value error jitter_flag offset_flag
  2453611.419560   -2229.0   5.1   0   0
  2453613.430197   -2156.0   5.1   0   0
  2453614.434525   -2575.0   5.1   0   0
```

- 1<sup>st</sup> column: epoch of the observation, expressed in BJD_TDB
- 2<sup>nd</sup> column: radial velocity measurements at the given epoch, in metre per second
- 3<sup>rd</sup> column: error associated to the measurement, again in metre per second
- 4<sup>th</sup> column: flag to activate the jitter parameter(s)
- 5<sup>th</sup> column: flag to activate the offset parameter(s)

The **jitter** is a value added in quadrature to the error estimates to take into account underestimated error bars.

The **offset** is constant value added to the model - in this case it corresponds to the systemic velocity of the star.

Folllowing the Python notation,  **0** activates a flag, **-1** deactives it.

Other flags and ancillary data can be included, as detailed in the [Prepare a dataset file](prepare_dataset)

## Configuration file

The configuration file is again just a text file that you can modify with any editor, this time written in [YAML](https://yaml.org/) language.
Similarly to Python, indentation is quite important, and it is the major cause of errors for beginners. Whitespace is part of YAML's formatting, with *tabs* specifically forbidden by the language's specifications.

Given the simple model we are going to fit - it's just a sinusoid after all - the configuration file may seem rather complex!
This structure reflects the logical organization of PyORBIT, and it will become extremely useful with increasingly more complex models.

Look at the example below:

```{code-block} yaml
:lineno-start: 1

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

You can see that there are five main sections:

- ``inputs``: the properties of each dataset are enclosed within a dedicated subsection.
- ``common``: it encompasses two subsections, *planets* and *star*, where you list the characteristics of the objects that will be common to all the models.
- ``models``: any model used to fit your data, to be recalled in the ``model`` subsection of the dataset
- ``parameters``: values that are not system parameters per se, but are required to properly perform the analysis.
- ``solver``: the parameters required by the samplers are all listed here.

In this example, we will fit the ``Bouchy2005_RV01_noRML`` dataset with ``radial_velocity`` model encompassing only the planet ``b``. Note that the stellar mass is used only to convert the RV semi-amplitude into physical mass, while keeping into account the error associated to the stellar mass itself.

For more details on the configuration files, don't forget to check the [Prepare a configuration file](prepare_yaml) page.

## Launch PyORBIT

Analyzing your dataset with `PyORBIT` is a two-step process: you first run the actual analysis, and then you extract the relevant information (posteriors distributions of the derived parameters, corner plots...).

If you have installed `PyORBIT` through pip (from PyPI or from the local repository), just run this command on your terminal:

```{code} bash
pyorbit_run emcee configuration_file.yaml
```

The first argument of `pyorbit_run` is the sampler to be used in the analysis (see [Samplers](samplers)), the second argument is the configuration file.

PyORBIT will produce a very extensive terminal output, detailed in the section [Interpreting the output](terminal_output).

To save the terminal output to a file, add  `> configuration_file_emcee_run.log` to the command and the terminal output will be saved to a file.

After the analysis is complete, you can run convergence diagnostics and write model files for plotting by running the command:

```{code} bash
pyorbit_results emcee configuration_file.yaml
```

The confidence intervals of the main parameters you are interested to are printedin the section called *Statistics on the physical parameters obtained from the posteriors samples*

```text
====================================================================================================
     Statistics on the physical parameters obtained from the posteriors samples
====================================================================================================

----- dataset:  Bouchy2005_RV01_noRML
offset_0               -2363.5         -4.8          4.4    (15-84 p)
jitter_0                  16.0         -3.2          5.0    (15-84 p)


----- common model:  b
P                      2.21511     -0.00027      0.00028    (15-84 p)
K                        204.4         -7.4          7.4    (15-84 p)
e                         0.00                              [fixed]
omega                90.000000                              [fixed]
mean_long                  264         -117          127    (15-84 p)
```

More details on the interpretation of the results can be found in the page [Interpreting the output](results_interpretation).

(note_on_orbital_parametrization)=
## Note on the orbital parametrization

Unless otherwise specified, `PyORBIT`  assumes that the orbital period of the planet is not knwon.
If you want to explore a large range of period, you should avoid the parametrization introduced by [Eastman et al. 2013](https://ui.adsabs.harvard.edu/abs/2013PASP..125...83E/abstract) involving the time of inferior conjuction $T_C$ (sometimes confused with the time of central transit) as you may end up with several peaks in the posterior of $T_C$ spaced according to the orbital period of the planet (formally correct, but not optimal for the computational point of view).
In this case I suggest to use the parametrization introduced by [Ford 2006](https://ui.adsabs.harvard.edu/abs/2006ApJ...642..505F/abstract) in section 4.3. In `PyORBIT`, the sum of the *argument of pericenter of the planet $\omega_p$* and the *mean anomaly at the reference time $\M_0$* is called *mean longitude*, as the sum of the two angles is indeed equal to the formal defitiion of mean longitude when the *ascending node* is equal to zero. With this parametrization, you are not forced to use ad-hoc boundaries for the $T_C$ to avoid multiple peaks in its posterior, as the problem is solved at its roots.

If you want to fit for an eccentric orbit, you can do it by changing the orbit keyword to keplerian. You

```{code-block} yaml
:lineno-start: 1

common:
    planets:
        b:
            orbit: keplerian
            parametrization: Standard
            boundaries:
                P: [0.50, 5.0]
                K: [0.01, 300.0]
                e: [0.00, 0.95]
```

There are three available parametrization for eccentricity $e$ and argument of pericenter of the planet $\omega_p$:
* `Standard`:  $e$ and $\omega_p$ are fitted directly,
* `Ford2006`: the free parameters are $e \sin{\omega_p}$ and
  $e \cos{\omega_p}$, as proposed by [Ford 2006](https://ui.adsabs.harvard.edu/abs/2006ApJ...642..505F/abstract). 
* `Eastman2013`: the free parameters are $\sqrt{e} \sin{\omega_p}$ and 
  $\sqrt{e} \cos{\omega_p}$. This is the default choice for non-circular orbits. 
In all cases it is possible to put a prior on the eccentricity, even if it is a derived parameter.

## Multiple datasets

To add more radial velocity datasets, you just need to include them in the `input` sections alongside with the model you want to use to analyze the, - in this case, still `radial_velocities`. The rest of the configuration file is unchanged.

In this example below, we have added RV measurements from [Winn et al. 2006](https://ui.adsabs.harvard.edu/abs/2006ApJ...653L..69W/abstract) and [Boisse et al. 2009](https://ui.adsabs.harvard.edu/abs/2009A%26A...495..959B/abstract).

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
common:
   ...

```

`PyORBIT` will automatically include a *jitter* term and an *offset* term for each dataset, if the data file has been formatted as described in the [Dataset](quickstart_rv_dataset) section.


```text
====================================================================================================
     Statistics on the model parameters obtained from the posteriors samples
====================================================================================================

----- dataset:  Bouchy2005_RV01_noRML
offset_0               -2362.8         -4.6          4.5    (15-84 p)
jitter_0                  16.2         -3.1          4.3    (15-84 p)


----- dataset:  Winn2006_RV02_noRML
offset_0                 -11.0         -2.1          2.3    (15-84 p)
jitter_0                  10.3         -1.4          1.8    (15-84 p)


----- dataset:  Boisse2009_RV03_noRML
offset_0               -2273.6         -1.6          1.7    (15-84 p)
jitter_0                   8.6         -1.2          1.5    (15-84 p)


----- common model:  b
omega                90.000000                              [fixed]
P                     2.218551    -0.000022     0.000022    (15-84 p)
e                         0.00                              [fixed]
K                        202.0         -1.8          1.8    (15-84 p)
mean_long                140.3         -8.9          8.9    (15-84 p)
```

The inclusion of the two datasets has lowered the uncertainty on the radial velocity semiamplitude of planet b from $K_\mathrm{b} = 204.4 \pm 7.4 \, [\mathrm{m/s}] $ to $K_\mathrm{b} = 202.0 \pm 1.8 \, [\mathrm{m/s}] $. The conversion in Jupiter masses is also provided in the `Statistics on the derived parameters obtained from the posteriors samples`, if a prior on the stellar mass is provided in the configuration file:

```text
====================================================================================================
     Statistics on the derived parameters obtained from the posteriors samples
====================================================================================================

----- common model:  b
Inclination fixed to 90 deg!
Computing exact mass of the planet (average approximate mass larger than 30.0 Me)
M_Mj                     1.122       -0.045        0.047    (15-84 p)
M_Me                       357          -14           15    (15-84 p)
Tc                 2459501.908       -0.055        0.055    (15-84 p)
a_AU_(M)               0.03099     -0.00063      0.00061    (15-84 p)

```

```{warning}
The error on the mass wil ltake into account the error associated to the stellar mass.
However, If no information on the orbital inclination is provided, `PyORBIT` will return the *minimum mass* of the planet.
```

Orbital inclination with the associated error can be included as a `fixed` parameter - in the sense that the parameter is not involved in the computation fo the log-likelihood and it will not be optimized

```{code-block} yaml
:lineno-start: 1
common:
  planets:
    b:
      orbit: circular
      boundaries:
        P: [0.50, 5.0]
        K: [0.01, 300.0]
      fixed:
        i: [85.580, 0.060]
  star:
    ...
```

You can see that the mass as *slightly* increased as a result of the updated inclination

```
====================================================================================================
     Statistics on the derived parameters obtained from the posteriors samples
====================================================================================================

----- common model:  b
Inclination randomized to 85.58 +- 0.06 deg
Computing exact mass of the planet (average approximate mass larger than 30.0 Me)
i                       85.579       -0.059        0.059    (15-84 p)
M_Mj                     1.128       -0.046        0.045    (15-84 p)
M_Me                       359          -15           14    (15-84 p)
Tc                 2459501.908       -0.055        0.055    (15-84 p)
a_AU_(M)               0.03100     -0.00062      0.00060    (15-84 p)

```

```{tip}
When `fixed`, you can change the value of the orbital inclination and obtain updated masses by just running `pyorbit_results` again. 
```


