(radial_velocities)=

# Planetary RV signal

A Keplerian model for non-interacting planets.
The effect of each planet on the radial velocity of the host star are modelled as non-interacting Keplerians.

## Model definition and requirements

- model name: ``radial_velocities``
- required common objects : ``planets`` , ``star``

## Keywords

There are no keywords defined for this model, all the relevant information are inherited from the planet common model

## Examples

Minimal configuration for a model

```yaml
models:
  radial_velocities:
    planets:
     - b
```

```{tip}
Planets are automatically assigned to the `star` model when there is only one stellar object specified in the model
```

## Model parameters

The following parameters will be inherited from the common model (column *Common?: common*) or a different value will be assigned for each dataset (column *Common?: dataset*)

| Name        | Parameter | Common?  | Definition  | Notes |
| :---        | :-------- | :-------------  | :-----  | :---- |
| P      | orbital period of the planet | common | ``planets``     | |
| K      | Radial velocity semiamplitude | common | ``planets``     | |
| mean_long | mean longitude of the orbit, assuming $\Omega=0$ | common | ``planets`` | (1) |
| Tc     | Time of inferior conjunction              | common | ``planets``     | (2)|
| e      | eccentricity of the orbit  | common | ``planets`` | (3) |
| omega  | argument of pericenter of the *planet* $\omega_p$  | common |  ``planets`` | (3) |
| e_coso | $ e \cos{\omega}$ in Ford 2006 parametrization | common |  ``planets`` | (4) |
| e_sino | $ e \sin{\omega}$ in Ford 2006 parametrization | common |  ``planets`` | (4) |
| sre_coso | $ \sqrt{e} \cos{\omega}$ in Eastman 2013 parametrization | common |  ``planets`` | (5)|
| sre_sino | $ \sqrt{e} \sin{\omega}$ in Eastman 2013 parametrization | common |  ``planets`` | (5) |

Notes:

  1. replaced by ``Tc`` when ``use_time_inferior_conjunction: True``
  2. ``use_time_inferior_conjunction: True`` parametrization only
  3. ``Standard`` parametrization only
  4. ``Ford2006`` parametrization only
  5. ``Eastman2013`` parametrization only