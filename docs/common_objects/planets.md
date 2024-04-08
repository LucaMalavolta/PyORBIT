(planets)=

# Planets

## Parameters definition

Each *planet* is defined by the following parameters

| Name        | Parameter | Unit     |
| :---        | :----     | :---     |
| P      | Orbital period of the planet                      | days     |
| K      | Radial velocity semiamplitude        | m/s      |
| Tc     | Central time of transit              | days      |
| mean_long | Mean longitude of the orbit, assuming $\Omega=0$ | degrees |
| e      | eccentricity of the orbit | adimensional |
| omega  | argument of periastron of the *planet* $\omega_p$  | degrees |
| e_coso | $ e \cos{\omega_p}$ in [Ford 2006](https://ui.adsabs.harvard.edu/abs/2006ApJ...642..505F/abstract) parametrization | adimensional |
| e_sino | $ e \sin{\omega_p}$ in [Ford 2006](https://ui.adsabs.harvard.edu/abs/2006ApJ...642..505F/abstract) parametrization | adimensional |
| sre_coso | $ \sqrt{e} \cos{\omega_p}$ in [Eastman et al. 2013](https://ui.adsabs.harvard.edu/abs/2013PASP..125...83E/abstract) parametrization | adimensional |
| sre_sino | $ \sqrt{e} \sin{\omega_p}$ in [Eastman et al. 2013](https://ui.adsabs.harvard.edu/abs/2013PASP..125...83E/abstract) parametrization | adimensional |
| M_Me | planet mass in Earth masses | $M_\oplus$ |
| R_Rs | planet radius in stellar radii | $R_\star$ |
| a_Rs | semi-major axis of th eorbit in stellar radii | $R_\star$ |
| b | impact parameter | $R_\star$ |
| i | orbital inclination with respect to the plane of the sky | degrees |
| Omega  | longitude of the ascending note $\Omega$| degrees |
| lambda | Spin-orbit alignment angle $\lambda$  | degrees |
| delta_occ | Depth of the occulation | Normalized flux |
| phase_amp | Amplitude of the phase curve, *if no occultation was present* | Normalized flux |
| phase_off | Offset of the peak of the phase curve (light  travel time effect not included) | degrees |


```{warning}
`PyORBIT` uses the argument of pericenter of the **planet** $\omega_p$, while other packages may use the argument of pericenter of the **star** $\omega _\star$. Some papers report the latter without specifying the $\star$ pedix.
```

## Keywords

The default keyword is highlighted in boldface.

**orbit**
* accepted values: `circular` | **`keplerian`** | `dynamical`
* define if the planet is moving on a circular orbit ($e=0$, $\omega=90°$), a
  standard Keplerian, or on an orbit computed through N-body integration

**parametrization**
* accepted values: `Standard` | `Standard_Tcent` | `Ford2006` |
  `Ford2006_Tcent` | **`Eastman2013`** |  `Eastman2013_Tcent`
* define the parametrization for eccentricity and argument of periastron:
  ($e$, $\omega_p$) for `Standard`, ($e \sin{\omega_p}$,
  $e \cos{\omega_p}$) for `Ford2006` ([Ford 2006](https://ui.adsabs.harvard.edu/abs/2006ApJ...642..505F/abstract)), ($\sqrt{e} \sin{\omega_p}$,
  $\sqrt{e} \cos{\omega_p}$) for `Eastman2013` ([Eastman et al. 2013](https://ui.adsabs.harvard.edu/abs/2013PASP..125...83E/abstract)).
* Appending `_Tcent` to the
  parametrization label will replace the mean longitude `mean_long` with the
  central time of transit `Tc`.

**use_inclination**
* accepted values:  `True` | **`False`**
* if `True`, the inclination of the planet `i` replaces the impact parameter `b`.

**use_semimajor_axis**
* accepted values: `True` |  **`False`**
* if `True`, the scaled semimajor axis of the planet `a_Rs` replaces the stellar
  density (defined in the `star` section).

**use_time_inferior_conjunction**
* accepted values: `True` |  **`False`** / overridden by `parametrization`
* alternative way to replace the mean longitude `mean_long` with the
  central time of transit `Tc` when set to `True`. It is overridden by `parametrization` if its value ends with `_Tcent`

**use_mass_for_planets**
* accepted values:  `True` |  **`False`**
* if `False`, the mass of the planet replaces the radial velocity semiamplitude.
  It should be activated with those methods that allow the determination of the
  true mass of a planet, e.g., TTVs


## Examples

The most basic example is provided by a non-transiting planet in a circular orbit. For purely coding reasons, the dictionary of a planet (in this case, `b`) cannot be empty, so we provide just one of the keywords:

```yaml
common:
  planets:
    b:
      orbit: circular
```

The code will automatically pick the parametrization to be employed (here, the *mean longitude* as the planet is not transiting), the way the parameter space is explored (here, Base-2 Logarithm for period and RV semi-amplitude), the boundaries for the parameters (in sampler space):

```text
----- common model:  b
mean_long     id:   0  s:Linear      b:[      0.0000,     360.0000]   p:Uniform   []
P             id:   1  s:Log_Base2   b:[     -1.3219,      16.6096]   p:Uniform   []
K             id:   2  s:Log_Base2   b:[     -9.9658,      10.9658]   p:Uniform   []
omega         derived (no id, space, bound)                           p:None   []
e             derived (no id, space, bound)                           p:None   []
```

The default behavior of `PyORBIT` is to set boundaries and parametrizations according to hard-coded when not explicitly provided in the configuration files. While in most cases the default values are fine (e.g., angles, offsets, jitter parameters), sometimes it may be useful to restrict the range explored by the parameters:


```yaml
common:
  planets:
    b:
      orbit: circular
      boundaries:
        P: [0.50, 500.0]
        K: [0.01, 300.0]
```

```{warning}
Boundaries and priors must be expressed in the *physical*/*natural* space, even if the parameters space will be explored in logarithmic space. If the parameters are explored in logarithmic space, the boundaries must be strictly positive.
```

The logarithmic parametrization allows an efficient exploration of the parameter space across several orders of magnitude. When the parameter space is restricted by tight boundaries or a prior is provided, as in the case of transiting planets (i.e., well-known period and RV semi-amplitude within an expected range), then the space exploration can be switched to Linear (if it is not already in that space):

```yaml
common:
  planets:
    b:
      orbit: circular
      use_time_inferior_conjunction: True
      boundaries:
        P: [2.20, 2.25]
        Tc: [2456194.00, 2456194.1.0]
        K: [190.0, 220.0]
      priors:
        P: ['Gaussian', 2.218574944, 0.000000030]
        Tc: ['Gaussian', 2456194.067619, 0.000034]
      spaces:
        P: Linear
        K: Linear
```

Finally, keep in mind that spaces, boundaries, and priors may change depending on the datasets you are fitting. In the case of transit photometry, for example, you may want to specify boundaries for the scaled planetary radius, while dropping the boundaries for the RV semi-amplitude. 

```yaml
common:
  planets:
    b:
      orbit: keplerian
      parametrization: Eastman2013
      use_time_inferior_conjunction: True
      use_inclination: False    # can be omitted when default value is used
      use_semimajor_axis: False # can be omitted when default value is used
      boundaries:
        P: [2.20, 2.25]
        Tc: [2456194.00, 2456194.1.0]
        R_Rs: [0.00, 1.00]
        e: [0.00, 0.90]
      priors:
        P: ['Gaussian', 2.218574944, 0.000000030]
        Tc: ['Gaussian', 2456194.067619, 0.000034]
      spaces:
        P: Linear
        K: Linear
```

In this case, we also required to use a keplerian orbit with the [Eastman et al. 2013](https://ui.adsabs.harvard.edu/abs/2013PASP..125...83E/abstract) parametrization.

<!---
PyORBIT has been used in the following works
      use_inclination: False    # can be omitted when default value is used
      use_semimajor_axis: False # can be omitted when default value is used
        P: [1.210, 1.240]
        Tc: [59144.60, 59144.63]
        b: [0.0, 1.0]
        R_Rs: [0.00, 1.00]

--->