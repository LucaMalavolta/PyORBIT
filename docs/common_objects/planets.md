(planets)=

# Planets

### Parameters definition

Each *planet* is defined by the following parameters

| Name        | Parameter | Unit     |
| :---        | :----     | :---     |
| P      | Orbital period of the planet                      | days     |
| K      | Radial velocity semiamplitude        | m/s      |
| Tc     | Central time of transit              | days      |
| mean_long | Mean longitude of the orbit, assuming $\Omega=0$ | degrees |
| e      | eccentricity of the orbit | adimensional |
| omega  | argument of periastron of the *planet* $\omega_p$  | adimensional |
| e_coso | $ e \cos{\omega}$ in Ford 2006 parametrization | adimensional |
| e_sino | $ e \sin{\omega}$ in Ford 2006 parametrization | adimensional |
| sre_coso | $ \sqrt{e} \cos{\omega}$ in Eastman 2013 parametrization | adimensional |
| sre_sino | $ \sqrt{e} \sin{\omega}$ in Eastman 2013 parametrization | adimensional |
| M_Me | planet mass in Earth masses | $M_\oplus$ |
| R_Rs | planet radius in stellar radii | $R_\star$ |
| a_Rs | semi-major axis of th eorbit in stellar radii | $R_\star$ |
| b | impact parameter | $R_\star$ |
| i | orbital inclination with respect to the plane of the sky | degrees |
| Omega  | longitude of the ascending note $\Omega$| degrees |
| lambda | $\lambda$ | degrees |
| phase_amp | **check** Amplitude of the phase curve| |
| phase_off | **check** Offset of the peack of the phase curve| degrees |
| redist | **check**  | adimensional |
| insol | **check**  | adimensional |


```{warning}
`PyORBIT` uses the argument of pericenter of the **planet** $\omega_p$, while other packages may use the argument of pericenter of the **star** $\omega _\star$. Some papers report the former without specifying the $\star$ pedix.
```

### Keywords

**orbit**
* accepted values: `circular` | `keplerian` | `dynamical`. default: `keplerian`
* define if the planet is moving on a circular orbit ($e=0$, $\omega=90°$), a
  standard Keplerian, or on an orbit computed through N-body integration

**parametrization**
* accepted values: `Standard` | `Standard_Tcent` | `Ford2006` |
  `Ford2006_Tcent` | `Eastman2013` |  `Eastman2013_Tcent`. default: `Eastman2013`
* define the parametrization for eccentricity and argument of periastron:
  ($e$, $\omega$) for `Standard`, ($e \sin{\omega}$,
  $e \cos{\omega}$) for `Ford2006`, ($\sqrt{e} \sin{\omega}$,
  $\sqrt{e} \cos{\omega}$) for `Eastman2013`.
* Appending `_Tcent` to the
  parametrization label will replace the mean longitude `mean_long` with the
  central time of transit `Tc`.

**use_inclination**
* accepted values:  `True` | `False`. default: `False`
* if `True`, the inclination of the planet `i` replaces the impact parameter `b`.

**use_semimajor_axis**
* accepted values: `True` | `False`. default: `False`
* if `True`, the scaled semimajor axis of the planet `a_Rs` replaces the stellar
  density (defined in the `star` section).

**use_time_of_transit**
* accepted values: `True` | `False`. default: `False` / overridden by `parametrization`
* alternative way to replace the mean longitude `mean_long` with the
  central time of transit `Tc` when set to `True`. It is overridden by `parametrization` if its
  value is ending with `_Tcent`

**use_mass_for_planets**
* accepted values:  `True` | `False`. default: `False`
* If `False`, the mass of the planet replaces the radial velocity semiamplitude.
  It should be activated with those methods that allow the determination of the
  true mass of a planet, e.g., TTVs

Example for a transiting planet:

```yaml
common:
  planets:
    b:
      orbit: circular
      parametrization: Eastman2013_Tcent
      use_inclination: False    # can be omitted when default value is used
      use_semimajor_axis: False # can be omitted when default value is used
      boundaries:
        P: [1.210, 1.240]
        Tc: [59144.60, 59144.63]
        b: [0.0, 1.0]
        R_Rs: [0.00, 1.00]
```