(planets)=

# Planets

### Parameters definition

Each *planet* is defined by the following parameters

| Name        | Parameter | Unit     |
| :---        | :----     | :---     |
| P      | Orbital planets                      | days     |
| K      | Radial velocity semiamplitude        | m/s      |
| Tc     | Central time of transit              | days      |
| mean_long | Mean longitude of the orbit, assuming $\Omega=0$ | degrees |
| e      | eccentricity of the orbit | adimensional |
| omega  | argument of pericenter of the *planet* $\omega_p$  | adimensional |
| e_coso | $ e \cos{\omega}$ in Ford 2006 parametrization | adimensional |
| e_sino | $ e \sin{\omega}$ in Ford 2006 parametrization | adimensional |
| sre_coso | $ \sqrt{e} \cos{\omega}$ in Eastman 2013 parametrization | adimensional |
| sre_sino | $ \sqrt{e} \sin{\omega}$ in Eastman 2013 parametrization | adimensional |
| M_Me | planet mass in Earth masses | $M_\bigoplus$ |
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
``PyORBIT`` uses the argument of pericenter of the **planet** $\omega_p$, while other packages may use the argument of pericenter of the **star** $\omega _\star$, but some papers report the second without specifying the $\star$ pedix. 
```

### Keywords

orbit | ``circular`` | ``keplerian`` | ``dynamical`` |
: default: ``keplerian``
: define if the planet is moving on a circular orbit ($e=0$, $\omega=90°$), on a
computed through N-body integration

parametrization | ``Standard`` | ``Standard_Tcent`` | ``Ford2006`` | ``Ford2006_Tcent`` |
: default: ``Eastman2013``
: define something

use_inclination | ``True`` | ``False`` |
: default: ``False``
: define something and blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah 

use_semimajor_axis | ``True`` | ``False`` |
: default: ``False``
: define something

use_time_of_transit | ``True`` | ``False`` |
: default: ``False``
: define something

use_mass_for_planets | ``True`` | ``False`` |
: default: ``False``
: define something

