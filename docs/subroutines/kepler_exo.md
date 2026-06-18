# Keplerian orbital subroutines

The module `pyorbit.subroutines.kepler_exo` contains low-level utilities for
Keplerian orbits.  They are used by the radial-velocity, transit and derived
parameter models to convert between orbital parametrizations, compute radial
velocities, and derive masses or orbital positions.

Unless stated otherwise, periods and times are in days, angles supplied by the
user are in degrees, radial velocities are in m/s, and masses are in Solar
mass units.  The functions accept scalar values and, where used internally,
NumPy arrays.

## Backward-compatible aliases

The module still exports several legacy function names.  They call the newer
`kepler_compute_*` functions and are kept for backward compatibility.

| Alias | Current function |
| --- | --- |
| `kepler_K1` | `kepler_compute_rv_semiamplitude` |
| `kepler_RV` | `kepler_compute_rv` |
| `kepler_RV_T0P` | `kepler_compute_rv_deltabjd` |
| `kepler_phase2Tc_Tref` | `kepler_compute_deltaTc_from_meanlong` |
| `kepler_Tc2phase_Tref` | `kepler_compute_meanlong_from_deltaTc` |
| `kepler_Tc2Tperi_Tref` | `kepler_compute_deltaTperi_from_deltaTc` |
| `kepler_phase2Tperi_Tref` | `kepler_compute_deltaTperi_from_meanlong` |
| `get_planet_mass` | `kepler_get_planet_mass` |
| `kepler_true_anomaly_orbital_distance` | `kepler_compute_trueanomaly_orbitaldistance` |

New code should prefer the current function names.

## Solving Kepler's equation

`kepler_E(M_in, ecc)` solves Kepler's equation and returns the eccentric
anomaly `E` from the mean anomaly `M_in` and eccentricity `ecc`.

The solver wraps the input mean anomaly to the `[0, 2 pi)` interval and uses an
iterative fourth-order correction scheme.  This function is the internal
building block used by the radial-velocity and orbital-distance routines for
eccentric orbits.

The helper functions `f0_keplerE` and `f1_keplerE` implement the residual and
first derivative of Kepler's equation.  They are support routines and are not
part of the preferred public interface.

## Radial-velocity calculations

`kepler_compute_rv_semiamplitude(mass_primary, mass_secondary, period,
inclination, eccentricity)` returns the radial-velocity semi-amplitude of the
primary star.  The primary and secondary masses are in Solar masses, the period
is in days, the inclination is in degrees, and the result is in m/s.

`kepler_compute_rv(bjd, Tperi, period, rv_semiamplitude, ecc, omega_deg)`
computes the Keplerian radial velocity at absolute times `bjd` when the time of
periastron passage `Tperi` is known.

`kepler_compute_rv_deltabjd(bjd_tref, rv_semiamplitude, period, mean_long, ecc,
omega_deg, Omega_deg=0.0)` computes the same radial velocity using times
relative to the reference epoch, `bjd_tref = BJD - Tref`, and the mean longitude
at the reference epoch.

For nearly circular orbits, `abs(ecc) < 1e-3`, the routines use the mean anomaly
directly as the true anomaly.  Negative eccentricities are interpreted by
flipping the sign and adding 180 degrees to the argument of periastron.

## Converting between `Tc`, mean longitude and periastron time

These routines convert between the orbital phase conventions used in PyORBIT:

| Function | Output |
| --- | --- |
| `kepler_compute_deltaTc_from_meanlong(period, mean_long, ecc=0., omega_deg=90., Omega_deg=0.0)` | Minimum positive time difference between inferior conjunction and the reference epoch. |
| `kepler_compute_meanlong_from_deltaTc(period, delta_Tc, ecc, omega_deg, Omega_deg=0.0)` | Mean longitude at the reference epoch, in degrees in the range `[0, 360)`. |
| `kepler_compute_deltaTperi_from_deltaTc(period, delta_Tc, ecc, omega_deg)` | Time difference between periastron passage and the reference epoch, starting from `delta_Tc`. |
| `kepler_compute_deltaTperi_from_meanlong(period, mean_long, ecc, omega_deg, Omega_deg=0.0)` | Time difference between periastron passage and the reference epoch, starting from mean longitude. |

The inferior conjunction is defined as the instant when the planet is closest to
the line of sight between the observer and the star.  For circular orbits this
corresponds to the central transit time.

## Planet mass from radial velocity

`kepler_get_planet_mass(period, rv_semiamplitude, ecc, mass_star,
approximation_limit=30.)` estimates the planet mass from the period, RV
semi-amplitude, eccentricity and stellar mass.

The result is returned in Solar mass units.  For low-mass planets the function
uses the approximation `M_planet << M_star`.  If the approximate mass is larger
than `approximation_limit` Earth masses, the routine switches to a numerical
solution of the full mass equation.

The support routines are:

| Function | Role |
| --- | --- |
| `get_approximate_mass` | Approximate solution in Solar masses under `M_planet << M_star`. |
| `f_get_mass` | Residual used by the numerical root finder. |

To convert the returned Solar-mass value:

```python
import pyorbit.subroutines.constants as constants

mass_planet_earth = mass_planet_solar * constants.Msear
mass_planet_jupiter = mass_planet_solar * constants.Msjup
```

## True anomaly and orbital distance

`kepler_compute_trueanomaly_orbitaldistance(bjd_tref, semimajor_axis,
delta_Tc, period, ecc, omega_deg, Omega_deg=0.0)` returns the true anomaly and
the orbital distance at the requested time relative to the reference epoch.

The orbital distance is returned in the same units as `semimajor_axis`, so the
function can be used both with scaled semi-major axes and physical distances.
For circular orbits the orbital distance is equal to `semimajor_axis`.

## Example

```python
from pyorbit.subroutines.kepler_exo import (
    kepler_compute_deltaTc_from_meanlong,
    kepler_compute_rv_deltabjd,
)

period = 10.0
mean_long = 45.0
ecc = 0.1
omega = 90.0

delta_tc = kepler_compute_deltaTc_from_meanlong(period, mean_long, ecc, omega)
rv = kepler_compute_rv_deltabjd(
    bjd_tref=0.0,
    rv_semiamplitude=5.0,
    period=period,
    mean_long=mean_long,
    ecc=ecc,
    omega_deg=omega,
)
```

These routines are intentionally small and procedural.  In normal PyORBIT
configuration files they are used through the higher-level models rather than
called directly.
