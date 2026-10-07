(planets)=

# Planets

The `planets` container holds one common object for each planet. The object
name (for example, `b` or `c`) is the label used by the models that include
that planet. Only parameters required by the selected models need to be
specified; other bounds, priors, spaces, and fixed values come from
`CommonPlanets`.

The reference epoch $T_{\rm ref}$ is set in the top-level `parameters`
section. Orbital elements such as `mean_long` refer to that epoch. The mean
longitude $L_0$ is $\omega_p + \Omega + M_0$, where $M_0$ is the mean anomaly
at $T_{\rm ref}$. For radial velocities alone, the ascending node normally
cannot be measured and its default value is fixed at $180^{\circ}$.

## Model definition and requirements

- container name: `planets`
- common object name: the planet label, such as `b`
- source class: `CommonPlanets`
- requirements: declare the planet under `common: planets:` and refer to its
  label in the `planets` list of each model that uses it. Transit geometry can
  also require stellar density or stellar mass and radius from `star_parameters`.

```yaml
common:
  planets:
    b:
      orbit: circular
models:
  radial_velocities:
    planets: [b]
```

## Model parameters

The parameters available to a planet depend on the selected orbital and
photometric models. Orbital and geometric quantities are:

| Name | Symbol | Parameter | Unit |
| :--- | :----- | :-------- | :--- |
| `P` | $P$ | Orbital period | days |
| `K` | $K$ | Radial-velocity semiamplitude | m/s |
| `Tc` | $T_c$ | Time of inferior conjunction | days |
| `mean_long` | $L_0$ | Mean longitude at $T_{\rm ref}$ | degrees |
| `e` | $e$ | Orbital eccentricity | unitless |
| `omega` | $\omega_p$ | Argument of periastron of the planet | degrees |
| `e_coso` | $e\cos\omega_p$ | [Ford 2006](https://ui.adsabs.harvard.edu/abs/2006ApJ...642..505F/abstract) eccentricity component | unitless |
| `e_sino` | $e\sin\omega_p$ | Ford 2006 eccentricity component | unitless |
| `sre_coso` | $\sqrt{e}\cos\omega_p$ | [Eastman et al. 2013](https://ui.adsabs.harvard.edu/abs/2013PASP..125...83E/abstract) eccentricity component | unitless |
| `sre_sino` | $\sqrt{e}\sin\omega_p$ | Eastman et al. 2013 eccentricity component | unitless |
| `M_Me` | $M_{\rm p}$ | Planet mass | Earth masses |
| `R_Rs` | $R_{\rm p}/R_\star$ | Planet-to-star radius ratio | unitless |
| `a_Rs` | $a/R_\star$ | Scaled semimajor axis | unitless |
| `b` | $b$ | Transit impact parameter | unitless |
| `i` | $i$ | Orbital inclination to the plane of the sky | degrees |
| `Omega` | $\Omega$ | Longitude of the ascending node | degrees |
| `lambda` | $\lambda$ | Sky-projected spin-orbit angle | degrees |
| `delta_occ` | — | Occultation depth | normalized flux |
| `phase_amp` | — | Phase-curve amplitude without occultation | normalized flux |
| `phase_off` | — | Phase-curve peak offset, excluding light-travel time | degrees |

```{warning}
`PyORBIT` uses the argument of periastron of the **planet**, $\omega_p$.
Other packages and papers may use the stellar argument $\omega_\star$
without specifying the subscript.
```

## Keywords

The default value is highlighted in boldface.

**orbit**
* accepted values: `circular` | **`keplerian`** | `dynamical` | `apodized`
* selects a circular orbit ($e=0$, $\omega_p=90^{\circ}$), a Keplerian orbit,
  an orbit computed by N-body integration, or an apodized Keplerian signal.

**parametrization**
* accepted values: `Standard` | `Ford2006` | **`Eastman2013`**; each also
  accepts a `_Tcent` or `_Tc` suffix.
* selects ($e$, $\omega_p$) for `Standard`, ($e\cos\omega_p$,
  $e\sin\omega_p$) for `Ford2006`, or ($\sqrt{e}\cos\omega_p$,
  $\sqrt{e}\sin\omega_p$) for `Eastman2013`. A `_Tcent` or `_Tc` suffix
  selects `Tc` in place of `mean_long`.

**use_inclination**
* accepted values: `True` | **`False`**
* if `True`, samples orbital inclination `i` instead of impact parameter `b`
  in models that require transit geometry.

**use_scaled_semimajor_axis**
* accepted values: `True` | **`False`**
* if `True`, samples `a_Rs` instead of using the stellar density from the
  `star_parameters` object to determine the scaled semimajor axis.

**use_time_inferior_conjunction**
* accepted values: `True` | **`False`**
* if `True`, samples the time of inferior conjunction `Tc` instead of
  `mean_long`; a `parametrization` label ending in `_Tcent` or `_Tc` also
  selects `Tc`.

**use_mass**
* accepted values: `True` | **`False`**
* if `True`, samples the planet mass `M_Me` instead of radial-velocity
  semiamplitude `K`. The inclination or impact parameter and stellar mass
  must then be available; this is useful for astrometry and dynamical fits.

**use_longitude_of_nodes**
* accepted values: `True` | **`False`**
* if `True`, samples `Omega` instead of holding it at its default value.


## Examples

A minimal non-transiting, circular planet can be declared with one keyword:

```yaml
common:
  planets:
    b:
      orbit: circular
```

With an RV model that includes `b`, the default phase parameter is
`mean_long`. `P` and `K` use base-2 logarithmic sampler spaces. The printed
bounds below are therefore in sampler space:

```text
----- common model:  b
mean_long     id:   0  s:Linear      b:[      0.0000,     360.0000]   p:Uniform   []
P             id:   1  s:Log_Base2   b:[     -1.3219,      16.6096]   p:Uniform   []
K             id:   2  s:Log_Base2   b:[     -9.9658,      10.9658]   p:Uniform   []
omega         derived (no id, space, bound)                           p:None   []
e             derived (no id, space, bound)                           p:None   []
```

The circular RV fit in
[HD189733_example01_onedataset.yaml](https://github.com/LucaMalavolta/PyORBIT_examples/blob/main/quickstart/HD189733_example01_onedataset.yaml)
restricts the physical period and semiamplitude bounds:


```yaml
common:
  planets:
    b:
      orbit: circular
      boundaries:
        P: [0.50, 5.0]
        K: [0.01, 300.0]
```

```{warning}
Express bounds and priors in physical units even when sampling in logarithmic
space. Bounds for logarithmic parameters must be strictly positive.
```

For a transiting planet, the fit in
[HD189733_example05_TESSandRV.yaml](https://github.com/LucaMalavolta/PyORBIT_examples/blob/main/quickstart/HD189733_example05_TESSandRV.yaml)
uses the transit time `Tc` as its phase parameter and switches `P` to a
linear sampler space over a narrow range:

```yaml
common:
  planets:
    b:
      orbit: circular
      use_time_inferior_conjunction: True
      boundaries:
        P: [2.2185600, 2.2185800]
        Tc: [2459770.4100, 2459770.4110]
        K: [0.01, 300.0]
      spaces:
        P: Linear
```

The combined multiband transit and RV fit in
[HD189733_example06_multiband_RVs.yaml](https://github.com/LucaMalavolta/PyORBIT_examples/blob/main/quickstart/HD189733_example06_multiband_RVs.yaml)
instead uses a Keplerian orbit with the `Eastman2013` eccentricity
parametrization:

```yaml
common:
  planets:
    b:
      orbit: keplerian
      parametrization: Eastman2013
      use_time_inferior_conjunction: True
      boundaries:
        P: [2.2185600, 2.2185800]
        Tc: [2459770.4100, 2459770.4110]
        K: [0.01, 300.0]
        e: [0.00, 0.95]
      spaces:
        P: Linear
```

For astrometry and RVs, the
[HD5388_test03.yaml](https://github.com/LucaMalavolta/PyORBIT_examples/blob/main/astrometry/HD5388_test03.yaml)
example samples the planet's true mass and orbital inclination:

```yaml
common:
  planets:
    b:
      orbit: keplerian
      use_mass: True
      use_inclination: True
      boundaries:
        P: [100.0, 10000.0]
        M_Me: [100, 3000.0]
        e: [0.00, 0.90]
        i: [0.0, 90.0]
```
