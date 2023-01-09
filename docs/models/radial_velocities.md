(radial_velocities)=

# Planetary RV signal

- model name: ``radial_velocities``
- required common objects : ``planets`` , ``star``

### Model parameters 

| Name        | Parameter | Common?  | Definition  | Notes |
| :---        | :-------- | :-------------  | :-----  | :---- |
| P      | orbital period of the planet | common | ``planets``     | |
| K      | Radial velocity semiamplitude | common | ``planets``     | |
| mean_long | mean longitude of the orbit, assuming $\Omega=0$ | common | ``planets`` | (1) |
| Tc     | Central time of transit              | common | ``planets``     | (2)|
| e      | eccentricity of the orbit  | common | ``planets`` | (3) |
| omega  | argument of pericenter of the *planet* $\omega_p$  | common |  ``planets`` | (3) |
| e_coso | $ e \cos{\omega}$ in Ford 2006 parametrization | common |  ``planets`` | (4) |
| e_sino | $ e \sin{\omega}$ in Ford 2006 parametrization | common |  ``planets`` | (4) |
| sre_coso | $ \sqrt{e} \cos{\omega}$ in Eastman 2013 parametrization | common |  ``planets`` | (5)|
| sre_sino | $ \sqrt{e} \sin{\omega}$ in Eastman 2013 parametrization | common |  ``planets`` | (5) |

Notes:

  1. replaced by ``Tc`` when ``_Tcent`` parametrization is used
  2. ``_Tcent`` parametrization only
  3. ``Standard`` parametrization only
  4. ``Ford2006`` parametrization only
  5. ``Eastman2013`` parametrization only

### Keywords

to be defined