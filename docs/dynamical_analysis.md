(dynamical_analysis)=

# Dynamical Analysis 

`PyORBIT` integrates `TRADES` (TRAnsits and Dynamics of Exoplanetary Systems), a Fortran 90 code developed by [Borsato et al. 2014](https://ui.adsabs.harvard.edu/abs/2014A%2526A...571A..38B/abstract) to perform dynamical modelling of exoplanetary signals. The code has been expanded with a Python interface ([Borsato et al. 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.484.3233B/abstract)) and to perform photodynamical modelling ([Borsato et al. 2024](https://ui.adsabs.harvard.edu/abs/2024A%2526A...689A..52B/abstract)).

By using `PyORBIT` and `TRADES` together, you can model the TTVs or the transit light curves of interacting planets while account for stellar activity, for example with  


```{toctree}
:maxdepth: 1
dynamical_analysis/pyorbit_trades_integration.md
dynamical_analysis/TTV_and_RV_fit.md
```
