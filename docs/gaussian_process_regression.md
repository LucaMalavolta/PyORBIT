(gaussian_process_regression)=

# Gaussian process regression

Gaussian process regression (GPR) is a nonparametric, Bayesian approach to regression that has been very successful for the analysis of radial velocity datasets in the presence of stellar activity, e.g., [Haywood at al. 2014](https://ui.adsabs.harvard.edu/abs/2014MNRAS.443.2517H/abstract),  [Grunblatt et al. 2015](https://ui.adsabs.harvard.edu/abs/2015ApJ...808..127G/abstract). The first stable and tested implementation of GPR in `PyORBIT` dates back to 2018. 

In this section, we will discuss only models encompassing independent covariance matrix among datasets, with some hyperparameters in common. I will refer to these approaches as *classic* or *standard* Gaussian processes.

```{toctree}
:maxdepth: 1
gaussian_process/quasiperiodic_kernel.md
```
