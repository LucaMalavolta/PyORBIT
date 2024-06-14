(gaussian_process_regression)=

# Gaussian Process regression

Gaussian Process regression (GPR) is a nonparametric, Bayesian approach to regression that has been very successful in the analysis of radial velocity datasets in the presence of stellar activity, e.g., [Haywood et al. 2014](https://ui.adsabs.harvard.edu/abs/2014MNRAS.443.2517H/abstract),  [Grunblatt et al. 2015](https://ui.adsabs.harvard.edu/abs/2015ApJ...808..127G/abstract). The first stable and tested implementation of GPR in `PyORBIT` dates back to 2018 and was showcased in [Malavolta et al. 2018](https://ui.adsabs.harvard.edu/abs/2018AJ....155..107M/abstract). Over the years, several kernels and new packages have been implemented.

In the following, I will assume that you are already familiar with the mathematical basis of GP. The [Gaussian Process website](https://gaussianprocess.org/) and the review by [Aigrain & Foreman-Mackey 2022](https://ui.adsabs.harvard.edu/abs/2023ARA%26A..61..329A/abstract) represent good starting points to learn more about GPs.

In this section, we will discuss only models encompassing independent covariance matrix among datasets, with some hyperparameters in common. I will refer to these approaches as *classic* or *standard* Gaussian processes regression. Multi-dimensional GP, formerly known as GP framework, will be presented in a dedicated section.

In `PyORBIT` unlimited number of additional datasets can be included for the simultaneous training of the hyperparameters. Unless differently specified, all the hyperparameters will be shared (if referred to the same *model* and *common model*) except the amplitude of covariance matrix, which is dataset-dependent. Each dataset will be characterized by its own covariance matrix.



```{toctree}
:maxdepth: 1
gaussian_process/quasiperiodic_kernel.md
gaussian_process/quasiperiodic_cosine_kernel.md
gaussian_process/multidimensional_gaussian_processes
```
