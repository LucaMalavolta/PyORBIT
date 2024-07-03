(home)=

# PyORBIT
## a code for exoplanet orbital parameters and stellar activity.
<!---
## the ultimate tool for exoplanet characterization
--->

[![GitHub][github-badge]][github-link]
[![PyPI version fury.io](https://badge.fury.io/py/pyorbit-package.svg)](https://pypi.python.org/pypi/pyprbot-package/)


``PyORBIT`` is a robust, versatile framework for the characterization of planetary systems.
With ``PyORBIT`` you can model light curves, radial velocities, activity indexes, and transit time variations.
In addition to the exoplanet signature, you can add to your model instrumental systematics and stellar activity,
including Gaussian processes regression with a variety of kernels. Both Markov Chain Monte-Carlo samplers and Nested Sampling algorithms are supported.
Different parametrization can be used for orbital parameters, and orbits can be computed using dynamical integration or just non-interacting Keplerians.
For every parameter in the model, it is possible to set a prior, explore it in linear or logarithmic space,
or keep it fixed if necessary. Thanks to abstraction, it is virtually possible to implement any physical model you can think of.

Every information you need to analyze your data is contained in an easy-to-use and portable configuration file in the ``yaml`` format.
Alternatively, for easy automatization, PyORBIT can be called as a Python function by passing a dictionary instead of a configuration file.

<!---
..
  One of the main strength of PyORBIT is in the abstraction of the models. For
  example,  you can use more than one Gaussian process, with independent
  hyper-parameters, just by using different labels for the abstract class.
  Alternatively, for your Gaussian process regression you can use "celerite" on
  your light curve and "george " on your radial velocities, by having them sharing
  the same hyper-parameters.
-->

[**Check the poster**](https://k-poster.kuoni-congress.info/eas-2024/poster/bcefa082-f961-4539-a5f8-3b13783b520c) presented at the European Astronomical Society Annual Meeting in Padova

## Updates on version 10

```{warning}
Starting from version 10.3, `PyORBIT` has been upgraded to support `tinygp` (version 0.3.0), which in turns requires Python **3.10** to work properly.
If you are using `PyORBIT` \=> 10.3, follow the installation instructions to create an environment with Python 3.10
```

- *Improved speed*
After several failed attempts, I finally managed to apply the advice from the [emcee parallelization page](https://emcee.readthedocs.io/en/stable/tutorials/parallel/) to the rather complex structure of PyORBIT. The speed-up is noticeable for large datasets (e.g., photometry).

- *Rossiter McLaughlin*
Rossiter McLaughlin effect can now be precisely modelled using the CCF simulation approach employed in [Covino et al. 2013](https://ui.adsabs.harvard.edu/abs/2013A%26A...554A..28C/abstract).
When the rotation period is known - together with the stellar radius - the stellar inclination can be derived avoiding the bias reported by [Masuda & Winn 2020](https://ui.adsabs.harvard.edu/abs/2020AJ....159...81M/abstract).
This model has been successfully employed in [Mantovan et al. 2024b](https://ui.adsabs.harvard.edu/abs/2024A%26A...684L..17M/abstract) for the characterization of TOI-5398b.

- *Multidimensional Gaussian Process*
The model was introduced a few years back, but now it can finally take advantage of improved parallelization.
Recent examples of multidimensional Gaussian Processes through `PyORBIT` can be found in [Nardiello et al. 2023](https://ui.adsabs.harvard.edu/abs/2022A%26A...664A.163N/abstract) and [Mantovan et al. 2024a](https://ui.adsabs.harvard.edu/abs/2024A%26A...682A.129M/abstract).

- *tinyGP for better performances*
Working on both *classic* and _multidimensional_ Gaussian Processes, although the former is showing some JAX problems when producing the output results.

```{admonition} No back-compatibility
Version 10 is not compatible with the results obtained with version 9.
If you have been using the development version of V10, you may run into incompatibility issues as well.
```

```{warning}
In the previous version of the documentation and in several papers, it was reported that `PyORBIT` was relying on the quasi-periodic kernel definiton by [Grunblatt et al. 2015](https://ui.adsabs.harvard.edu/abs/2015ApJ...808..127G/abstract). I came to realize only recently that the factor $2$ accompanying the decay time-scale of the activity regions $P_\mathrm{dec}$ was implicitely included in the exponenitial-squared kernel as defined in `george` and `tinygp`, making the kernel definiton identical to the one reported in [Rajpaul et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.452.2269R/abstract). Please keep this problem in mind when comparing your results with other analysis.
```

# Documentation updates

Documentation is being updated slowly but steadily, with new pages appearing every week. Below you can check which pages have been updated recently (last update: 28/07/2023).

- *Quickstart*
  + Radial velocities
  + Lightcurve analysis
  + Combined fit
- *Data modeling*
  + Planetary RV signal
  + Correlated datasets
  + Gaussian processes regression **in progress**
  + Multidimensional GPs **new**
  + TinyGP caveats **new**
  + CHEOPS detrending **coming soon**
  + Polynomial **coming soon**

```{tip}
The configuration files used in this documentation can be found at [this repository](https://github.com/LucaMalavolta/PyORBIT_examples)
```


## References

If you find `PyORBIT` useful for your work, please cite [Malavolta et al. (2016)](https://ui.adsabs.harvard.edu//#abs/2016A&A...588A.118M/abstract) and [Malavolta et al. (2018)](https://ui.adsabs.harvard.edu//#abs/2018AJ....155..107M/abstract).

<!---
PyORBIT has been used in the following works
--->


## Some background

``PyORBIT`` started in 2015 as an exercise to learn Python (2), and also because at
the time the only publicly available code for the Bayesian analysis of RVs was
written in IDL, for which I didn't have a license. Since then I've been adding new options every time I needed, and I kept updating the code while improving my Python skills.

If you are wondering what ``PyORBIT`` stands for: I am bad at acronym creation so
I decided to write it with capital _ORBIT_ just because I liked how it looked.
Feel free to submit your retrofitting acronym!

``PyORBIT`` has been now converted and tested for Python 3 for a while, back compatibility with Python 2 is not guaranteed anymore at least since version 8.0. If you haven't done it already, I strongly suggest you move to Python 3.

[github-badge]: https://img.shields.io/badge/GitHub-PyORBIT-blue
[github-link]: https://github.com/LucaMalavolta/PyORBIT


## Table of contents

```{toctree}
:maxdepth: 2
installation
quickstart
prepare_datasets
prepare_yaml
common_objects
models
gaussian_process_regression
multidimensional_gps
running_pyorbit
GitHub Repository <https://github.com/LucaMalavolta/PyORBIT>
```

<!---
```{eval-rst}
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   prepare_datasets
   prepare_yaml
   api
```
--->

<!---

..
  Following `PEP 8 Style Guide for Python Code <https://www.python.org/dev/peps/pep-0008/>`_  ,
  `PEP 257 Docstring Conventions <https://www.python.org/dev/peps/pep-0257/>`_ and `Google Python Style Guide <http://google.github.io/styleguide/pyguide.html>`_

..
  Indices and tables
  ==================

  * :ref:`genindex`
  * :ref:`modindex`
  * :ref:`search`

--->

[![forthebadge made-with-python](http://ForTheBadge.com/images/badges/made-with-python.svg)](https://www.python.org/)

This documentation has been rendered using the [Sphinx Book Theme](https://sphinx-book-theme.readthedocs.io/) and the [Myst parser](https://myst-parser.readthedocs.io/)
