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


[**Check the poster**](https://k-poster.kuoni-congress.info/eas-2024/poster/bcefa082-f961-4539-a5f8-3b13783b520c) presented at the European Astronomical Society Annual Meeting in Padova


## Documentation update

I have extensively updated the documentation, and several examples are now available in a new dedicated repository. Check them out:
- [`PyORBIT` 11 Documentation](https://pyorbit.readthedocs.io/)
- [`PyORBIT` examples repository](https://github.com/LucaMalavolta/PyORBIT_examples)

For nostalgic people, `PyORBIT` 10 is available as a downloadable package from the [`PyORBIT` Releases page](https://github.com/LucaMalavolta/PyORBIT/releases). 

## Updates on version 11

- The mean longitude now strictly follows its definition as the sum of the argument of pericenter $\omega$, the longitude of the nodes $\Omega$, and the mean anomaly $\mathcal{M}_0$ at the reference time $T_{\rm ref}$. The default value for the longitude of the nodes $\Omega$ has been changed to 180 degrees. All the orbital parameters refer to the *planet*. The code should be able to recognise analysis performed with an older version of PyORBIT and use the corresponding default value of $\Omega$

- Possibility to perform  simultaneous photometric and spectroscopic multivariate Gaussian Process regression, while selecting the hyperparameters to be shared

- New kernels for stellar/photometric variability in photometric data


:::{admonition} Keep `PyORBIT` updated!
:class: tip

Current version of `PyORBIT` is [![PyPI version fury.io](https://badge.fury.io/py/pyorbit-package.svg)](https://pypi.python.org/pypi/#pyprbot-package/), keep the software updated for improved performances and new models by running `pip install --upgrade pyorbit-package` on your terminal
:::

## Relevant updates in previous versions 

**Additional BIC/AIC/AICc information** in output and in dictionaries when running `pyorbit_results`. Information criteria computed through the ln-posterior (in addition to the ln-likelihood) have been dropped, as definitely wrong. There is also a small bugfix regarding datasets not in chronological order which was preventing analysis with `spleaf_esp` model from starting altogether.

**S+LEAF exponential-sine periodic kernel now supported**: The exponential-sine periodic (ESP) kernel is a fast approximation of the quasi-periodic (QP) kernel, implemented in the `S+LEAF` software by [Delisle et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...638A..95D/abstract) and [Delisle et al. 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...659A.182D/abstract). The kernel has been implemented and tested in `PyORBIT`, and it can be used as a faster replacement of `tinyGP`

**Outputs available as python dictionaries**: textual output is now also saved as dictionaries in the corresponding *dictionaries* folder.

**Updated TTV measurement and harmonic models** TTV measurement models have all been revised and improved to support a variety of cases. 

**Improved speed**: after several failed attempts, I finally managed to apply the advice from the [emcee parallelization page](https://emcee.readthedocs.io/en/stable/tutorials/parallel/) to the rather complex structure of PyORBIT. The speed-up is noticeable for large datasets (e.g., photometry).

**Rossiter McLaughlin effect** can now be precisely modelled using the CCF simulation approach employed in [Covino et al. 2013](https://ui.adsabs.harvard.edu/abs/2013A%26A...554A..28C/abstract).
When the rotation period is known - together with the stellar radius - the stellar inclination can be derived avoiding the bias reported by [Masuda & Winn 2020](https://ui.adsabs.harvard.edu/abs/2020AJ....159...81M/abstract).
This model has been successfully employed in [Mantovan et al. 2024b](https://ui.adsabs.harvard.edu/abs/2024A%26A...684L..17M/abstract) for the characterization of TOI-5398b.

**Multidimensional Gaussian Process** The model was introduced a few years back, but now it can finally take advantage of improved parallelization. Recent examples of multidimensional Gaussian Processes through `PyORBIT` can be found in [Nardiello et al. 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...664A.163N/abstract) and [Mantovan et al. 2024a](https://ui.adsabs.harvard.edu/abs/2024A%26A...682A.129M/abstract).


:::{admonition} No back-compatibility
:class: caution

Version 11 may not not be compatible with the results obtained with version 10.
Version 10 and11 are not compatible with the results obtained with previous versions .
:::



# Documentation updates

Documentation is being updated slowly but steadily, with new pages appearing every week. Below you can check which pages have been updated recently (last update: August 2024).

- *Quickstart*
  + Radial velocities
  + Lightcurve analysis
  + Combined fit
- *Data modeling*
  + Planetary RV signal
  + Correlated datasets
  + Gaussian processes regression **updated**
  + Multidimensional GPs **updated**
  + TinyGP caveats
  + CHEOPS detrending **coming soon**
  + Polynomial **coming soon**

```{tip}
The configuration files used in this documentation can be found at [this repository](https://github.com/LucaMalavolta/PyORBIT_examples)
```


## References

If you find `PyORBIT` useful for your work, please cite [Malavolta et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016A%26A...588A.118M/abstract) and [Malavolta et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018AJ....155..107M/abstract).
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
advanced_use
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
