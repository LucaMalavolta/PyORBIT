PyORBIT: the ultimate tool for exoplanet characterization!
==========================================================

``PyORBIT`` is a robust, versatile framework for the characterization of planetary systems.
With ``PyORBIT`` you can model light curves, radial velocities, activity indexes and transit time variations.
In addition to the exoplanet signature, you can add to your model instrumental systematics and stellar activity,
including Gaussian processes regression with a variety of kernels. Both Markov Chain Monte-Carlo samplers or Nested Sampling algortihms are supported.
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

``PyORBIT`` started in 2015 as an exercise to learn Python (2), and also because at
the time the only publicly available code for the Bayesian analysis of RVs was
written in IDL, for which I didn't have a license. Since then I've been added new options
every time I needed so, and I kept updating the code while improving my Python skills.

``PyORBIT`` has been now converted and tested for Python 3, but it should still
be back-compatible with Python 2 (at your own risk).

If you are wondering what ``PyORBIT`` stands for: I am bad at acronym creation so
I decided to write it with capitol *ORBIT* just because I liked how it looked.
Feel free to submit your retrofitting acronym!  

[![GitHub][github-badge]][github-link]

[github-badge]: https://img.shields.io/badge/GitHub-PyORBIT-blue
[github-link]: https://github.com/LucaMalavolta/PyORBIT

..
 ```{eval-rst}
 .. image::  https://img.shields.io/badge/GitHub-PyORBIT-blue
  :target: https://github.com/LucaMalavolta/PyORBIT
 ```


```{admonition} Be careful
Documentation is still work in progress
```

```{eval-rst}
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   prepare_datasets
   prepare_yaml
   api



..
  Following `PEP 8 Style Guide for Python Code <https://www.python.org/dev/peps/pep-0008/>`_  ,
  `PEP 257 Docstring Conventions <https://www.python.org/dev/peps/pep-0257/>`_ and `Google Python Style Guide <http://google.github.io/styleguide/pyguide.html>`_

..
  Indices and tables
  ==================

  * :ref:`genindex`
  * :ref:`modindex`
  * :ref:`search`


..
 https://sphinx-book-theme.readthedocs.io/
 https://myst-parser.readthedocs.io/


```
