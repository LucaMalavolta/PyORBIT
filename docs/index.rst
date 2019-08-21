.. PyORBIT documentation master file, created by
   sphinx-quickstart on Sat Dec 16 17:20:15 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyORBIT: the ultimate tool for exoplanet characterization!
===================================

``PyORBIT`` is a robust, versatile framework for the characterization of planetary systems.
With ``PyORBIT`` you can model light curves, radial velocities, activity indexes and transit time variations.
In addition to the exoplanet signature, you can add to your model instrumental systematics and stellar activity,
including Gaussian processes regression, among other things. Different parametrization can be used for orbital
parameters, and orbits can be computed using dynamical integration or just non-interacting Keplerians.
For every parameter in the model, it is possible to set a prior, explore it in linear or logarithmic space,
or keep it fixed if necessary. Thanks to abstraction, it is virtually possible to implement any physical model you can think of.

Setting-up a configuration file and running the analysis is super easy, barely an inconvenience, even if you have never wrote a line of Python (or any language at all) in your life, thanks to the versatility of YAML.
Alternatively, for easy automatization, PyORBIT can be called as a Python function by passing a dictionary instead of a configuration file.

..
  One of the main strength of PyORBIT is in the abstraction of the models. For 
  example,  you can use more than one Gaussian process, with independent
  hyper-parameters, just by using different labels for the abstract class.
  Alternatively, for your Gaussian process regression you can use "celerite" on
  your light curve and "george " on your radial velocities, by having them sharing
  the same hyperparameters.

``PyORBIT`` is written and tested for Python 3, but is back-compatible with Python 2 (at your own risk)

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   prepare_datasets
   prepare_yaml
   api



Following `PEP 8 Style Guide for Python Code <https://www.python.org/dev/peps/pep-0008/>`_  ,
`PEP 257 Docstring Conventions <https://www.python.org/dev/peps/pep-0257/>`_ and `Google Python Style Guide <http://google.github.io/styleguide/pyguide.html>`_


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
