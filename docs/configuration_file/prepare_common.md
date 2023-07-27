(prepare_common)=

# Include common objects

The physical *objects* of your system must be included in the ``common`` section. The name of the section derives from the fact each specific model required to model a given dataset will rely on physical parameters of planets and star (and more...) that are shared with all the other models eventually employed. For example, a simultaneous fit of a radial velocity curve and a transit light curve will be characterized by a single period for the planet of interest.

Let's look at this example where we have included a transiting planet in a circular orbit (*b*), a non-transiting planet (*c*), some properties for the stellar activity, and of course the stellar parameters.

```{eval-rst}
.. code-block:: yaml
   :linenos:

   common:
     planets:
       b:
         orbit: circular
         use_time_inferior_conjunction: True
         boundaries:
           P: [0.25, 0.31]
           K: [0.01, 100.0]
           Tc: [7744.06, 7744.08]
         priors:
           P: ['Gaussian', 0.2803244, 0.0000015]
           Tc: ['Gaussian', 7744.07160, 0.00022]
         spaces:
           P: Linear
         fixed:
           i: [86.3, 3.00]
       c:
         orbit: keplerian
         boundaries:
           P: [1.0, 1000.0]
           K: [0.01, 100.0]
           e: [0.00, 0.95]
     stellar_activity:
       kind: activity
       boundaries:
         Prot: [5.00, 16.00]
         Pdec: [10.0, 2000.00]
         Oamp: [0.01, 0.60]
       priors:
         Prot: ['Gaussian',13.941613, 0.66]
         Pdec: ['Gaussian',12.203273, 3.21]
         Oamp: ['Gaussian',0.342221, 0.054]
     star:
       star_parameters:
         priors:
           radius: ['Gaussian', 0.681, 0.018]
           mass: ['Gaussian', 0.708, 0.028]
           density: ['Gaussian', 2.244, 0.161]

```


Quite a lot to process, right? Let's start with the main sections. ``planets`` and ``star`` are kind of special because the section names are also the reference name of the objects, i.e., these names are hard-coded and if you try to put planet or star parameters in sections with different names you will break everything. The reason is that ``planets``  and ``star`` are actually containers for the true objects, which are ``b``, ``c`` (see the relative documentation for more details). ``stellar_activity`` instead is a label for the object with reference names ``activity``, e.g., if you want to know more you have to look for the object named ``activity`` in the documentation and for the file ``activity.py`` in the source code. Note that if you are including just a single object of a kind, you can use its reference name as a label and omit the ``type`` keyword, as in the actual example file in the repository.

.. code-block:: yaml
   :linenos:

   common:
     ...
     activity:
       boundaries:
         Prot: [5.00, 16.00]
         Pdec: [10.0, 2000.00]
         Oamp: [0.01, 0.60]
       priors:
         Prot: ['Gaussian',13.941613, 0.66]
         Pdec: ['Gaussian',12.203273, 3.21]
         Oamp: ['Gaussian',0.342221, 0.054]

Note that also ``b`` and ``c`` are labels, i.e., you can give the planets any random name, but in their case it's not necessary to specify what kind of object they are because it's implicitly assumed for being inside the ``planet`` section (**Note:** this is true only for the planets! ). The use of labels will be clear when dealing with multi-band photometry.

Let's focus on planet ``b`` as an example of the properties of an object. There are some keywords that are specific of the planet model, such as ``orbit`` and ``parametrization``. Check the documentation of each model to know more about the specific keywords. Then there are four sections that every model must possess:

- ``boundaries``: all the Bayesian samplers implemented in ``PyORBIT`` require the definition of lower and upper boundaries. If not specified, very broad boundaries will be set by default according to the source code.
- ``priors``: is it a Bayesian analysis if you don't specify at least one prior? When not included, the prior is assumed to be Uniform over the range specified in ``boundaries`` (or the one assumed by default)
- ``spaces``: the choice here is between ``Linear`` and ``Logarithmic``. Note that the logarithmic space is in *base 2*. Once again, the default choice is listed in the source code.
- ``fixed``: You listed here every parameter you want to keep fixed. The value must be accompanied by his errors, because in some cases it will be used in the computation of the derived parameters, e.g., the real mass of a transiting planet.

Under this section you need to list only the parameters that will be actually used by the models in ``models`` section. For example, a circular orbit does not require boundaries or priors for the eccentricity.

The default choices for each possible parameter and for each section listed above are declared in the source file of the object.

<!--- 

Include the models
------------------

In this section of the configuration file, called ``models``, we specify the properties of the models that we want to employ to analyze the data.

.. code-block:: yaml
   :linenos:

   models:
     rv_model:
       kind: radial_velocities
       planets:
         - b
         - c
     gp_regression:
       kind: gp_quasiperiodic
       common: stellar_activity
       RVdata:
         boundaries:
           Hamp: [0.01, 100.00]

In this example, our complete model comprises a ``radial_velocities`` model to model the orbital motion of the star due to the presence of planets ``b`` and ``c``, and a ``gp_quasiperiodic`` to model the stellar activity with Gaussian process regression and a quasi-periodic kernel. Note that ``rv_model`` and ``gp_regression`` are the labels assigned to the two models, and they are the string that has to be referenced in the ``models`` section under each ``dataset``.

There are two main sections:

- ``kind``: the model you want to employ, e.g., how the physical parameters are converted into theoretical predictions for the observations.
- ``common``: the list of labels referring to the common objects you want to be used in the model. For RVs and TTVs the keyword ``planet`` can be used as well.

In the following sections, e.g. ``RVdata`` in this example, the properties of parameters that depend specifically on the dataset are listed. The properties are ``boundaries``, ``priors``, ``spaces``, and ``fixed``, similarly as in :ref:`common-label`. Here for example we are specifying the boundaries of the amplitude of the covariance matrix in the GP regression when applied to the radial velocity data.

Additional keywords may be present depending on the model, see the documentation for more details.

**Note:**: the ``star_parameters`` object is included by default whenever needed, so you don't need to list it in the common section.

Additional parameters
---------------------

System-wide parameters that did not find place in any other section below are included in the ``parameters`` section.

.. code-block:: yaml
   :linenos:

   parameters:
     Tref: 7800.0

In this example, ``Tref`` is the epoch of reference, one of the most neglected orbital elements ever. For non-transiting planets, the argument of periapsis and the mean anomaly will be referred to this value. When not explicitly stated, it will be computed internally as the average of all the observational epochs.

Sampler parameters
------------------

Each sampler comes with its set of parameters, which fine-tuning depends on both the size of the datasets and the complexity of the model, among other things. These parameters can be specified in the configuration file under the ``solver`` section.

.. code-block:: yaml

  solver:
    pyde:
      ngen: 4000
      npop_mult: 8
    emcee:
      npop_mult: 8
      nsteps: 20000
      nburn: 5000
      thin: 100
      nsave: 10000
    nested_sampling:
      nlive: 1000
      num_repeats_mult: 5
      sampling_efficiency: 0.30
      shutdown_jitter: True
    recenter_bounds: True

This is a brief explanation of the parameters associated to each keyword, please refer to the sampler documentation for their proper usage.
- ``pyde``: parameters for the global optimization code `PyDE`_.
  - ``ngen``: number of generations.
  - ``npop_mult``: the size of the parameter vector population is given by the dimensionality of the problem multiplied by this number
- ``emcee``: parameters for the ensemble sampling toolkit for affine-invariant MCMC `emcee`_).
  - ``npop_mult``: the number of walkers in the ensemble is given by the dimensionality of the problem multiplied by this number. If PyDE and emcee are used sequentially, this keyword must have the same value in both sections (they are named in the same way a s a reminder).
  - ``nsteps``: number of steps of each chain.
  - ``nburn``: number of 'burn-in' steps.
  - ``thin``: thinning factor, should be at least equal to the autocorrelation time (before thinning). **Note:**: the chains will be saved with the thinning factor already applied
  - ``nsave``: results are saved every (unthinned) ``nsave`` steps, so that it is possible to perform a preliminary analysis while the code is still running
- ``nested_sampling``: these parameters are shared between all the implemented nested sampling algorithms, which are `MultiNest`_,  `PolyChordLite`_, and `dynesty`_.
  - ``nlive``: total number of live points.
  - ``num_repeats_mult``: The number of slice slice-sampling steps to generate a new point  (PolyChord only).
  - ``sampling_efficiency``: sampling efficiency (MultiNest only)
  - ``shutdown_jitter``: if True (default value), the jitter parameters are removed from the model (even if the flag in the dataset is active).
- ``recenter_bounds``: after the first run with (global) optimization, the boundaries of circular parameters (e.g. angles) will be recenter around the most likely value, in order to avoid border effects associated with some samplers.


.. _PyDE: https://github.com/hpparvi/PyDE
.. _emcee: https://github.com/dfm/emcee
.. _MultiNest: https://github.com/farhanferoz/MultiNest
.. _PolyChordLite: https://github.com/PolyChord/PolyChordLite
.. _dynesty: https://github.com/joshspeagle/dynesty


.. todo::
  Add links to the list of priors
  Add links to abstract_common model
  Add links to star and planet models
  Add documentation
--->
