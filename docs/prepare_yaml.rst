.. _prepare_yaml:

Prepare a configuration file
============================

The configuration file is the most important part of your analysis, because it
will contain the specifics of your model, the priors you want to assign to the
parameters, and the way you want to analyze your data. Thanks to the readability
of the YAML language, it will also serve as a useful memento at a glance for the
set-up of your model.

File structure
--------------

The file is divided in several sections, each one with a clear purpose
- ``inputs``: here you list the files containing the datasets
- ``common``: the physical parameters of your planetary system that are
going to be shared between the models (and are independent from them): orbital
parameters of the planets, mass and radius of the star, rotational period and
decay timescale of the active regions...
- ``models``: this is where you specific how the physical elements of the system
in the *common* section are going to influence your data: which planets include
in the radial velocity computation, what kind of model/kernel you want to use to
model the activity...
- ``parameters``: values that are not system parameters per se, but are required
for the correctness of the analysis
- ``solver``: the parameters required by PyDE, emcee, or MultiNest are all listed
here

To have a glance at how a configuration file looks like, check the .. _documentation_example.yaml: http://cnn.com/:

Adding a dataset
----------------

Datasets are grouped under the input section:

.. code-block:: yaml
   :linenos:

   inputs:
     RVdata:
       file: K2-141_dataset/K2-141_RV_HARPN.dat
       kind: RV
       models:
         - radial_velocities
         - gp_quasiperiodic

- ``RVdata``: the label to identify the dataset. The same label must be used
 later in the file if we need to specify a property of the dataset.
- ``file``: the file including the dataset.
- ``kind``: the kind of dataset provided. This label is used only by specific
models that requires a special treatment with respect to standard datasets,
i.e. central times of transits must be tagged with the ``Tcent`` data type.
- ``models``: list of models' labels to be used to analyze the dataset. If
only one model is specified, it can be written on the same line without
the ``-`` sign.

Dataset kinds
+++++++++++++
The following datasets are recognized by the code. In the majority of the cases, the way the dataset is treated depends on the specified models. A list of aliases is defined to circumvent the most common typos (and to avoid checking the documentation every time...).

- ``RV``: radial velocities. Aliases: ``RV``, ``RVs``, ``rv``, ``rvs``.
- ``Phot``: photometry. Aliases: ``P``, ``Ph``, ``p``, ``ph``, ``PHOT``, ``Phot``, ``phot``, ``Photometry``, ``photometry``.
- ``Tcent``: central times of transits. Aliases: ``Tcent``, ``TCent``, ``Tc``, ``TC``, ``T0``, ``TT``.
- ``FWHM``: Full Width Half Maximum of the Cross-Correlation Function (CCF). Aliases: ``FWHM``, ``fwhm``.
- ``BIS``: Bisector Inverse Span of the CCF. Aliases: ``BIS``, ``bis``.
- ``H-alpha``: measurements of the emission in the core of the H-alpha line. Aliases: ``H``, ``HA``, ``h``, ``ha``, ``Halpha``, ``H-alpha``, ``halpha``, ``h-alpha``.
- ``S_index``: (uncalibrated) measurements of the emission in the core of the CA H&K  lines. Aliases: ``S``, ``S_index``.
- ``Ca_HK``: as for the S index, but with the photometric contribution removed. Aliases: ``Ca_HK``, ``logR``.

Include common models
---------------------

The physical *objects* of your system must be included here. Let's have a look
at this model where we have included a transiting planet in a circular orbit (*b*),
a non-transiting planet (*c*), some properties for the stellar activity, and of
of course the stellar parameters.

.. code-block:: yaml
   :linenos:

   common:
     planets:
       b:
         orbit: circular
         parametrization: Eastman2013_Tcent
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
       type: activity
       boundaries:
         Prot: [5.00, 16.00]
         Pdec: [10.0, 2000.00]
         Oamp: [0.01, 0.60]
       priors:
         Prot: ['Gaussian',13.941613, 0.66]
         Pdec: ['Gaussian',12.203273, 3.21]
         Oamp: ['Gaussian',0.342221, 0.054]
     star:
       label_name:
         kind: star_parameters
         priors:
           radius: ['Gaussian', 0.681, 0.018]
           mass: ['Gaussian', 0.708, 0.028]
           rho: ['Gaussian', 2.244, 0.161]

Quite a lot to process, right? Let's start with the main sections. ``planets``
and ``star`` are kind special because the section names are also the reference
name of the objects, i.e., these names are hard coded and if you try to put
planet or star parameters in sections with different names you will break
everything. The reason is that ``planets``  and ``star`` are actually containers
for the true objects, which are ``b``, ``c`` (see the
relative documentation for more details)

``stellar_activity`` and ``label_name`` instead are labels for the object with
reference names ``activity`` and ``star_parameters`` respectively, e.g., if you
 want to know more you have to look for the object named ``activity`` in the
 documentation  and for the file ``activity.py`` in the source. Note that if
you are including just a single object of a kind, you can use its reference name
as label and omit the ``type`` keyword, as in actual example file in the repository.

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
     star:
       star_parameters:
         priors:
           radius: ['Gaussian', 0.681, 0.018]
           mass: ['Gaussian', 0.708, 0.028]
           rho: ['Gaussian', 2.244, 0.161]

Note that also ``b`` and ``c`` are labels, i.e., you can give the planets any
random name, but in their case it's not necessary to specify what kind of object
they are because it's implicitly assumed for being inside the ``planet`` section
 (**Note** this is true only for the planets! )

Let's focus on planet ``b`` as an example of the properties of an object. There
are some keywords that are specific of the planet model, such as ``orbit`` and
``parametrization``. Check the
documentation of each model to know more about the specific keywords. Then there
 are four sections that every model must possess:
- ``boundaries``: all the Bayesian samplers implemented in ``PyORBIT`` require
the definition of lower and upper boundaries. If not specified,
very broad boundaries will be set by default according to the source code.
- ``priors``: is it a Bayesian analysis if you don't specify at least one prior?
When not included, the prior is assumed to be Uniform over the range specified
 in ``boundaries`` (or the one assumed by default)
- ``spaces``: the choice here is between ``Linear`` and ``Logarithmic``.
Note that the logarithmic space is in *base 2*. Once again, the default choice is
listed in the source code.
- ``fixed``: You listed here every parameter you want to keep fixed. The value
 must be accompanied by his errors, because in some cases it will be used in the
  computation of the derived parameters, e.g., the real mass of a transiting planet.

Under this section you need to list only the parameters that will be actually used by the
models in ``models`` section. For example, a circular orbit does not require boundaries or priors for the eccentricity.

The default choices for each possible parameter and for each section listed above
are declared in the relative file of the object

.. todo::
  Add links to the list of priors
  Add links to abstrac_common model
  Add links to star and planet models
  Add documentation
