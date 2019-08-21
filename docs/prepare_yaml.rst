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

To have a glance at how a configuration file looks like, check this example: 

Adding a dataset
----------------

Datasets are grouped under the input category:

.. code-block:: yaml
   :linenos:

   inputs:
     RVdata:
       file: K2-141_dataset/K2-141_RV_HARPN.dat
       kind: RV
       models:
         - radial_velocities
         - gp_quasiperiodic

- ``RVdata``: the label to identify the dataset. The same label must be used later in the file if we need to specify a property of the dataset.
- ``file``: the file including the dataset.
- ``kind``: the kind of dataset provided. This label is used only by specific models that requires a special treatment with respect to standard datasets, i.e. central times of transits must be tagged with the ``Tcent`` data type.
- ``models``: list of models' labels to be used to analyze the dataset. If only one model is specified, it can be written on the same line without the ``-`` sign.

Dataset kinds
+++++++++++++
The following datasets are recognized by the code. In the majority of the cases, the way the dataset is treated depends on the specified models. A list of aliases is defined to circumvent the most common typos (and to avoid checking the documentation every time...).

- ``RV``: radial velcoities. Aliases: ``RV``, ``RVs``, ``rv``, ``rvs``.
- ``Phot``: photometry. Aliases: ``P``, ``Ph``, ``p``, ``ph``, ``PHOT``, ``Phot``, ``phot``, ``Photometry``, ``photometry``.
- ``Tcent``: central times of transits. Aliases: ``Tcent``, ``TCent``, ``Tc``, ``TC``, ``T0``, ``TT``.
- ``FWHM``: Full Width Half Maximum of the Cross-Correlation Function (CCF). Aliases: ``FWHM``, ``fwhm``.
- ``BIS``: Bisector Inverse Span of the CCF. Aliases: ``BIS``, ``bis``.
- ``H-alpha``: measurements of the emission in the core of the H-alpha line. Aliases: ``H``, ``HA``, ``h``, ``ha``, ``Halpha``, ``H-alpha``, ``halpha``, ``h-alpha``.
- ``S_index``: (uncalibrated) measurements of the emission in the core of the CA H&K  lines. Aliases: ``S``, ``S_index``.
- ``Ca_HK``: as for the S index, but with the photometric contribution removed. Aliases: ``Ca_HK``, ``logR``.
