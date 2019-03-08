.. _prepare_yaml:

Prepare a YAML file
===================

Let's suppose we want to use a Gaussian Process to determine the rotational period and the timescale of decay of active regions.
PyORBIT already implements such kind of analysis, with the formalism used by XXX


Adding a dataset
----------------

Datasets are grouped under the input category:

.. code-block:: yaml
   :linenos:

   inputs:
     RVdata:
       file: RVdata.dat
       kind: RV
       models:
        - gp_quasiperiodic
        - radial_velocities

- ``RVdata``: the label to identify the dataset. The same label must be used later in the file if we need to specify a property of the dataset (e.g. a prior on the offset).
- ``file``: the file including the dataset.
- ``kind``: the kind of dataset provided. This label is used only by specific models that requires a special treatment with respect to standard datasets, i.e. central times of transit must be tagged with the ``Tcent`` data type.
- ``models``: list of models' labels to be used to analyze the dataset. If only one model is specificed, it can be written on the same line without the ``-`` sign.

Dataset kinds
+++++++++++++
The following datasets are recognized by the code. In the majority of the cases, the way the dataset is treated depends on the specified models. A list of aliases is defined to circumvent the most common typos (and to avoid checking the documentation every time).

- ``RV``: radial velcoities. Aliases: ``RV``, ``RVs``, ``rv``, ``rvs``.
- ``Phot``: photometry. Aliases: ``P``, ``Ph``, ``p``, ``ph``, ``PHOT``, ``Phot``, ``phot``, ``Photometry``, ``photometry``.
- ``Tcent``: central times of transits. Aliases: ``Tcent``, ``TCent``, ``Tc``, ``TC``, ``T0``, ``TT``.
- ``FWHM``: Full Width Half Maximum of the Cross-Correlation Function (CCF). Aliases: ``FWHM``, ``fwhm``.
- ``BIS``: Bisector Inverse Span of the CCF. Aliases: ``BIS``, ``bis``.
- ``H-alpha``: measurements of the emission in the core of the H-alpha line. Aliases: ``H``, ``HA``, ``h``, ``ha``, ``Halpha``, ``H-alpha``, ``halpha``, ``h-alpha``.
- ``S_index``: (uncalibrated) measurements of the emission in the core of the CA H&K  lines. Aliases: ``S``, ``S_index``.
- ``Ca_HK``: as for the S index, but with the photometric contribution removed. Aliases: ``Ca_HK``, ``logR``.
