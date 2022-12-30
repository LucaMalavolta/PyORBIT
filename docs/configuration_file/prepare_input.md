(prepare_input)=

# Adding a dataset

Datasets are grouped under the input section:

```{eval-rst}
.. code-block:: yaml
   :linenos:

   inputs:
     RVdata:
       file: K2-141_dataset/K2-141_RV_HARPN.dat
       kind: RV
       models:
         - rv_model
         - gp_regression
```

- ``RVdata``: the label to identify the dataset. The same label must be used later in the file if we need to specify a property of the dataset.
- ``file``: the file including the dataset.
- ``kind``: the kind of dataset provided. This label is used only by specific models that requires a special treatment with respect to standard datasets, i.e. central times of transits must be tagged with the ``Tcent`` data type.
- ``models``: list of models' labels to be used to analyze the dataset. If only one model is specified, it can be written on the same line without the ``-`` sign.

## Dataset kinds

The following datasets are recognized by the code. In the majority of the cases, the way the dataset is treated depends on the specified models. A list of aliases is defined to circumvent the most common typos (and to avoid checking the documentation every time...).

- ``RV``: radial velocities. Aliases: ``RV``, ``RVs``, ``rv``, ``rvs``.
- ``Phot``: photometry. Aliases: ``P``, ``Ph``, ``p``, ``ph``, ``PHOT``, ``Phot``, ``phot``, ``Photometry``, ``photometry``.
- ``Tcent``: central times of transits. Aliases: ``Tcent``, ``TCent``, ``Tc``, ``TC``, ``T0``, ``TT``.
- ``FWHM``: Full Width Half Maximum of the Cross-Correlation Function (CCF). Aliases: ``FWHM``, ``fwhm``.
- ``BIS``: Bisector Inverse Span of the CCF. Aliases: ``BIS``, ``bis``.
- ``H-alpha``: measurements of the emission in the core of the H-alpha line. Aliases: ``H``, ``HA``, ``h``, ``ha``, ``Halpha``, ``H-alpha``, ``halpha``, ``h-alpha``.
- ``S_index``: (uncalibrated) measurements of the emission in the core of the CA H&K  lines. Aliases: ``S``, ``S_index``.
- ``Ca_HK``: as for the S index, but with the photometric contribution removed. Aliases: ``Ca_HK``, ``logR``.

.. _common-label:
