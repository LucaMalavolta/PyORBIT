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
     K2:
       file: K2_lighcurve.dat
       kind: Phot
       models: ['gp_quasiperiodic']

- ``K2``: the label to identify the dataset. The same label must be used later in the file if we need to specify a property of the dataset (e.g. a prior on the offset).
- ``file``: the file including the dataset.
- ``kind``: the kind of dataset provided. This label is used only by specific models that requires a special treatment with respect to standard datasets, i.e. central times of transit must be tagged with the ``Tcent`` data type.
- ``models``: list of models' labels to be used to analyze the dataset.
