.. _prepare:

Prepare a YAML file
===================

Let's suppose we want to use a Gaussian Process to determine the rotational period and the timescale of decay of active regions.
PyORBIT already implements such kind of analysis, with the formalism used by XXX


Prepare the dataset
-------------------

Generic input data file must have this structure:

- 1st column: epoch of the observation
- 2nd column: independent measurement
- 3rd column: error associated to the measurement
- 4th column: flag to activate the jitter parameter(s)
- 5th column: flag to activate the offset parameter(s)
- 6th column: flag to activate the linear trend parameter(s)

In the case of central time of transits, the first column identify the number of the transit (to keep into account missing T0s)

The flags of the last three columns bust be espressed as integers. The value of a flag must be set ``-1`` if you don't want to include the corresponding parameter in the model, otherwise to increasing number starting from ``0``.
The flags can be used to divide the dataset in groups where different parameters are used for a specific model. For example, it is possible to use different offset parameters for data taken before and after a given epoch (for example, if the instrument has been modified in some way). To do so, set to ``0`` the offset flag of the data taken before the chosen epoch, and to ``1`` the data taken after that epoch. Just increase by another unit if you want to add an additional offset parameter.
Jitter, offset and linear parameters cannot be shared between datasets. For linear trends of physical origin, a model that can be shared between datasets is avaialble.

The code will assume that the number of parameters is equal to the maximum value of the flag plus one, so pay attention in increasing the flag sequentially and without jumps.


This is an example of a K2 lightcurve, with epochs expressed as BJD_{UTC}-2450000.0 Juliand days:

::

  7139.6107   -13.6377   0.0031   0   0   -1
  7139.6311   -13.6375   0.0031   0   0   -1
  7139.6515   -13.6375   0.0031   0   0   -1
  .........   ........   ......   .   .   ..

When provided to the pyorbit_emcee() subroutine instead of being read from a file, the input_dataset must be a dictionary where each keyword corresponds to the label of the dataset. For each keyword, a n*6 numpy array must be
included, where n is the number of observations.  If the keyword for a dataset is present, it will have priority on the dataset file unless the keyword is empty.

For example:

::

  input_dataset{'RV'} = np.zeros([n,6])
  input_dataset{'RV'}[:,0] = epoch of observation
  input_dataset{'RV'}[:,1] = value of observation
  input_dataset{'RV'}[:,2] = associated error
  input_dataset{'RV'}[:,3] = jitter flag
  input_dataset{'RV'}[:,4] = offset flag
  input_dataset{'RV'}[:,5] = linear trend flag

For `Tcent` kind of data, the structure is slightly different, since the first columns denote the

::

  0   4959.70736   0.00145   0   -1   -1
  1   4968.99347   0.00225   0   -1   -1
  2   4978.28014   0.00202   0   -1   -1
  .   ..........   .......   .   ..   ..

and

::

  input_dataset{'Tcent_b'} = np.zeros([n,6])
  input_dataset{'Tcent_b'}[:,0] = number of transit (e.g. if some transit is missing)
  input_dataset{'Tcent_b'}[:,1] = transit time
  input_dataset{'Tcent_b'}[:,2] = associated error
  input_dataset{'Tcent_b'}[:,3] = jitter flag
  input_dataset{'Tcent_b'}[:,4] = should be set to -1 to avoid unphysical solution
  input_dataset{'Tcent_b'}[:,5] = should be set to -1 to avoid unphysical solution


What's the point of using the flags instead of creating different datasets? Here a few examples:
- Suppose your instrument undergoes some slight modifications, and the zero point of the RV is shifted but the overall instrument is the same: you can account for this zero-point difference without using a different jitter parameter.
- Again your instrument undergoes major changes and both the zero-point and jitter are affected. However, observational parameters that depend on the characteristics of the instrument will be the same (e.g. the amplitude of stellar activity signals observed at optical wavelength), so you want to use only one parameter for this dataset and a different one for another dataset (e.g. observations gathered in the infrared).
Shortly, the flags represent a way to separate instrumental issues from the physical problems.


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
