.. _prepare_datasets:

Prepare the dataset
===================

Sometimes your instrument may be affected by some weird systematics (e.g. offsets
between measurement taken at different epochs, trends with time), these can be
easily encoded in the dataset file without burden the configuration file with extra parameters,
or even worse by splitting the data into several files.


Generic dataset
+++++++++++++++

Generic input data file must have this structure:

- 1st column: epoch of the observation
- 2nd column: independent measurement
- 3rd column: error associated to the measurement
- 4th column: flag to activate the jitter parameter(s)
- 5th column: flag to activate the offset parameter(s)
- 6th column: flag to activate the linear trend parameter(s)

In the case of central time of transits, the first column identify the number
of the transit (to keep into account missing T0s), skip to `Central transit times`_
for more information .

The flags of the last three columns must be expressed as integers. Jitter,
offset and linear parameters cannot be shared between datasets. For linear
trends of physical origin, a model that can be shared between datasets is
available. Check the `Jitter, offset and linear flags`_ subsection for more
information.

This is an extract from the file ``TestCase01_RV.dat``, in the example folder,
with epochs expressed as BJD_{TDB}-2450000.0 Julian days:

::

  6000.010493  4764.73  1.00  0  0  -1
  6000.719975  4766.58  1.00  0  0  -1
  6001.967132  4779.52  1.00  0  0  -1
  .........   ........  ....  .  .  ..

The choice of the time standard is arbitrary, just pay attention to always be
self-consistent.

When provided to the ``pyorbit_emcee()`` subroutine instead of being read from
a file, the input_dataset must be a dictionary
where each keyword corresponds to the label of the dataset. For each keyword, a
[n,6] ``numpy`` array must be
provided, where *n* is the number of observations.  If the keyword for a dataset is present, it will have priority on the dataset file unless the keyword is empty.

For example:

::

  input_dataset{'RV'} = np.zeros([n,6])
  input_dataset{'RV'}[:,0] = epochs of the observations
  input_dataset{'RV'}[:,1] = values of the observations
  input_dataset{'RV'}[:,2] = associated errors
  input_dataset{'RV'}[:,3] = jitter flag
  input_dataset{'RV'}[:,4] = offset flag
  input_dataset{'RV'}[:,5] = linear trend flag


Central transit times
+++++++++++++++++++++

For central time of transit (`Tcent`) file data, the structure is slightly
different, since the first column denote the number of the transit.
This number will help in identifying missing transits (but honestly I don't
remember right now what happens if you start from a random number...)

::

  0   4959.70736   0.00145   0   -1   -1
  1   4968.99347   0.00225   0   -1   -1
  2   4978.28014   0.00202   0   -1   -1
  .   ..........   .......   .   ..   ..

**IMPORTANT** always set the 5th (offset) and 6th column to ``-1`` (or don't include
them at all) to avoid unphysical solution (drift and jumps in time are not allowed)

::

  input_dataset{'Tcent_b'} = np.zeros([n,6])
  input_dataset{'Tcent_b'}[:,0] = number of transit (e.g. if some transit is missing)
  input_dataset{'Tcent_b'}[:,1] = transit time
  input_dataset{'Tcent_b'}[:,2] = associated error
  input_dataset{'Tcent_b'}[:,3] = jitter flag
  input_dataset{'Tcent_b'}[:,4] = should be set to -1 to avoid unphysical solution
  input_dataset{'Tcent_b'}[:,5] = should be set to -1 to avoid unphysical solution

Jitter, offset and linear flags
+++++++++++++++++++++++++++++++

Jitter, offset and linear parameters cannot be shared between datasets. For linear trends of physical origin, a model that can be shared between datasets is avaialble.

Activating the jitter flag will introduce a new parameter which will be added in quadrature to the error bars of the measurement for which the flag has been activated.
The offset flag will add a constant offset (or zero-point) to all the measurement for which the flag has been activated.
The linear flag will include a linear trend in the dataset. Note that this flag will add only the slope as additional parameter, while the intercept of the linear trend must be set using the offset flat. Only a linear trend is supported, higher order terms have not been implemented simply because I never encpuntered such an extreme case, but on request it may be included as additional columns in the file.

The flags of the last three columns must be expressed as integers. The value of a flag must be set ``-1`` if you don't want to include the corresponding parameter in the model, otherwise to increasing number starting from ``0``.
The flags can be used to divide the dataset in groups where different parameters are used for a specific model. For example, it is possible to use different offset parameters for data taken before and after a given epoch (for example, if the instrument has been modified in some way). To do so, set to ``0`` the offset flag of the data taken before the chosen epoch, and to ``1`` the data taken after that epoch. Just increase by another unit if you want to add an additional offset parameter.

The code will assume that the number of parameters is equal to the maximum value of the flag plus one, so pay attention in increasing the flag sequentially and without jumps.

For a given kind of flag:

- Flags must be given in consecutive order starting from zero (Python notation).
- Inactive flag must be set to -1.
- All the parameters that share the same flag value will have that parameter in common.
- Different parameters will be used for measurements with different value of flag.
- Flags in different columns are independent.

Let's look at the following example:

::

  epoch_00  meas_00  err_00  0  0  -1
  epoch_01  meas_01  err_01  0  0  -1
  epoch_02  meas_02  err_02  1  0  -1
  epoch_03  meas_03  err_03  1  1  -1
  epoch_04  meas_04  err_04  1  1   0
  epoch_05  meas_05  err_05  2  1   0
  epoch_06  meas_06  err_06  2  0   0


- `epoch_00` and `epoch_01` share the same jitter term, so they do `(epoch_02, epoch_03, epoch_04)` and `(epoch_05, epoch_06)`, for a total of 3 jitter parameters.
- `epoch_00`, `epoch_01`, `epoch_02` and `epoch_06` share the same offset. `epoch_03`, `epoch_04`, `epoch_05` share a different offset parameter.
- `epoch_04`, `epoch_05`, `epoch_06` are modeled using a linear trend. `epoch_00`, `epoch_01`, `epoch_02` and `epoch_03` are not.


What's the point of using the flags instead of creating different datasets? Here a few examples:

- Suppose your instrument undergoes some slight modifications, and the zero point of the RV is shifted but the overall instrument is the same: you can account for this zero-point difference while sharing the same jitter parameter.
- Again your instrument undergoes major changes and both the zero-point and jitter are affected. However, observational parameters that depend on the characteristics of the instrument will be the same (e.g. the amplitude of stellar activity signals observed at optical wavelength), so you want to use only one parameter for this dataset and a different one for another dataset (e.g. observations gathered in the infrared).

Shortly, the flags represent a way to separate instrumental issues from the physical problems.

..
 References
 ----------

 Later
