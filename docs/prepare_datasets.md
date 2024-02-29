(prepare_dataset)=

# Prepare a dataset file

## How to properly format a dataset

`PyORBIT` accepts its specific dataset format, and it can be very picky about it!

In addition to the usual stuff (time, measurement, measurement error) a few additional columns must be specified in order to include instrument-specific parameters to your model. These parameters are:

- **jitter**: a value added in quadrature to the error estimates, useful if the error estimates are underestimated or if their estimation did not include addisional sources of (white) noise
- **offset**: the baseline value of the dataset on top of which our model is added, for example, the systematic radial velocity of the star or the average value of the Full Width Half Maximum of the spectral lines
- **any other signal of instrumental origin** for example a trend in the RV or the FWHM due to instrumental problems during a specific time interval

A generic input data file must have this structure:

- 1<sup>st</sup> column: time of the observation
- 2<sup>nd</sup> column: independent measurement (e.g. RV)
- 3<sup>rd</sup> column: error associated with the measurement
- 4<sup>th</sup> column: flag to activate the jitter parameter(s)
- 5<sup>th</sup> column: flag to activate the offset parameter(s)
- 6<sup>th</sup> column: flag to activate the subset modeling


```{tip}
  Usually it's always a good idea to include a *jitter* term, while the *offset* column may be required or not depending on the kind of dataset, while the last column applies only in specific cases.
```

Flags must be expressed as integer numbers following the *Pythonic* way of counting numbers: **0** to activate a flag, **-1** to deactivate it.
A generic radial velocity dataset should look like this:

```
  2456000.010493  4764.73  1.20  0  0  -1
  2456000.719975  4766.58  1.35  0  0  -1
  2456001.967132  4779.52  1.23  0  0  -1
  ..............  ........ ....  .  .  ..
```

Flags can be used to divide the dataset into groups with different jitter parameters or offset parameters, not necessarily correlated.
A common case is the change in the RV offset of HARPS observations after the intervention in 201X. A new offset parameter can be assigned to the observations after the intervention simply by increasing (+1) the value of the flag.

```
  2456000.010493  4764.73  1.20  0  0  -1
  2456000.719975  4766.58  1.35  0  0  -1
  2456001.967132  4779.52  1.23  0  1  -1
  2456002.447132  4779.52  1.23  0  1  -1
  2456002.337132  4779.52  1.23  0  1  -1
  ..............  ........ ....  .  .  ..
```

In the example above, we have decided to use the same jitter term regardless of the intervention.

Generally speaking, `PyORBIT` will assume that the number of parameters is equal to the maximum value of the flag plus one, so pay attention to increasing the flag sequentially and without jumps. Follow these guidelines for a simple and happy life:

- Flags must be given in consecutive order starting from zero (Python notation).
- The inactive flags must be set to -1.
- All the observations that share the same flag value will have the corresponding parameter in common *within the dataset*.
- Different parameters will be used for measurements with different values of the corresponding flag.
- Flags in different *columns* are independent.
- Flags in different *files* are independent

```{warning}
  If a column is missing, `PyORBIT` will assume that the corresponding flag is deactivated. However, columns are not labelled, so it is not possible to deactivate the *jitter* column by removing it without deactivating the *offset* column as well.
```

<!---
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

 --->

## Dealing with several datasets of the same type

If you are dealing with observations of the same type from different sources (for example, RVs taken with different instruments) you can follow two roads:

1) Put everything together in a single file, taking care of setting the jitter and offset flag properly
2) Write a file for each dataset

Many codes expect you to follow the first road. `PyORBIT` can work with both, although in some cases you *have* to use different files when different models must be employed (for example, photometric transit observed with different instruments thus requiring different limb darkening parameters). In general, my advice is to use a file for each dataset because it will make the configuration file more self-explicative and in the long term it will make your life much easier - especially when you are looking back at the analysis after some time!

## Exceptions to standard formatting

### Central transit times

For central time of transit (`Tcent`) file data, required by TTV analysis, the structure is slightly
different. The first column identifies the number of the transit (to keep into account missing T0s). This number will help in identifying missing transits (but honestly I don't remember right now what happens if you start from a random number...)

```
  0   2454959.70736   0.00145   0   -1   -1
  1   2454968.99347   0.00225   0   -1   -1
  2   2454978.28014   0.00202   0   -1   -1
  .   .............   .......   .   ..   ..
```

```{warning}
Always set the 5th (offset) and 6th column to ``-1`` (or don't include
them at all) to avoid unphysical solution (drift and jumps in time are not allowed)
```

### Astrometry

Working on it!

### Ancillary data

Some models require one or more additional datasets to work. These datasets are used as independent variables, as such they don't need to be compared to a model and they do not enter into the calculation of the likelihood. For example, when correlating variations of the flux with the position of the star on the CCD, the latter is the independent variable. These datasets do not require jitter or offset terms, so the structure is more relaxed, but the inclusion of a header with the appropriate dataset names - detailed in the documentation of each model - is a **fundamental requirement**.


```
# time flux flux_err xoff yoff bg contam smear deltaT roll_angle
9052.138151     1.000289     0.000259     0.230865     -1.670593     0.015518     0.023160     0.000013     0.669403     194.377112
9052.138846     1.000069     0.000258     0.447083     -1.553406     0.015485     0.023059     0.000012     0.648865     192.682123
9052.139541     1.000413     0.000259     0.459320     -1.494080     0.015407     0.023058     0.000012     0.628357     191.009003
9052.140235     1.000144     0.000258     0.791229     -0.604126     0.015293     0.023094     0.000012     0.566956     189.346915
9052.140930     1.000356     0.000259     0.864838     -0.110840     0.015305     0.023187     0.000012     0.566956     187.684830
9052.141625     1.000824     0.000259    -0.025452     -0.497986     0.015294     0.023091     0.000012     0.587402     186.012203
```


<!---

::

  input_dataset{'Tcent_b'} = np.zeros([n,6])
  input_dataset{'Tcent_b'}[:,0] = number of transit (e.g. if some transit is missing)
  input_dataset{'Tcent_b'}[:,1] = transit time
  input_dataset{'Tcent_b'}[:,2] = associated error
  input_dataset{'Tcent_b'}[:,3] = jitter flag
  input_dataset{'Tcent_b'}[:,4] = should be set to -1 to avoid unphysical solution
  input_dataset{'Tcent_b'}[:,5] = should be set to -1 to avoid unphysical solution

-->

## Standard units for input data and model parameters

- **Time**: **day**. The code assumes that all epochs, time of observations, and timescales, are expressed in BJD-TDB, although no check is enforced. You can remove a constant from the epochs/time of observations without consequences (e.g. use BJD_TDB - 2450000, Kepler BJD, Tess BJD...), just be sure to do it consistently on all the datasets and the parameters in the configuration file.
- **Radial Velocities** (RV): **meter/second**. If you use kilometers/second, the fit may still work but all the derived values will be meaningless.
- **Inverse Bisector Span** (BIS): **meter/second**. As for the RVs.
- **Full Width Half Maximum** (FWHM) of the Cross-Correlation Function: **kilometer/second**. This is the standard measurement unit for this quantity.
- **Stellar mass, radius, density**: **Solar units**. Note: some catalogs report the stellar density as g/cm<sup>3</sup> or kg/m<sup>3</sup>, so be sure to convert it accordingly. 
- **Planetary radius, semimajor axis of the orbit**: **stellar radii**. These parameters are conventionally called *scaled planetary radius* and *scaled semimajor axis*, and within `PyORBIT` they are denoted with a `_Rs` subscript.
- **Angles of any kind**: **degrees** (starting from `PyORBIT` version 9)
- Physical planetary masses and radii are respectively denoted by `_Me` and `_Re` subscripts when in Earth units, and by  `_Mj` and `_Rj` subscripts when in Jupiter units.

<!---

# Jitter, offset and linear flags

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


-->
