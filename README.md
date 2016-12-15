# PyORBIT

**`PyORBIT` v3.1 by Luca Malavolta - 2016**    

PyORBIT handles several kinds of datasets, such as radial velocity (RV), activity indexes, and photometry, to simultaneously characterize the orbital parameters of exoplanets and the noise induced by the activity of the host star. RV computation is performed using either non-interacting Kepler orbits or n-body integration. Stellar activity can be modeled either with sinusoids at the rotational period and its harmonics or gaussian process. Offsets and systematics in measurements from several instruments can be modeled as well. Thanks to the modular approach, new methods for stellar activity modeling or parameter estimation can be easily incorporated into the code.

Most of the information can be found in the paper by [Malavolta et al. (2016)][Malavolta2016]

---

## PyORBIT user guide
1. [Introduction](### Introduction)
2. [Install and Compile](### Install and Compile)
3. [Configuration file](### Configuration file)
4. [Data files](### Data files)
5. [Run the code](### Run the code)
6. [Analyze the output](### Analyze the output)

### Introduction

I started writing this code because I needed a versatile RV fitting code, with the ability of including different models for stellar activity. At the same time I used this project to improve my knowledge of Python, which was very limited at the time I started this project. If you dig into the code you may notice some inconsistencies in the way I treat data and variables, I'm trying to homogenize the code as long I'm learning new Pythonic tricks but it will take some time.  

### Install and Compile

 `PyORBIT` is avaialble at [https://github.com/LucaMalavolta/PyORBIT](https://github.com/LucaMalavolta/PyORBIT)

  Download the .zip files or clone the repository.

  To take advantage of dynamical integration, you need to install either `TRADES` ([Borsato et al. 2014][Borsato2014], available at [here](https://github.com/lucaborsato/trades)) or `TTVFast` ([Deck et al. 2014][Deck2014]) using the Python wrapper interface available [here](https://github.com/mindriot101/ttvfast-python).
  After installing `TRADES`, you have to specify its path in `PyORBIT_V3_Classes.pyx` (before importing the `pytrades_lib` library).
  `PyORBIT` can work without these codes installed, however some editing is required. Furthermore, the file `pytrades/constants.py` from the `TRADES` repository must be copied in the `PyORBIT` main folder

  To compile the program using `Cython`, simply execute `./compile.bash` in a `Bash` terminal. `kepler_exo.py` and `PyORBIT_V3_Classes.py` provide symbolic links to the respective `*.pyx` files, so that the program can be executed without using Cython. Be aware that the `*.so` libraries will be preferred to the `*.py` symbolic links, so remember to recompile or delete the `*.so` files if you are modifying the files.

  Older versions of the code can be found in the `V1` and `V2` folders.

### Configuration file

This file contains all the information regarding:
* datasets to be analyzed (e.g. RV, photometry, stellar activity indicators)
* number of planets and their characteristics (e.g. priors on orbital parameters, )
* models to be used to fit the data (e.g. gaussian processes for activity, long-term trends)


```yaml
Name: PyORBIT_test1
Output: PyORBIT_test1
```
The first two keywords identify the name of the planet and the name of the directory where all the output files are stored.

```yaml
Inputs:
  0:
      File: test1_PH.dat
      Kind: Phot
      Models: ['sinusoids']
  1:
      File: test1_RV.dat
      Kind: RV
      Models: ['kepler', 'sinusoids']
```
The ```Input``` section includes all the data files to be analyzed. The only exceptions are the files containing the central time of transits, which must be included in the ```Planets``` section and not here.  Following the Python standard, file must be listed starting from ```0``` without gaps.

The main goal of ```Kind``` is to distinguish between Radial Velocities (```RV```) and other kinds of data ( ```Phot```, while ```BIS```, ```FWHM``` and ```Act```, although only in the case of the ```sinusoid``` model a real distinction is made).

```Models``` specificies how you want to fit the data. The avaialble models are:
  * ```kepler```: keplerian curve, avaialble only for RVs.
  * ```curvature```: polynomial trend of order _n_ . The same trend is applied to each dataset with this keyword listed among the models. Instrumental trend (specific of a given dataset) can be corrected with a keyword inside the dataset file
  * ```sinusoid```: activity modelling with sinusoids at the rotational period of the star and its harmonics, following [Boisse et al. (2011)][Boisse2011]
  * ```gaussian```: activity modeling with gaussian processes

```yaml
Planets:
  0:
    Orbit: dynamical
    Boundaries:
      P: [20.0, 50.0]
      M: [10.0, 1000.0]
      e: [0.00, 0.40]
    Fixed:
      lN : 0.0
      i: 90.0
    Priors:
      P: ['Gaussian', 35.50, 0.01]
    Starts:
      o: 0.345
      f: 2.145
    Tcent: obsT0.dat
    Radius: [11.200, 0.100]
    Inclination: [90.0, 0.01]
```
In this section the characteristics of the planets can be specified. As before, the planets must be numbered starting from zero and without gaps.

*  ```Orbit```: the model to be used for the RV planetary signal:
  * ```circular```: circular orbit with the period of the planet ```P```, the RV semi-amplitude ```K``` and the phase ```f``` (not to be confused with the _true anomaly_) as free parameters.
  * ```keplerian```: keplerian (eccentric) orbit, same parameters as *circular* with the addition of the eccentricity ```e``` and the argument of pericenter ```o```
  * ```dynamical```: RVs are modelled using a N-body integration instead of non-interacting keplerians. The parameters in the fit are ```P```, ```e```, ```f```, ```o``` as in the previous models, the mass of the planet ```M```, the inclination of the planet ```i```, the longitude of the ascending node ```lN```.

  Note that ```f = o + mA``` where _mA_ is the mean anomaly at the epoch of reference  ```Tref``` (to be introduced later)
* ```Boundaries```: if not otherwise specified, an uninformative prior within the specified boundaries is chosen. Since _P_ and _K_ are explored in the logarithmic space, their lower limits must be always greater than zero.  _P_ is in days, _K_ in meter/second,  all the angles are in radians except _i_ which is in degrees.

* ```Fixed```: use this keyword to fix a parameter to a given value. This is useful if a parameter cannot be contrained by the data.

* ```Priors```: use this keyword to specify a prior on a variable. The selected function for the  probability distribution must be followed by the appropriate numbers of parameters, which depend in the implementation of such function. Right known only ```Gaussian``` and ```Uniform``` function have been implemented.

* ```Starts```: starting points for the MCMC chains (when required). If not defined, the starting points of the emcee chains are determined by running the global optimization code [PyDE](https://github.com/hpparvi/PyDE). ***Warning:*** if this keyword has been used for some of the variables, then all the unspecified variables will use the mid-point of the interval range as starting point

* ```Tcent```: if central time of transits have been measured, the file containing the measurements must be specified here and not in the ```Inputs```. In this way it is easier to associate the dataset to its planet.

*  ```Radius``` and ```Inclination``` are additional information on the planet used by ```TRADES``` (for the radius) and to determine the true mass of the planet (when _circular_ or _keplerian_ model are used.). This value for the inclination included here is not considered when performing the dynamical integration.

```yaml
Sinusoids:
  Prot: [5.00, 12.00]
  Seasons:
    0: [5000, 6250, 4, 4]
    1: [6250, 7000, 4, 4]
  Phase_coherence: 0
  Phase_dataset: False
```
Use this section  if the sinusoid model has been included for any of the datasets.
* ```Prot```: lower and upper limits on rotational period of the star (days)
* ```Seasons```: use this keyword to split the full observational range in several "Seasons". The first two values identify the range of each season (here BJD-2450000.0 following the convention in the datasets). The following integer values identify the number of harmonics  to be used in the sinusoidal model (0=only one sinusoid at the rotational period). The number of integer values on each row must be the same of the number of dataset, counted using the main keyword in the ```Inputs``` section.
* ```Phase_coherence```: turn on if you want all the harmonics to be at the same phase
* ```Phase_dataset```: turn on if you do not want any phase offset between photometry and other observations
*** Warning *** ```Phase_coherence``` and ```Phase_dataset``` have been introduced for testing purpose and they should be kept off for physical reasons

```yaml
emcee:
  Ngen: 2000
  Nsteps: 10000
  Nburn:   5000
  Npop_mult: 4
  Thin: 1
  Recenter_Bounds: True
```

```emcee``` parameters ([Foreman-Mackey et al. 2013](Foreman-Mackey2013), available [here](https://github.com/dfm/emcee)).
* ```Nburn```:  can be modified after running the sampler
* ```Npop_mult```: multiplicative factor, the minimum acceptable value is ```2```. The number of walkers is given by ```Npop_mult``` _x_ number of parameters.
* ```Thin```: thinning factor for posterior sampling. If no thinning is wanted, must be set to ```1``` and not to ```0```

### Data files

Datasets files must be built in the following way:

```terminal
$ more test1_RV.dat
6379.428786     34880.000         5.200   0   0  -1
6380.432890     35438.100         8.600   0   0  -1
6381.399091     34807.500         6.400   0   0  -1
```
For each row:

1. BJD (either BJD-2450000.0 or just BJD)
2. RV (in m/s)
3. error on RV (in m/s)
4. flag on jitter: -1 to not use any jitter term, increase +1 for each dataset with a different jitter
5. flag on offset: as for jitter
6. flag for linear trend: as for jitter

Last flag can be omitted if no linear trend due to the instrument is present. For trends of physical nature (i.e. RV trends due to the presence of a long-period binary star) the ```Curvature``` model should be used instead.

### Run the code
The main program is ```PyORBIT_V3_emcee.py```, you can start it by executing in a terminal:
```text
$ python PyORBIT_V3_emcee.py
```
followed by the ```yaml``` configuration file.
This program run the sampler and the outcome is store in the ```Outpout``` directory. Other samplers are currently under testing.

### Analyze the output

To analyze the outcome of the sampler run:
```text
$ python PyORBIT_V3_GetResults.py
```
followed by the ```yaml``` configuration file and the appropriate flag. Use the ```-h ``` flag to get a basic help on the feature of this script. The script does not overwrite the output of the sampler.


[Malavolta2016]: http://adsabs.harvard.edu/abs/2016A%26A...588A.118M
[Borsato2014]: http://adsabs.harvard.edu/abs/2014A%26A...571A..38B
[Deck2014]: http://adsabs.harvard.edu/abs/2014ApJ...787..132D
[Boisse2011]: http://adsabs.harvard.edu/abs/2011A%26A...528A...4B
[Foreman-Mackey2013]: http://adsabs.harvard.edu/abs/2013PASP..125..306F
