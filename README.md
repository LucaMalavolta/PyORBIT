# PyORBIT

**`PyORBIT` v3.1 by Luca Malavolta - 2016**    

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
The ```Input``` section includes all the data files to be analyzed. Following th ePython standard, file must be listed starting from ```0``` without gaps.

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
      lN : 90.0
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
* ```Boundaries```: to be continued....

### Data files


### Run the code

### Analyze the output

The main program is in PyORBIT_MultiPhase_Sinfit.py, you can start it by executing in a terminal:

    *python PyORBIT_V2_Sinfit.py*

and then entering the name of the *yaml* file

V1 is the old version of the code (the one used in Malavolta et al. 2016)


#Help on config file
This is a brief guide on how to write the file requested by PyORBIT (see the example files)

    Name:      PyORBIT_test1
Name of the star

    Output:   PyORBIT_test1
Directory where all the output files are stored

    Inputs:  
List here all the input files, starting from 0

      0:
        File: test1_PH.dat
        Kind: Phot
        Models: ['sinusoids']

List the input file name, its kind (RV, Phot, BIS etc) and the models you want to use for this dataset (currently sinusoids, kepler or gaussian for GP)


    Planets:
List here the characteristics of each planet, starting from 0, as in the example below

    0:
      Orbit: kepler
      Boundaries:
        P: [2.0, 50.0]
        K: [25.0, 60.0]
      Fixed:
        e: 0.13
      Tcent: Tcent_0.dat

    Sinusoids:
Use this set of options if the sinusoid model is used for any of the datasets

      Prot: [5.00, 12.00]
      Seasons:
        0: [5000, 6250, 4, 4]
        1: [6250, 7000, 4, 4]

Prot: Lower and upper limits on rotational period of the star (days)
Seasons: Use this to split the full observational range in several "Seasons"
(here BJD-2450000.0 following the convention in the datasets). Then for each dataset the number of harmonics to be used to correct for activity.

      Phase_coherence: 0
Turn on if you want all the harmonics to be at the same phase

    Phase_dataset: False
Turn on if you do not want any phase offset between photometry and other observations
  Phase_coherence and Phase_synchro have been introduced for testing purpose and they should be kept to 0 for physical reasons

    emcee:
      Ngen: 5000
      Nsteps: 5000
      Nburn: 4000
      Npop_mult: 2
      Thin: 100
      Recenter_Bounds: True

*emcee* parameters. Pay attention to those ones
Npop_mult: number of walkers is Npop_mult x number of parameters
Thin: thinning factor for posterior sampling

#Help on dataset files

Datasets files must be built in the following way:

_Py_OC102_RV.dat:_

    6379.428786     34880.000000         5.200000   0   0  -1   0
    6380.432890     35438.100000         8.600000   0   0  -1   0
    6381.399091     34807.500000         6.400000   0   0  -1   0

_For each row:

1. BJD (either BJD-2450000.0 or just BJD)

2. RV (in m/s)

3. error on RV (in m/s)

4. flag on jitter: -1 to not use any jitter term, increase +1 for each dataset with a different jitter

5. flag on offset: as for jitter

6. flag for linear trend: as for jitter

Last flag can be omitted if no linear trend due to the instrument is present


[Malavolta2016]: http://adsabs.harvard.edu/abs/2016A%26A...588A.118M
[Borsato2014]: http://adsabs.harvard.edu/abs/2014A%26A...571A..38B
[Deck2014]: http://adsabs.harvard.edu/abs/2014ApJ...787..132D
[Boisse2011]: http://adsabs.harvard.edu/abs/2011A%26A...528A...4B
