# PyORBIT
PyORBIT_MultiPhase*.py files are from the old version of the code (the one used in Malavolta et al. 2016)\
The main program is in PyORBIT_MultiPhase_Sinfit.py

You can start it by executing in a terminal:

*python PyORBIT_MultiPhase_Sinfit.py*

and then entering the name of the *config* file

I'm working on a new version (PyORBIT_V2*.py) with the same functionalities but better coding style.


#Help on config file
This is a brief guide on how to write the file requested at PyORBIT  file (see OC102_HT_Act_2pEE.config in the example directory) for V1 of PyORBIT:
(OC102 is just the GAPS name of Pr0211)

Name      OC102

  _Name of the star_

Input_RV  Py_OC102_RV.dat

  _File with all the radial velocities from different instruments_

Input_P1  Py_OC102_PH_B.dat

  _First photometry dataset_

Input_P2  Py_OC102_PH_R.dat

  _Second photometry dataset_

Output    OC102_HT_Act_2pEE

  _String pre-pended to output files_

Nplanets  2

  _Number of planets in the model_

Planet1   5    2.14     2.15  300.0  317.0    0.00 0.4

Planet2   5  900.0   25000.0  126.0  150.0    0.20 1.0

  _For each planet:

    _Number of parameters (3 for circular orbit, 5 for eccentric one)_

    _Lower and upper limits on Period (in days)_

    _Lower and upper limits on K semi-amplitude (m/s)_

    _Lower and upper limits on eccentricity (not read if first parameter is 3)_

Prot      4    7.83     8.00

  _4 to include sinusoidal modeling of activity , 0 otherwise_

  _Lower and upper limits on rotational period of the star (days)_

Nperiods  3

  _Use this to split the full observational range in several "periods"_

Period0  6300 6600

Period1  6600 6900

Period2  6900 9000

  _Lower and upper limits of each temporal "period" (here BJD-2450000.0 following the convention in the datasets_

Spec_dataset 2 2 1

  _The first parameter is the number of spectroscopic datasets in Input_RV_

  _Then for each dataset the number of harmonics to be used to correct for activity_

Phase_coherence 0

  _Turn on if you want all the harmonics to be at the same phase_

Phase_synchro   0

  _Turn on if you do not want any phase offset between photometry and other observations_

  _Phase_coherence and Phase_synchro have been introduced for testing purpose and they should be kept to 0 for physical reasons_

Ngen      16000

Nsteps   400000

Nburn    200000

_emcee parameters_

Npop_mult     8

_number of walkers is Npop_mult x number of parameters_

Thin        100

_Thinning factor for posterior sampling_

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

7. flag on activity (RV only!) as for jitter_
