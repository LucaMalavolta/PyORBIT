# PyORBIT
PyORBIT_MultiPhase*.py files are from the old version of the code (the one used in Malavolta et al. 2016)
The main program is in PyORBIT_MultiPhase_Sinfit.py
You can start it by executing in a terminal:
python PyORBIT_MultiPhase_Sinfit.py
and then entering the name of the *config* file
I'm working on a new version (PyORBIT_V2*.py) with the same functionalities but better coding style.


#Help on config file
This is a brief guide on how to write the file requested at PyORBIT  file (see OC102_HT_Act_2pEE.config in the example directory) for V1 of PyORBIT:
(OC102 is just the GAPS name of Pr0211)

Name      OC102
__  Name of the star__
Input_RV  Py_OC102_RV.dat
__  File with all the radial velocities from different instruments__
Input_P1  Py_OC102_PH_B.dat
__  First photometry dataset__
Input_P2  Py_OC102_PH_R.dat
__  Second photometry dataset__
Output    OC102_HT_Act_2pEE
__  String pre-pended to output files__
Nplanets  2
__  Number of planets in the model
Planet1   5    2.14     2.15  300.0  317.0    0.00 0.4
Planet2   5  900.0   25000.0  126.0  150.0    0.20 1.0
__ For each planet:
__  Number of parameters (3 for circular orbit, 5 for eccentric one)__
__  Lower and upper limits on Period (in days)__
__  Lower and upper limits on K semi-amplitude (m/s)__
__  Lower and upper limits on eccentricity (not read if first parameter is 3)__
Prot      4    7.83     8.00
__  4 to include sinusoidal modeling of activity , 0 otherwise__
__  Lower and upper limits on rotational period of the star (days)__
Nperiods  3
__ Use this to split the full observational range in several "periods"__
Period0  6300 6600
Period1  6600 6900
Period2  6900 9000
__ Lower and upper limits of each temporal "period" (here BJD-2450000.0 following the convention in the datasets__
Spec_dataset 2 2 1
__ The first parameter is the number of spectroscopic datasets in Input_RV__
__ Then for each dataset the number of harmonics to be used to correct for activity__
Phase_coherence 0
__ Turn on if you want all the harmonics to be at the same phase__
Phase_synchro   0
__ Turn on if you do not want any phase offset between photometry and other observations__
__ Phase_coherence and Phase_synchro have been introduced for testing purpose and they should be kept to 0 for physical reasons__
Ngen      16000
Nsteps   400000
Nburn    200000
__  emcee parameters__
Npop_mult     8
__ number of walkers is Npop_mult x number of parameters__
Thin        100
__ Thinning factor for posterior sampling__


#Help on dataset files
Datasets files must be built in the following way:
Py_OC102_RV.dat:
6379.428786     34880.000000         5.200000   0   0  -1   0
6380.432890     35438.100000         8.600000   0   0  -1   0
6381.399091     34807.500000         6.400000   0   0  -1   0
__For each row:
BJD (either BJD-2450000.0 or just BJD)
RV (in m/s)
error on RV (in m/s)
flag on jitter: -1 to not use any jitter term, increase +1 for each dataset with a different jitter
flag on offset: as for jitter
flag for linear trend: as for jitter
flag on activity (RV only!) as for jitter__
