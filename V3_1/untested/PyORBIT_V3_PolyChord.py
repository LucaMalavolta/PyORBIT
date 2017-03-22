from PyORBIT_V3_Classes import *
import numpy as np
from numpy import pi,log,sqrt
import os
import sys
import argparse
import json
import threading, subprocess
if os.path.isdir('/Users/malavolta/Astro/CODE/'):
    sys.path.insert(0, '/Users/malavolta/Astro/CODE/others/PolyChord/')
else:
    sys.path.insert(0, '/home/malavolta/CODE/others/PolyChord/')
import PyPolyChord.PyPolyChord as PolyChord


def show(filepath):
    """ open the output (pdf) file for the user """
    if os.name == 'mac': subprocess.call(('open', filepath))
    elif os.name == 'nt': os.startfile(filepath)

parser = argparse.ArgumentParser(prog='PyORBIT_V3_PolyChord.py', description='PyDE+emcee runner')
parser.add_argument('config_file', type=str, nargs=1, help='config file')


args = parser.parse_args()

file_conf = args.config_file[0]

mc = ModelContainerPolyChord()

yaml_parser(file_conf, mc)


dir_output = './' + mc.planet_name + '/'
if not os.path.exists(dir_output):
    os.makedirs(dir_output)
if not os.path.exists(dir_output+'clusters'):
    os.makedirs(dir_output+'clusters')

mc.create_bounds()
print 'Dimensions = ', mc.ndim
print '   '
print 'Variable list:', mc.variable_list
print
print 'Variable bounds:', mc.bounds
print

mc.nDerived = 0

'''
    On Linux system (BASH):
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/lib
    export LD_PRELOAD=/usr/lib/openmpi/lib/libmpi.so:$LD_PRELOAD
    export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libgfortran.so.3
    mpirun -np 4 python run_PyPolyChord.py

    on Mac:
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/lib
    export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libgfortran.so.3
    export LD_PRELOAD=/opt/local/lib/openmpi//lib/libmpi.so:$LD_PRELOAD
    mpirun -np 4 python run_PyPolyChord.py

'''
PolyChord.mpi_notification()
PolyChord.run_nested_sampling(mc.polychord_call, mc.ndim, mc.nDerived, prior=mc.polychord_priors, base_dir=dir_output)
