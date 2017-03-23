from classes.model_container import ModelContainerPolyChord
from classes.input_parser import yaml_parser
import numpy as np
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

parser = argparse.ArgumentParser(prog='PyORBIT_V4_PolyChord.py', description='PyDE+emcee runner')
parser.add_argument('config_file', type=str, nargs=1, help='config file')


args = parser.parse_args()

file_conf = args.config_file[0]

mc = ModelContainerPolyChord()

yaml_parser(file_conf, mc)

if mc.polychord_parameters['shutdown_jitter']:
    for dataset in mc.dataset_list:
        dataset.shutdown_jitter()

mc.initialize_model()

dir_output = './' + mc.planet_name + '/'
os.system("mkdir -p " +dir_output+ mc.polychord_parameters['base_dir'] + "/clusters")

if bool(mc.pcv.dynamical):
        mc.dynamical_model.prepare(mc, mc.pcv)

print 'Dimensions = ', mc.ndim
print '   '
print 'Variable list:', mc.variable_list
print
print 'Variable bounds:', mc.bounds
print


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


# number of dimensions our problem has
print '--->', mc.pam_names

parameters = mc.pam_names
num_repeats = mc.polychord_parameters['num_repeats_mult'] * mc.ndim

nlive = mc.ndim * mc.polychord_parameters['nlive_mult']
if 'nlive' in mc.polychord_parameters:
    nlive = mc.polychord_parameters['nlive']

os.chdir(dir_output)

PolyChord.mpi_notification()
PolyChord.run_nested_sampling(mc.polychord_call, nDims=mc.ndim, nDerived=0,
                              feedback=mc.polychord_parameters['feedback'],
                              base_dir=mc.polychord_parameters['base_dir'],
                              precision_criterion=mc.polychord_parameters['precision_criterion'],
                              max_ndead=mc.polychord_parameters['max_ndead'],
                              boost_posterior=mc.polychord_parameters['boost_posterior'],
                              read_resume=mc.polychord_parameters['read_resume'],
                              file_root=mc.planet_name,
                              prior=mc.polychord_priors, nlive=nlive, num_repeats=num_repeats)

print
print 'PolyChord COMPLETED'
print