from classes.common import *
from classes.model_container_multinest import ModelContainerMultiNest
from classes.input_parser import yaml_parser, pars_input
from classes.io_subroutines import polychord_save_to_cpickle, polychord_load_from_cpickle, polychord_create_dummy_file
import os
import sys
import argparse

__all__ = ["pyorbit_multinest", "yaml_parser"]

""" 
def show(filepath):
    # open the output (pdf) file for the user
    if os.name == 'mac': subprocess.call(('open', filepath))
    elif os.name == 'nt': os.startfile(filepath)
"""


def pyorbit_multinest(config_in, input_datasets=None, return_output=None):


    multinest_dir_output = './' + config_in['output'] + '/multinest/'

    reloaded_mc = False


    try:
        mc = polychord_load_from_cpickle(multinest_dir_output, prefix='')
        reloaded_mc = True
    except:
        pass

    if reloaded_mc:
        mc.model_setup()
        mc.initialize_logchi2()
        #mc.results_resumen(flatchain)
    else:
        mc = ModelContainerMultiNest()
        pars_input(config_in, mc, input_datasets)

        if mc.polychord_parameters['shutdown_jitter']:
            for dataset in mc.dataset_dict.itervalues():
                dataset.shutdown_jitter()

        mc.model_setup()
        mc.create_variables_bounds()
        mc.initialize_logchi2()

        mc.create_starting_point()

        mc.results_resumen(None, skip_theta=True)

        mc.polychord_dir_output = multinest_dir_output

    os.system("mkdir -p " + multinest_dir_output + mc.polychord_parameters['base_dir'] + "/clusters")
    #os.system("mkdir -p " +polychord_dir_output + "chains/clusters")

    print
    print 'Reference Time Tref: ', mc.Tref
    print
    print '*************************************************************'
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

    parameters = mc.get_theta_dictionary()

    num_repeats = mc.polychord_parameters['num_repeats_mult'] * mc.ndim

    nlive = mc.ndim * mc.polychord_parameters['nlive_mult']
    if 'nlive' in mc.polychord_parameters:
        nlive = mc.polychord_parameters['nlive']

    #os.chdir(multinest_dir_output)

    print ' Sampling efficiency: ', mc.polychord_parameters['sampling_efficiency']

    import pymultinest
    mnest_kwargs = dict(n_live_points=nlive, outputfiles_basename=multinest_dir_output + './',
                        sampling_efficiency=mc.polychord_parameters['sampling_efficiency'],
                        verbose=True)

    pymultinest.run(LogLikelihood= mc.multinest_call, Prior=mc.multinest_priors, n_dims=mc.ndim, **mnest_kwargs)

    polychord_save_to_cpickle(mc)

    print
    print 'MultiNest COMPLETED'
    print

    """ A dummy file is created to let the cpulimit script to proceed with the next step"""
    polychord_create_dummy_file(mc)

    if return_output:
        return mc
    else:
        return
