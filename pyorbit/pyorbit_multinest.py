from classes.common import *
from classes.model_container_multinest import ModelContainerMultiNest
from classes.input_parser import yaml_parser, pars_input
from classes.io_subroutines import nested_sampling_save_to_cpickle, nested_sampling_load_from_cpickle, nested_sampling_create_dummy_file
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


    output_directory = './' + config_in['output'] + '/multinest/'
    reloaded_mc = False

    try:
        mc = nested_sampling_load_from_cpickle(output_directory, prefix='')
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

        if mc.nested_sampling_parameters['shutdown_jitter']:
            for dataset in mc.dataset_dict.itervalues():
                dataset.shutdown_jitter()

        mc.model_setup()
        mc.create_variables_bounds()
        mc.initialize_logchi2()

        mc.create_starting_point()

        mc.results_resumen(None, skip_theta=True)

        mc.output_directory = output_directory

    os.system("mkdir -p " + output_directory + mc.nested_sampling_parameters['base_dir'] + "/clusters")
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
        export LD_PRELOAD=/opt/local/lib/openmpi/lib/libmpi.so:$LD_PRELOAD
        mpirun -np 4 python run_PyPolyChord.py
    
    '''

    # parameters = mc.get_theta_dictionary()

    if 'nlive' in mc.nested_sampling_parameters:
        nlive = mc.nested_sampling_parameters['nlive']
    elif 'nlive_mult' in mc.nested_sampling_parameters:
        nlive = mc.ndim * mc.nested_sampling_parameters['nlive_mult']

    print ' Sampling efficiency: ', mc.nested_sampling_parameters['sampling_efficiency']

    import pymultinest
    mnest_kwargs = dict(n_live_points=nlive, outputfiles_basename=output_directory + './')

    for key_name, key_value in mc.nested_sampling_parameters.items():
        if key_name in mc.pymultinest_signature:
            mnest_kwargs[key_name] = key_value

    pymultinest.run(LogLikelihood= mc.multinest_call, Prior=mc.multinest_priors, n_dims=mc.ndim, **mnest_kwargs)

    nested_sampling_save_to_cpickle(mc)

    print
    print 'MultiNest COMPLETED'
    print

    """ A dummy file is created to let the cpulimit script to proceed with the next step"""
    nested_sampling_create_dummy_file(mc)

    if return_output:
        return mc
    else:
        return
