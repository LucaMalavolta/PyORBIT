from classes.common import *
from classes.model_container import ModelContainerPolyChord
from classes.input_parser import yaml_parser, pars_input
from classes.io_subroutines import polychord_save_to_cpickle, polychord_load_from_cpickle
import os
import sys
import argparse

__all__ = ["pyorbit_polychord", "yaml_parser"]

""" 
def show(filepath):
    # open the output (pdf) file for the user
    if os.name == 'mac': subprocess.call(('open', filepath))
    elif os.name == 'nt': os.startfile(filepath)
"""


def pyorbit_polychord(config_in, input_datasets=None, pl_version=''):


    polychord_dir_output = './' + config_in['output'] + '/polychord' + pl_version + '/'

    reloaded_mc = False



    try:
        mc = polychord_load_from_cpickle(polychord_dir_output, prefix='')
        reloaded_mc = True
    except:
        pass

    if reloaded_mc:
        mc.model_setup()
        mc.initialize_logchi2()
        #mc.results_resumen(flatchain)
    else:
        mc = ModelContainerPolyChord()
        pars_input(config_in, mc, input_datasets)

        #if mc.polychord_parameters['shutdown_jitter']:
        #    for dataset in mc.dataset_dict.itervalues():
        #        dataset.shutdown_jitter()

        mc.model_setup()
        mc.create_variables_bounds()
        mc.initialize_logchi2()

        mc.create_starting_point()

        mc.results_resumen(None, skip_theta=True)

        mc.polychord_dir_output = polychord_dir_output



    os.system("mkdir -p " + polychord_dir_output + mc.polychord_parameters['base_dir'] + "/clusters")
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

    #os.chdir(polychord_dir_output)

    if pl_version == '':

        import PyPolyChord
        from PyPolyChord.settings import PolyChordSettings

        settings = PolyChordSettings(nDims=mc.ndim, nDerived=0)
        settings.feedback=mc.polychord_parameters['feedback']
        settings.base_dir = polychord_dir_output
        settings.precision_criterion=mc.polychord_parameters['precision_criterion']
        settings.max_ndead=mc.polychord_parameters['max_ndead']
        settings.boost_posterior=mc.polychord_parameters['boost_posterior']
        settings.read_resume=mc.polychord_parameters['read_resume']
        settings.file_root = 'pyorbit'

        settings.nlive=nlive
        settings.num_repeats=num_repeats
        settings.do_clustering = True

        output = PyPolyChord.run_polychord(mc.polychord_call, nDims=mc.ndim, nDerived=0, settings=settings,
                                           prior=mc.polychord_priors)

        paramnames = [('p%i' % i, r'\theta_%i' % i) for i in range(mc.ndim)]
        paramnames += [('r*', 'r')]
        output.make_paramnames_files(paramnames)
        quit()

        print output

    elif pl_version=='v1.9':
        if os.path.isdir('/Users/malavolta/Astro/CODE/pyorbit/'):
            sys.path.insert(0, '/Users/malavolta/Astro/CODE/others/PolyChord_1.9/')
        else:
            sys.path.insert(0, '/home/malavolta/CODE/others/PolyChord_1.9/')
        import PyPolyChord.PyPolyChord as PyPolyChord_1_9


        PyPolyChord_1_9.mpi_notification()
        PyPolyChord_1_9.run_nested_sampling(mc.polychord_call, nDims=mc.ndim, nDerived=0,
                                      feedback=mc.polychord_parameters['feedback'],
                                      base_dir = polychord_dir_output,
                                      precision_criterion=mc.polychord_parameters['precision_criterion'],
                                      max_ndead=mc.polychord_parameters['max_ndead'],
                                      boost_posterior=mc.polychord_parameters['boost_posterior'],
                                      read_resume=mc.polychord_parameters['read_resume'],
                                      file_root='pyorbit',
                                      prior=mc.polychord_priors, nlive=nlive, num_repeats=num_repeats)

    polychord_save_to_cpickle(mc)

    print
    print 'PolyChord COMPLETED'
    print


