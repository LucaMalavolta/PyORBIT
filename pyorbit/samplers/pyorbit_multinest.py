from __future__ import print_function
from pyorbit.subroutines.common import *
from pyorbit.classes.model_container_multinest import ModelContainerMultiNest
from pyorbit.subroutines.input_parser import yaml_parser, pars_input
from pyorbit.subroutines.io_subroutines import nested_sampling_save_to_cpickle, nested_sampling_load_from_cpickle, nested_sampling_write_dummy_file
import pyorbit.subroutines.results_analysis as results_analysis
import os
import sys
import argparse

__all__ = ["pyorbit_multinest"]

"""
def show(filepath):
    # open the output (pdf) file for the user
    if os.name == 'mac': subprocess.call(('open', filepath))
    elif os.name == 'nt': os.startfile(filepath)
"""


def pyorbit_multinest(config_in, input_datasets=None, return_output=None):

    output_directory = './' + config_in['output'] + '/multinest/'

    mc = ModelContainerMultiNest()
    pars_input(config_in, mc, input_datasets)

    if mc.nested_sampling_parameters['shutdown_jitter']:
        for dataset_name, dataset in mc.dataset_dict.items():
            dataset.shutdown_jitter()

    mc.model_setup()
    mc.boundaries_setup()
    mc.initialize_logchi2()

    mc.starting_points_setup()

    results_analysis.print_bayesian_info(mc)

    sys.stdout.flush()

    mc.output_directory = output_directory

    os.system("mkdir -p " + output_directory + mc.nested_sampling_parameters['base_dir'] + "/clusters")
    #os.system("mkdir -p " +polychord_dir_output + "chains/clusters")

    print()
    print('Reference Time Tref: ', mc.Tref)
    print()
    print('*************************************************************')
    print()

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

    print(' Sampling efficiency: ', mc.nested_sampling_parameters['sampling_efficiency'])
    print(' N live points:', nlive)

    import pymultinest

    mnest_kwargs = dict(n_live_points=nlive, outputfiles_basename=output_directory + './')

    for key_name, key_value in mc.nested_sampling_parameters.items():
        if key_name in mc.pymultinest_signature:
            mnest_kwargs[key_name] = key_value

    print('Including priors to log-likelihood calculation (must be False):', mc.include_priors)

    sys.stdout.flush()

    pymultinest.run(LogLikelihood=mc.multinest_call, Prior=mc.multinest_priors, n_dims=mc.ndim, **mnest_kwargs)

    nested_sampling_save_to_cpickle(mc)

    analyzer = pymultinest.Analyzer(mc.ndim, outputfiles_basename=output_directory)
    stats = analyzer.get_stats()
    samples = analyzer.get_equal_weighted_posterior()[:, :-1]

    result = dict(logZ=stats['nested sampling global log-evidence'],
         logZerr=stats['nested sampling global log-evidence error'],
         samples=samples,
         )

    nested_sampling_save_to_cpickle(mc, 'result')

    sys.stdout.flush()
    print()
    print('MultiNest COMPLETED')
    print()

    #result = pymultinest.solve(LogLikelihood=mc.multinest_call, Prior=mc.multinest_priors,
    #                           n_dims=mc.ndim, outputfiles_basename=output_directory + './',
    #                           n_live_points=1000, sampling_efficiency=0.3, multimodal=True,
    #                           verbose=True, resume=True)

    print('evidence: %(logZ).1f +- %(logZerr).1f' % result)
    print('evidenze in base-10 log: ',result['logZ']//np.log(10.00), result['logZerr']//np.log(10.00))

    if return_output:
        return mc
    else:
        return
