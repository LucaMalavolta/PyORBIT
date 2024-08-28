from __future__ import print_function
#from pyorbit.classes.common import *
from pyorbit.classes.model_container_nestle import ModelContainerNestle
from pyorbit.subroutines.input_parser import yaml_parser, pars_input
from pyorbit.subroutines.io_subroutines import nested_sampling_save_to_cpickle, \
    nested_sampling_load_from_cpickle, nested_sampling_write_dummy_file, \
    nestle_results_save_to_cpickle, nestle_results_load_from_cpickle, \
    nestle_results_maxevidence_save_to_cpickle, nestle_results_maxevidence_load_from_cpickle

import pyorbit.subroutines.results_analysis as results_analysis
import os
import sys
import argparse
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt

__all__ = ["pyorbit_nestle"]

"""
def show(filepath):
    # open the output (pdf) file for the user
    if os.name == 'mac': subprocess.call(('open', filepath))
    elif os.name == 'nt': os.startfile(filepath)
"""

#http://mattpitkin.github.io/samplers-demo/pages/samplers-samplers-everywhere/#Nestle

def pyorbit_nestle(config_in, input_datasets=None, return_output=None):

    mc = ModelContainerNestle()
    pars_input(config_in, mc, input_datasets)

    mc.output_directory = './' + config_in['output'] + '/nestle/'

    try:
        results = nestle_results_load_from_cpickle(mc.output_directory)
        print('Nestle results already saved in the respective directory, run pyorbit_results')
        if return_output:
            return mc
        else:
            return
    except FileNotFoundError:
        pass

    if not os.path.exists(mc.output_directory):
        os.makedirs(mc.output_directory)

    if mc.nested_sampling_parameters['shutdown_jitter']:
        'Jitter term not included for evidence calculation - ARE YOU SURE???'
        print()
        for dataset_name, dataset in mc.dataset_dict.items():
            dataset.shutdown_jitter()

    os.environ["OMP_NUM_THREADS"] = "1"
    try:
        num_threads = int(config_in['parameters'].get('cpu_threads', "1"))
    except:
        num_threads = multiprocessing.cpu_count()-1

    mc.model_setup()
    mc.boundaries_setup()
    mc.initialize_logchi2()

    mc.starting_points_setup()

    results_analysis.print_bayesian_info(mc)

    nthreads = mc.nested_sampling_parameters['nthreads']

    if 'nlive_mult' in mc.nested_sampling_parameters:
        nlive = mc.ndim * mc.nested_sampling_parameters['nlive_mult']
    else:
        nlive = mc.nested_sampling_parameters['nlive']


    print('*************************************************************')
    print()
    try:
        import nestle
        print('Nestle version: {}'.format(nestle.__version__))
        print()
    except (ModuleNotFoundError,ImportError):
        print("ERROR: nestle not installed, this will not work")
        quit()


    print('Number of live points:', nlive)
    print('Number of threads:', nthreads)

    print('Reference Time Tref: ', mc.Tref)
    print()
    print('*************************************************************')
    print()

    # "Standard" nested sampling.
    #print('Setting up the Standard Nested Sampling')
    #sampler = dynesty.NestedSampler(mc.dynesty_call, mc.dynesty_priors, mc.ndim)
    #print('Running Nested Sampling')
    # sampler.run_nested()
    #print('Getting the results')
    #results = sampler.results
    # print()


    results = nestle.sample(mc.nestle_call,
                        mc.nestle_priors,
                        mc.ndim,
                        method=mc.nested_sampling_parameters['bound'],
                        npoints=nlive,
                        callback=nestle.print_progress,
                        dlogz=mc.nested_sampling_parameters['dlogz'])



    #print('Getting the results')
    #print(results.summary())
    #print()

    #taken from dynesty/dynesty/results.py  but without the nlive point causing an error
    #res = ("niter: {:d}\n"
    #        "ncall: {:d}\n"
    #        "logz: {:6.3f} +/- {:6.3f}"
    #        .format(results.niter, sum(results.ncall),
    #                results.logz[-1], results.logzerr[-1]))
    #
    #print()
    #print('Summary\n=======\n'+res)

    #print()


    """ A dummy file is created to let the cpulimit script to proceed with the next step"""
    nested_sampling_write_dummy_file(mc)
    nested_sampling_save_to_cpickle(mc)
    nestle_results_save_to_cpickle(mc.output_directory, results)

    if return_output:
        return mc
    else:
        return
