from __future__ import print_function
#from pyorbit.classes.common import *
from pyorbit.classes.model_container_dynesty import ModelContainerDynesty
from pyorbit.classes.input_parser import yaml_parser, pars_input
from pyorbit.classes.io_subroutines import nested_sampling_save_to_cpickle, \
    nested_sampling_load_from_cpickle, nested_sampling_create_dummy_file, \
    dynesty_results_save_to_cpickle, dynesty_results_load_from_cpickle
import pyorbit.classes.results_analysis as results_analysis
import os
import sys
import argparse
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt

__all__ = ["pyorbit_dynesty", "yaml_parser"]

""" 
def show(filepath):
    # open the output (pdf) file for the user
    if os.name == 'mac': subprocess.call(('open', filepath))
    elif os.name == 'nt': os.startfile(filepath)
"""


def pyorbit_dynesty(config_in, input_datasets=None, return_output=None):

    output_directory = './' + config_in['output'] + '/dynesty/'

    mc = ModelContainerDynesty()
    pars_input(config_in, mc, input_datasets)

    if mc.nested_sampling_parameters['shutdown_jitter']:
        'Jitter term not included for evidence calculation'
        print()
        for dataset_name, dataset in mc.dataset_dict.items():
            dataset.shutdown_jitter()

    mc.model_setup()
    mc.create_variables_bounds()
    mc.initialize_logchi2()

    mc.create_starting_point()

    results_analysis.results_resumen(mc, None, skip_theta=True)

    mc.output_directory = output_directory

    if 'nlive' in mc.nested_sampling_parameters:
        nlive = mc.nested_sampling_parameters['nlive']
    elif 'nlive_mult' in mc.nested_sampling_parameters:
        nlive = mc.ndim * mc.nested_sampling_parameters['nlive_mult']

    print('Number of live points:', nlive) 

    print()
    print('Reference Time Tref: ', mc.Tref)
    print()
    print('*************************************************************')
    print()

    try:
        import dynesty
    except ImportError:
        print("ERROR: dynesty not installed, this will not work")
        quit()

    # "Standard" nested sampling.
    #print('Setting up the Standard Nested Sampling')
    #sampler = dynesty.NestedSampler(mc.dynesty_call, mc.dynesty_priors, mc.ndim)
    #print('Running Nested Sampling')
    # sampler.run_nested()
    #print('Getting the results')
    #results = sampler.results
    # print()

    with multiprocessing.Pool() as pool:

        # "Dynamic" nested sampling.
        print('Setting up the Dynamic Nested Sampling')
        dsampler = dynesty.DynamicNestedSampler(mc.dynesty_call,
                                                mc.dynesty_priors,
                                                mc.ndim,
                                                nlive=nlive,
                                                pool=pool,
                                                queue_size=16,
                                                use_pool={
                                                    'prior_transform': False},
                                                wt_kwargs={'pfrac': 0.0}
                                                )
        print('Running Dynamic Nested Sampling')
        dsampler.run_nested()

    print('Getting the results')
    results = dsampler.results

    #taken from dynesty/dynesty/results.py  but without the nlive point causing an error
    res = ("niter: {:d}\n"
            "ncall: {:d}\n"
            "eff(%): {:6.3f}\n"
            "logz: {:6.3f} +/- {:6.3f}"
            .format(results.niter, sum(results.ncall),
                    results.eff, results.logz[-1], results.logzerr[-1]))

    print('Summary\n=======\n'+res)

    print()
    

    """ A dummy file is created to let the cpulimit script to proceed with the next step"""
    nested_sampling_create_dummy_file(mc)
    nested_sampling_save_to_cpickle(results)

    if return_output:
        return mc
    else:
        return
