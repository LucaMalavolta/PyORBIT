from __future__ import print_function
#from pyorbit.classes.common import *
from pyorbit.classes.model_container_dynesty import ModelContainerDynesty
from pyorbit.subroutines.input_parser import yaml_parser, pars_input
from pyorbit.subroutines.io_subroutines import nested_sampling_save_to_cpickle, \
    nested_sampling_load_from_cpickle, nested_sampling_write_dummy_file, \
    dynesty_results_save_to_cpickle, dynesty_results_load_from_cpickle, \
    dynesty_results_maxevidence_save_to_cpickle, dynesty_results_maxevidence_load_from_cpickle

import pyorbit.subroutines.results_analysis as results_analysis
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

    mc = ModelContainerDynesty()
    pars_input(config_in, mc, input_datasets)

    mc.output_directory = './' + config_in['output'] + '/dynesty/'

    try:
        results = dynesty_results_load_from_cpickle(mc.output_directory)
        print('Dynesty results already saved in the respective directory, run PyORBIT_GetResults')
        if return_output:
            return mc
        else:
            return
    except FileNotFoundError:
        pass

    if not os.path.exists(mc.output_directory):
        os.makedirs(mc.output_directory)

    if mc.nested_sampling_parameters['shutdown_jitter']:
        'Jitter term not included for evidence calculation'
        print()
        for dataset_name, dataset in mc.dataset_dict.items():
            dataset.shutdown_jitter()

    mc.model_setup()
    mc.boundaries_setup()
    mc.initialize_logchi2()

    mc.starting_points_setup()

    results_analysis.results_resumen(mc, None, skip_theta=True)

    nthreads = mc.nested_sampling_parameters['nthreads']

    if 'nlive_mult' in mc.nested_sampling_parameters:
        nlive = mc.ndim * mc.nested_sampling_parameters['nlive_mult']
    else:
        nlive = mc.nested_sampling_parameters['nlive']


    print('Number of live points:', nlive)
    print('Number of threads:', nthreads)

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

    if mc.nested_sampling_parameters['pfrac'] > 1.0 or  mc.nested_sampling_parameters['pfrac'] < 0.0:
        print('Double run with 0/100 posterior/evidence split; then 100/0 posterior/evidence split')
        print()
        pfrac = 0.00

        with multiprocessing.Pool() as pool:

            # "Dynamic" nested sampling.
            print('Setting up the Dynamic Nested Sampling, posterior/evidence split = {0:4.3f}'.format(pfrac))
            print()
            dsampler_maxevidence = dynesty.DynamicNestedSampler(mc.dynesty_call,
                                                    mc.dynesty_priors,
                                                    mc.ndim,
                                                    nlive=nlive,
                                                    pool=pool,
                                                    queue_size=nthreads,
                                                    bound= mc.nested_sampling_parameters['bound'],
                                                    sample= mc.nested_sampling_parameters['sample'],
                                                    use_pool={
                                                        'prior_transform': False},
                                                    )
            print('Running Dynamic Nested Sampling')
            dsampler_maxevidence.run_nested(wt_kwargs={'pfrac': pfrac})
            print()
        print()

        print('Getting the results')
        results_maxevidence = dsampler_maxevidence.results

        #taken from dynesty/dynesty/results.py  but without the nlive point causing an error
        res = ("niter: {:d}\n"
                "ncall: {:d}\n"
                "eff(%): {:6.3f}\n"
                "logz: {:6.3f} +/- {:6.3f}"
                .format(results_maxevidence.niter, sum(results_maxevidence.ncall),
                        results_maxevidence.eff, results_maxevidence.logz[-1], results_maxevidence.logzerr[-1]))

        print()
        print('Summary\n=======\n'+res)

        print()
        print()

        """ A dummy file is created to let the cpulimit script to proceed with the next step"""
        dynesty_results_maxevidence_save_to_cpickle(mc.output_directory, results_maxevidence)

        pfrac = 1.00
    else:
        pfrac = mc.nested_sampling_parameters['pfrac']

    with multiprocessing.Pool() as pool:

        # "Dynamic" nested sampling.
        print('Setting up the Dynamic Nested Sampling, posterior/evidence split = {0:4.3f}'.format(pfrac))
        print()
        #dsampler = dynesty.DynamicNestedSampler(mc.dynesty_call,
        #                                        mc.dynesty_priors,
        #                                        mc.ndim,
        #                                        nlive=nlive,
        #                                        pool=pool,
        #                                        queue_size=nthreads,
        #                                        use_pool={
        #                                            'prior_transform': False},
        #                                        wt_kwargs={'pfrac': pfrac}
        #                                        )
        dsampler = dynesty.DynamicNestedSampler(mc.dynesty_call,
                                                mc.dynesty_priors,
                                                mc.ndim,
                                                nlive=nlive,
                                                pool=pool,
                                                #bounds='multi',
                                                bound= mc.nested_sampling_parameters['bound'],
                                                sample= mc.nested_sampling_parameters['sample'],
                                                queue_size=nthreads,
                                                use_pool={
                                                    'prior_transform': False},
                                                )

        print('Running Dynamic Nested Sampling')
        dsampler.run_nested(wt_kwargs={'pfrac': pfrac})
        print()
    print()

    print('Getting the results')
    results = dsampler.results

    #taken from dynesty/dynesty/results.py  but without the nlive point causing an error
    res = ("niter: {:d}\n"
            "ncall: {:d}\n"
            "eff(%): {:6.3f}\n"
            "logz: {:6.3f} +/- {:6.3f}"
            .format(results.niter, sum(results.ncall),
                    results.eff, results.logz[-1], results.logzerr[-1]))

    print()
    print('Summary\n=======\n'+res)

    print()


    """ A dummy file is created to let the cpulimit script to proceed with the next step"""
    nested_sampling_write_dummy_file(mc)
    nested_sampling_save_to_cpickle(mc)
    dynesty_results_save_to_cpickle(mc.output_directory, results)

    if return_output:
        return mc
    else:
        return
