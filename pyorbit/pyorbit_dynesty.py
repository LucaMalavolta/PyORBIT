from __future__ import print_function
#from pyorbit.classes.common import *
from pyorbit.classes.model_container_dynesty import ModelContainerDynesty
from pyorbit.classes.input_parser import yaml_parser, pars_input
from pyorbit.classes.io_subroutines import nested_sampling_save_to_cpickle, \
    nested_sampling_load_from_cpickle, nested_sampling_create_dummy_file
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
        for dataset in mc.dataset_dict.itervalues():
            dataset.shutdown_jitter()

    mc.model_setup()
    mc.create_variables_bounds()
    mc.initialize_logchi2()

    mc.create_starting_point()

    results_analysis.results_resumen(mc, None, skip_theta=True)

    mc.output_directory = output_directory

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
                                                pool=pool,
                                                queue_size=16,
                                                use_pool={
                                                    'prior_transform': False}
                                                )
        print('Running Dynamic Nested Sampling')
        dsampler.run_nested()

    print('Getting the results')
    results = dsampler.results
    print('Results: ', results)
    print()
    from dynesty import plotting as dyplot

    # Plot a summary of the run.
    print('Plot a summary of the run.')
    rfig, raxes = dyplot.runplot(results)
    rfig.savefig('plot01.pdf', bbox_inches='tight', dpi=300)
    plt.close(rfig)

    # Plot traces and 1-D marginalized posteriors.
    print('Plot traces and 1-D marginalized posteriors.')
    tfig, taxes = dyplot.traceplot(results)
    tfig.savefig('plot02.pdf', bbox_inches='tight', dpi=300)
    plt.close(tfig)

    # Plot the 2-D marginalized posteriors.
    print('Plot the 2-D marginalized posteriors.')
    cfig, caxes = dyplot.cornerplot(results)
    cfig.savefig('plot03.pdf', bbox_inches='tight', dpi=300)
    plt.close(cfig)

    from dynesty import utils as dyfunc

    # Extract sampling results.
    samples = results.samples  # samples
    weights = np.exp(results.logwt - results.logz[-1])  # normalized weights

    # Compute 5%-95% quantiles.
    quantiles = dyfunc.quantile(samples, [0.05, 0.95], weights=weights)

    # Compute weighted mean and covariance.
    mean, cov = dyfunc.mean_and_cov(samples, weights)

    # Resample weighted samples.
    samples_equal = dyfunc.resample_equal(samples, weights)

    # Generate a new set of results with statistical+sampling uncertainties.
    results_sim = dyfunc.simulate_run(results)

    """ A dummy file is created to let the cpulimit script to proceed with the next step"""
    nested_sampling_create_dummy_file(mc)

    if return_output:
        return mc
    else:
        return
