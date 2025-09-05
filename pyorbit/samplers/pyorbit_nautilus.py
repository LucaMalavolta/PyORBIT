from __future__ import print_function
from __future__ import print_function
#from pyorbit.classes.common import *
from pyorbit.classes.model_container_nautilus import ModelContainerNautilus
from pyorbit.subroutines.input_parser import yaml_parser, pars_input
from pyorbit.subroutines.io_subroutines import nested_sampling_save_to_cpickle, \
    nested_sampling_load_from_cpickle, nested_sampling_write_dummy_file, \
    nautilus_results_save_to_cpickle, nautilus_results_load_from_cpickle, \
    nautilus_sampler_save_to_cpickle, nautilus_sampler_load_from_cpickle

import pyorbit.subroutines.results_analysis as results_analysis
import os
import sys
import argparse
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt
from pyorbit.subroutines.common import np, nested_sampling_prior_compute

__all__ = ["pyorbit_nautilus"]

"""
def show(filepath):
    # open the output (pdf) file for the user
    if os.name == 'mac': subprocess.call(('open', filepath))
    elif os.name == 'nt': os.startfile(filepath)
"""


def pyorbit_nautilus(config_in, input_datasets=None, return_output=None, run_nested=True):

    multiprocessing.set_start_method('fork')
    os.environ["OMP_NUM_THREADS"] = "1"

    reloaded_nautilus = False
    output_directory = './' + config_in['output'] + '/nautilus/'
    save_checkpoint = output_directory + 'checkpoint.hdf5'
    save_sampler = output_directory + 'sample.hdf5'

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    try:
        mc = nested_sampling_load_from_cpickle(output_directory)
        reloaded_nautilus = True

        pars_input(config_in, mc, input_datasets, reload_emcee=True)
        mc.output_directory = output_directory

    except:
        mc = ModelContainerNautilus()

        pars_input(config_in, mc, input_datasets)
        mc.output_directory = output_directory

    try:
        results = nautilus_results_load_from_cpickle(mc.output_directory)
        print('nautilus results already saved in the respective directory, run pyorbit_results')
        if return_output:
            return mc
        else:
            return
    except FileNotFoundError:
        pass


    if mc.nested_sampling_parameters['shutdown_jitter']:
        'Jitter term not included for evidence calculation'
        print()
        for dataset_name, dataset in mc.dataset_dict.items():
            dataset.shutdown_jitter()

    os.environ["OMP_NUM_THREADS"] = "1"
    try:
        num_threads = int(config_in['parameters'].get('cpu_threads',  multiprocessing.cpu_count()))
    except:
        print(" Something happened when trying to setup multiprocessing, switching back to 1 CPU")
        num_threads = 1

    if reloaded_nautilus:
        mc.model_setup()
        #mc.boundaries_setup()
        mc.initialize_logchi2()
    else:
        mc.model_setup()
        mc.boundaries_setup()
        mc.initialize_logchi2()

    mc.starting_points_setup()

    results_analysis.print_bayesian_info(mc)

    discard_exploration = mc.nested_sampling_parameters['discard_exploration']
    nautilus_nlive_multiplier = mc.nested_sampling_parameters['nautilus_nlive_multiplier']

    if 'nlive_mult' in mc.nested_sampling_parameters:
        nlive = int(mc.ndim * mc.nested_sampling_parameters['nlive_mult'] * nautilus_nlive_multiplier)
    else:
        nlive = int(mc.nested_sampling_parameters['nlive'] * nautilus_nlive_multiplier)

    n_networks = mc.nested_sampling_parameters['n_networks']
    equal_weight_boost = mc.nested_sampling_parameters['equal_weight_boost']


    global nautilus_priors
    def nautilus_priors(cube):
        return mc.nautilus_priors(cube)

    global nautilus_loglikelihood
    def nautilus_loglikelihood(theta):
        return mc.nautilus_call(theta)



    print('Number of live points:', nlive)
    print()
    print('number of multiprocessing threads:', num_threads)

    print()
    print('Reference Time Tref: ', mc.Tref)
    print()
    print('*************************************************************')
    print()

    try:
        #import nautilus
        from nautilus import Sampler
    except (ModuleNotFoundError,ImportError):
        print("ERROR: nautilus not installed, this will not work")
        quit()


    #dlogz = mc.nested_sampling_parameters.get('dlogz', 0.01)
    use_threading_pool = mc.nested_sampling_parameters.get('use_threading_pool', True)
    #pfrac_value = mc.nested_sampling_parameters.get('pfrac', 0.000)
    use_default = (mc.nested_sampling_parameters.get('default', False) | mc.nested_sampling_parameters.get('use_default', False))

    print('Using threading pool for nautilus:', use_threading_pool)

    if not reloaded_nautilus:
        """ A dummy file is created to let the cpulimit script to proceed with the next step"""
        nested_sampling_write_dummy_file(mc)
        nested_sampling_save_to_cpickle(mc)

    if use_default:
        print('Setting up nautilus, number of live points = {0:6.0f}'.format(nlive))
        print()
    else:
        print('Setting up nautilus, number of live points = {0:6.0f}'.format(nlive))
        #print('                                        posterior/evidence split = {0:4.3f}'.format(pfrac_value))
        #print('                                        initial stopping criterion = {0:5.4f}'.format(dlogz))
        #print('                                        bound = ', mc.nested_sampling_parameters['bound'])
        #print('                                        sample = ', mc.nested_sampling_parameters['sample'])

    print('                     live points multiplier = {0:3.0f}'.format(nautilus_nlive_multiplier))
    print('                     number of networks = {0:6.0f}'.format(n_networks))
    print('                     equal weights boost = {0:6.0f}'.format(equal_weight_boost))
    print()

    if not use_threading_pool:
        num_threads = None

    sampler = Sampler(
        prior=nautilus_priors,
        likelihood=nautilus_loglikelihood,
        n_dim=mc.ndim,
        n_live=nlive,
        pool=num_threads,
        filepath = save_checkpoint,
        n_networks = n_networks
    )
    success = sampler.run(discard_exploration=discard_exploration, verbose=True)
    sampler.write(save_sampler, overwrite=True)


    print()


    unweighted_points, unweighted_log_w, unweighted_log_l = sampler.posterior()
    unweighted_log_z = sampler.log_z
    points, log_w, log_l = sampler.posterior(return_as_dict=True,
                                                equal_weight=True,
                                                equal_weight_boost=equal_weight_boost)
    effective_sample_size = sampler.effective_sample_size()
    log_z = sampler.log_z

    results = {
        'unweighted_points': unweighted_points,
        'unweighted_log_w': unweighted_log_w,
        'unweighted_log_l': unweighted_log_l,
        'unweighted_log_z': unweighted_log_z,
        'points': points,
        'log_w': log_w,
        'log_l': log_l,
        'log_z': log_z,
        'effective_sample_size': effective_sample_size,
        'success': success,
    }

    print('Getting the results')

    #taken from dynesty/dynesty/results.py  but without the nlive point causing an error
    res = ("success: {:}\n"
           "effective_sample_size: {:f}\n"
            "logz: {:6.3f}"
            .format(success, effective_sample_size, log_z))
    print()
    print('Summary\n=======\n'+res)

    print()
    print()


    """ A dummy file is created to let the cpulimit script to proceed with the next step"""
    #nested_sampling_write_dummy_file(mc)
    #nested_sampling_save_to_cpickle(mc)
    #nautilus_sampler_save_to_cpickle(mc.output_directory, sampler)
    nautilus_results_save_to_cpickle(mc.output_directory, results)

    #try:
    #    os.remove(save_checkpoint)
    #except:
    #    pass

    if return_output:
        return mc
    else:
        return
