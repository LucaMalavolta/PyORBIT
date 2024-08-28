from __future__ import print_function
from pyorbit.classes.model_container_ultranest import ModelContainerUltranest
from pyorbit.subroutines.input_parser import yaml_parser, pars_input
from pyorbit.subroutines.io_subroutines import nested_sampling_save_to_cpickle, \
    nested_sampling_load_from_cpickle, nested_sampling_write_dummy_file, \
    ultranest_sampler_save_to_cpickle, ultranest_sampler_load_from_cpickle

import pyorbit.subroutines.results_analysis as results_analysis
import os
import sys
import re
import argparse
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt

__all__ = ["pyorbit_ultranest"]

"""
def show(filepath):
    # open the output (pdf) file for the user
    if os.name == 'mac': subprocess.call(('open', filepath))
    elif os.name == 'nt': os.startfile(filepath)
"""


def pyorbit_ultranest(config_in, input_datasets=None, return_output=None):

    mc = ModelContainerUltranest()
    pars_input(config_in, mc, input_datasets)

    mc.output_directory = './' + config_in['output'] + '/ultranest/'

    os.environ["OMP_NUM_THREADS"] = "1"

    try:
        results = ultranest_sampler_load_from_cpickle(mc.output_directory)
        print('UltraNest results already saved in the respective directory, run PyORBIT_GetResults')
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

    results_analysis.print_bayesian_info(mc)

    theta_dictionary = results_analysis.get_theta_dictionary(mc)
    labels_array = [None] * len(theta_dictionary)
    for key_name, key_value in theta_dictionary.items():
        labels_array[key_value] = re.sub('_', '-', key_name)

    if 'nlive_mult' in mc.nested_sampling_parameters:
        nlive = mc.ndim * mc.nested_sampling_parameters['nlive_mult']
    else:
        nlive = mc.nested_sampling_parameters['nlive']

    print('Number of minimum live points:', nlive)
    print('Desired accuracy:', mc.nested_sampling_parameters['desired_accuracy'])
    print('Minimum number of effective samples:', mc.nested_sampling_parameters['min_ess'])

    print()
    print('Reference Time Tref: ', mc.Tref)
    print()
    print('*************************************************************')
    print()

    try:
        from ultranest import ReactiveNestedSampler
    except (ModuleNotFoundError,ImportError):
        print("ERROR: ultranest not installed, this will not work")
        quit()

    global ultranest_transform
    def ultranest_transform(cube):
        return mc.ultranest_transform(cube)

    global ultranest_call
    def ultranest_call(theta):
        return mc.ultranest_call(theta)

    sampler = ReactiveNestedSampler(
        labels_array,
        ultranest_call,
        transform=ultranest_transform,
        log_dir=mc.output_directory, # folder where to store files
        resume=True, # whether to resume from there (otherwise start from scratch)
    )

    sampler.run(
        min_num_live_points=nlive,
        dlogz=mc.nested_sampling_parameters['desired_accuracy'], # desired accuracy on logz
        min_ess=mc.nested_sampling_parameters['min_ess'], # number of effective samples
        max_num_improvement_loops=mc.nested_sampling_parameters['improvement_loops'] # how many times to go back and improve
    )

    sampler.print_results()

    print()
    sampler.plot()
    sampler.plot_trace()
    sampler.plot_run()

    """ A dummy file is written to let the cpulimit script to proceed with the next step"""
    nested_sampling_write_dummy_file(mc)
    nested_sampling_save_to_cpickle(mc)
    # ultranest_sampler_save_to_cpickle(mc.output_directory, sampler)

    if return_output:
        return mc
    else:
        return
