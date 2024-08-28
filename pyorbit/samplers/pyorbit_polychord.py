from __future__ import print_function
from pyorbit.subroutines.common import *
from pyorbit.classes.model_container_polychord import ModelContainerPolyChord
from pyorbit.subroutines.input_parser import yaml_parser, pars_input
from pyorbit.subroutines.io_subroutines import nested_sampling_save_to_cpickle, nested_sampling_load_from_cpickle, \
    nested_sampling_write_dummy_file
import pyorbit.subroutines.results_analysis as results_analysis
import os
import sys
import argparse

__all__ = ["pyorbit_polychord"]

"""
def show(filepath):
    # open the output (pdf) file for the user
    if os.name == 'mac': subprocess.call(('open', filepath))
    elif os.name == 'nt': os.startfile(filepath)
"""


# | Optional dumper function giving run-time read access to
# | the live points, dead points, weights and evidences
def dumper(live, dead, logweights, logZ, logZerr):
    print("Last dead point:", dead[-1])


def pyorbit_polychord(config_in, input_datasets=None, return_output=None):
    output_directory = './' + config_in['output'] + '/polychord/'

    mc = ModelContainerPolyChord()
    pars_input(config_in, mc, input_datasets)

    if mc.nested_sampling_parameters['shutdown_jitter']:
        for dataset_name, dataset in mc.dataset_dict.items():
            dataset.shutdown_jitter()

    mc.model_setup()
    mc.boundaries_setup()
    mc.initialize_logchi2()

    mc.starting_points_setup()

    results_analysis.print_bayesian_info(mc)

    mc.output_directory = output_directory

    # os.system("mkdir -p " + output_directory + "/clusters")
    # os.system("mkdir -p " +output_directory + "chains/clusters")

    print()
    print('Reference Time Tref: ', mc.Tref)
    print()
    print('*************************************************************')
    print()

    try:
        import pypolychord
        from pypolychord.settings import PolyChordSettings
    except (ModuleNotFoundError,ImportError):
        try:
            import PyPolyChord
            from PyPolyChord.settings import PolyChordSettings
            print('Consider updating to newer version of PolyChord: https://github.com/PolyChord/PolyChordLite ')
        except (ModuleNotFoundError,ImportError):
            print('ERROR: pypolychord not installed, this will not work')
            quit()

    settings = PolyChordSettings(nDims=mc.ndim, nDerived=0)

    settings.file_root = 'pyorbit'
    settings.base_dir = output_directory

    for key_name, key_value in mc.nested_sampling_parameters.items():

        if hasattr(settings, key_name):
            setattr(settings, key_name, key_value)

    if 'nlive_mult' in mc.nested_sampling_parameters:
        setattr(settings, 'nlive', mc.ndim * mc.nested_sampling_parameters['nlive_mult'])

    if 'num_repeats_mult' in mc.nested_sampling_parameters:
        setattr(settings, 'num_repeats', mc.ndim * mc.nested_sampling_parameters['num_repeats_mult'])

    if 'include_priors' in mc.nested_sampling_parameters:
        mc.include_priors = mc.nested_sampling_parameters['include_priors']

    output = pypolychord.run_polychord(mc.polychord_call, nDims=mc.ndim, nDerived=0, settings=settings,
                                       prior=mc.polychord_priors, dumper=dumper)

    paramnames = [('p%i' % i, r'\theta_%i' % i) for i in range(mc.ndim)]
    paramnames += [('r*', 'r')]
    output.make_paramnames_files(paramnames)

    nested_sampling_save_to_cpickle(mc)

    print()
    print('PolyChord COMPLETED')
    print()

    """ A dummy file is written to let the cpulimit script to proceed with the next step"""
    nested_sampling_write_dummy_file(mc)

    if return_output:
        return mc
    else:
        return
