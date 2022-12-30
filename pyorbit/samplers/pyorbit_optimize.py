from __future__ import print_function
from pyorbit.subroutines.common import *
from pyorbit.classes.model_container_optimize import ModelContainerOptimize
from pyorbit.subroutines.input_parser import yaml_parser, pars_input
from pyorbit.subroutines.io_subroutines import starting_point_load_from_cpickle, starting_point_save_to_cpickle
import pyorbit.subroutines.results_analysis as results_analysis
from scipy import optimize
import os


__all__ = ["pyorbit_optimize"]


def pyorbit_optimize(config_in, input_datasets=None, return_output=None):

    optimize_dir_output = './' + config_in['output'] + '/optimize/'

    print()

    mc = ModelContainerOptimize()

    pars_input(config_in, mc, input_datasets)

    mc.model_setup()
    mc.boundaries_setup()
    mc.initialize_logchi2()

    mc.starting_points_setup()

    mc.optimize_dir_output = optimize_dir_output

    if not os.path.exists(mc.optimize_dir_output):
        os.makedirs(mc.optimize_dir_output)

    print()
    print('Reference Time Tref: ', mc.Tref)



    # method keyword is removed to suppress optimize warning
    mc.optimize_options = mc.optimize_parameters.copy()
    mc.optimize_options.pop('method')

    print()
    print('Scipy.optimixe.minize method: ', mc.optimize_parameters['method'])
    print('Options: ')
    for key_name, key_val in mc.optimize_options.items():
        print('   {0:15s}:  {1}'.format(key_name, key_val))

    print()
    results_analysis.results_summary(mc, mc.starting_point, compute_lnprob=False, is_starting_point=True)


    optimize_results = optimize.minimize(
        mc.negative_log_priors_likelihood,  # we minimize the negative loglikelyhood (including priors)
        mc.starting_point,  # initial parameters parameters
        method=mc.optimize_parameters['method'],
        options=mc.optimize_options
    )
    output_results = optimize_results['x']
    print(optimize_results['success'])

    if optimize_results['success']:

        results_analysis.results_summary(mc, output_results, compute_lnprob=True)

        theta_dict = results_analysis.get_theta_dictionary(mc)
        starting_point_save_to_cpickle(optimize_dir_output, output_results, mc.bounds, theta_dict)

    else:
        print('Unsuccessful optimization - results have not been saved, please try again by changing method or options')
        print()
        print('Accepted options by scipy.optimize.minimize for method :', mc.optimize_parameters['method'])
        print('You can set them in the optimize section of the YAML input file.\n')
        optimize.show_options(solver='minimize', method=mc.optimize_parameters['method'])
        print()
