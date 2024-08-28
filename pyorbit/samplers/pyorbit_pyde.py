from __future__ import print_function
from pyorbit.subroutines.common import np
from pyorbit.classes.model_container_emcee import ModelContainerEmcee
from pyorbit.subroutines.input_parser import yaml_parser, pars_input
from pyorbit.subroutines.io_subroutines import pyde_save_to_pickle,\
    pyde_load_from_cpickle
import pyorbit.subroutines.results_analysis as results_analysis
import os
import sys
import multiprocessing

__all__ = ["pyorbit_pyde"]


def pyorbit_pyde(config_in, input_datasets=None, return_output=None):

    ## Check how many CPU threads (I guess) should be used
    # omp_num_threads = config_in['parameters'].get('cpu_threads', "1")
    # if type(omp_num_threads) == type("1"):
    #     os.environ["OMP_NUM_THREADS"] = omp_num_threads
    # else:
    #     os.environ["OMP_NUM_THREADS"] = "{0:.0f}".format(omp_num_threads)
    #
    # pyde_dir_output = './' + config_in['output'] + '/pyde/'
    #
    # reloaded_pyde = False


    try:
        mc, population, starting_point, theta_dict = pyde_load_from_cpickle(
            pyde_dir_output, prefix='')
        reloaded_pyde = True
    except FileNotFoundError:
        pass


    os.environ["OMP_NUM_THREADS"] = "1"
    try:
        num_threads = int(config_in['parameters'].get('cpu_threads', "1"))
    except:
        num_threads = multiprocessing.cpu_count()-1

    multiprocessing.set_start_method('fork')

    print()
    print('reloaded_pyde: ', reloaded_pyde)
    print()
    print('number of multiprocessing threads:', num_threads )
    print("multiprocessing method (should be fork): ", multiprocessing.get_start_method())
    print()

    if reloaded_pyde:
        print("ERROR: PyDE output already present, nothing will be done")

    mc = ModelContainerEmcee()
    pars_input(config_in, mc, input_datasets)

    if mc.pyde_parameters['shutdown_jitter'] or mc.emcee_parameters['shutdown_jitter']:
        for dataset_name, dataset in mc.dataset_dict.items():
            dataset.shutdown_jitter()

    # keep track of which version has been used to perform emcee computations

    mc.model_setup()
    mc.boundaries_setup()
    mc.initialize_logchi2()

    results_analysis.print_bayesian_info(mc)

    mc.pyde_dir_output = pyde_dir_output

    print('Include priors: ', mc.include_priors)
    print()
    print('Reference Time Tref: ', mc.Tref)
    print()
    print('Dimensions = ', mc.ndim)
    print()


    """ We have to run PyDE """
    try:
        from pyde.de import DiffEvol
    except (ModuleNotFoundError,ImportError):
        print('ERROR! PyDE is not installed, run first with optimize instead of emcee')
        quit()

    if not os.path.exists(mc.pyde_dir_output):
        os.makedirs(mc.pyde_dir_output)

    print('Using threading pool for PyDE:',
        mc.pyde_parameters['use_threading_pool'])

    print('PyDE running')
    sys.stdout.flush()

    if mc.pyde_parameters['use_threading_pool']:
        with multiprocessing.Pool(num_threads) as pool:

            de = DiffEvol(
                mc,
                mc.bounds,
                mc.emcee_parameters['nwalkers'],
                maximize=True,
                pool=pool)

            de.optimize(int(mc.pyde_parameters['ngen']))
    else:
        de = DiffEvol(
            mc,
            mc.bounds,
            mc.emcee_parameters['nwalkers'],
            maximize=True)

        de.optimize(int(mc.pyde_parameters['ngen']))

    population = de.population
    starting_point = np.median(population, axis=0)

    theta_dict = results_analysis.get_theta_dictionary(mc)

    """ bounds redefinition and fix for PyDE anomalous results """
    if mc.recenter_bounds_flag:
        pyde_save_to_pickle(
            mc, population, starting_point, theta_dict, prefix='orig')

        mc.recenter_bounds(starting_point)
        population = mc.fix_population(starting_point, population)
        starting_point = np.median(population, axis=0)

        print('Boundaries redefined after PyDE output')

    pyde_save_to_pickle(mc, population, starting_point, theta_dict)

    print('PyDE completed')
    print()
    sys.stdout.flush()

    results_analysis.results_summary(
        mc, starting_point, compute_lnprob=True, is_starting_point=True)
    sys.stdout.flush()

    print()
    print('*************************************************************')
