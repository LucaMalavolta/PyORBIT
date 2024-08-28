from __future__ import print_function
from pyorbit.subroutines.common import np
from pyorbit.classes.model_container_zeus import ModelContainerZeus
from pyorbit.subroutines.input_parser import yaml_parser, pars_input
from pyorbit.subroutines.io_subroutines import pyde_save_to_pickle,\
    pyde_load_from_cpickle,\
    zeus_save_to_cpickle, zeus_load_from_cpickle, zeus_flatchain,\
    zeus_write_dummy_file, starting_point_load_from_cpickle, zeus_simpler_load_from_cpickle
import pyorbit.subroutines.results_analysis as results_analysis
import os
import sys
import multiprocessing

__all__ = ["pyorbit_zeus_legacy"]


def pyorbit_zeus_legacy(config_in, input_datasets=None, return_output=None):

    try:
        import zeus
    except:
        print("ERROR: zeus not installed, this will not work")
        quit()

    os.environ["OMP_NUM_THREADS"] = "1"
    try:
        num_threads = int(config_in['parameters'].get('cpu_threads', "1"))
    except:
        num_threads = multiprocessing.cpu_count()-1

    optimize_dir_output = './' + config_in['output'] + '/optimize/'
    pyde_dir_output = './' + config_in['output'] + '/pyde/'

    zeus_dir_output = './' + config_in['output'] + '/zeus/'

    reloaded_optimize = False
    reloaded_pyde = False
    reloaded_zeus = False

    try:
        mc, population, starting_point, theta_dict = pyde_load_from_cpickle(
            pyde_dir_output, prefix='')
        reloaded_pyde = True
    except FileNotFoundError:
        pass

    try:
        mc, starting_point, population, prob, sampler_chain, \
            sampler_lnprobability, _, theta_dict = \
            zeus_load_from_cpickle(zeus_dir_output)
        state, sampler = zeus_simpler_load_from_cpickle(zeus_dir_output)
        reloaded_zeus = True
    except FileNotFoundError:
        pass

    try:
        starting_point, previous_boundaries, theta_dict = starting_point_load_from_cpickle(
            optimize_dir_output)
        reloaded_optimize = True
    except FileNotFoundError:
        pass

    print()
    print('reloaded_optimize: ', reloaded_optimize)
    print('reloaded_pyde: ', reloaded_pyde)
    print('reloaded_zeus: ', reloaded_zeus)
    print()

    if reloaded_pyde or reloaded_zeus:
        previous_boundaries = mc.bounds

    if reloaded_zeus:
        print('Requested steps:', mc.zeus_parameters['nsteps'])

        mc.zeus_parameters['completed_nsteps'] = \
            int(sampler_chain.shape[1] * mc.zeus_parameters['thin'])

        print('Completed:', mc.zeus_parameters['completed_nsteps'] )
        pars_input(config_in, mc, input_datasets, reload_zeus=True)
        print('Total:', mc.zeus_parameters['nsteps'])

        """ There's no need to do anything"""
        flatchain = zeus_flatchain(
            sampler_chain, mc.zeus_parameters['nburn'], mc.zeus_parameters['thin'])
        mc.model_setup()
        mc.initialize_logchi2()
        results_analysis.print_bayesian_info(mc)

        results_analysis.print_integrated_ACF(
            sampler_chain, theta_dict, mc.zeus_parameters['thin'])

        """ In case the current startin point comes from a previous analysis """
        mc.zeus_dir_output = zeus_dir_output

        if mc.zeus_parameters['nsteps'] <= mc.zeus_parameters['completed_nsteps']:

            print('Reference Time Tref: ', mc.Tref)
            print()
            print('Dimensions = ', mc.ndim)
            print('Nwalkers = ', mc.zeus_parameters['nwalkers'])
            print('Steps = ', mc.zeus_parameters['nsteps'])
            print()
            print('Original starting point of zeus:')
            print()

            results_analysis.results_summary(
                mc, starting_point, compute_lnprob=True, is_starting_point=True)

            print('zeus results:')
            print()

            results_analysis.results_summary(mc, flatchain)

            print()
            print('zeus completed')
            print()

            if return_output:
                return mc, sampler_chain, sampler_lnprobability
            else:
                return

        elif not sampler:
            print('Sampler file is missing - only analysis performed with PyORBIT >8.1 cn be resumed')
            if return_output:
                return mc, sampler_chain, sampler_lnprobability
            else:
                return

    if reloaded_zeus:
        sampled = mc.zeus_parameters['completed_nsteps']
        nsteps_todo = mc.zeus_parameters['nsteps'] \
            - mc.zeus_parameters['completed_nsteps']

        print('Resuming from a previous run:')
        print('Performed steps = ', mc.zeus_parameters['completed_nsteps'])
        print('Final # of steps = ', mc.zeus_parameters['nsteps'])
        print('Steps to be performed = ', nsteps_todo)
        print()

    else:

        mc = ModelContainerZeus()
        pars_input(config_in, mc, input_datasets)

        if mc.pyde_parameters['shutdown_jitter'] or mc.zeus_parameters['shutdown_jitter']:
            for dataset_name, dataset in mc.dataset_dict.items():
                dataset.shutdown_jitter()

        # keep track of which version has been used to perform zeus computations
        mc.zeus_parameters['version'] = zeus.__version__[0]

        mc.model_setup()
        mc.boundaries_setup()
        mc.initialize_logchi2()

        results_analysis.print_bayesian_info(mc)

        mc.pyde_dir_output = pyde_dir_output
        mc.zeus_dir_output = zeus_dir_output

        mc.zeus_parameters['nwalkers'] = mc.ndim * \
            mc.zeus_parameters['npop_mult']
        if mc.zeus_parameters['nwalkers'] % 2 == 1:
            mc.zeus_parameters['nwalkers'] += 1

        if not os.path.exists(mc.zeus_dir_output):
            os.makedirs(mc.zeus_dir_output)

        state = None
        sampled = 0
        nsteps_todo = mc.zeus_parameters['nsteps']

    print('Include priors: ', mc.include_priors)
    print()
    print('Reference Time Tref: ', mc.Tref)
    print()
    print('Dimensions = ', mc.ndim)
    print('Nwalkers = ', mc.zeus_parameters['nwalkers'])
    print()

    if mc.zeus_parameters['version'] == '2':
        mc.zeus_parameters['use_threading_pool'] = False

    if not mc.pyde_parameters.get('use_threading_pool', False):
        mc.pyde_parameters['use_threading_pool'] = False

    if not mc.zeus_parameters.get('use_threading_pool', False):
        mc.zeus_parameters['use_threading_pool'] = False


    if reloaded_zeus:
        sys.stdout.flush()
        pass
    elif reloaded_pyde:

        theta_dict_legacy = theta_dict.copy()
        population_legacy = population.copy()

        theta_dict = results_analysis.get_theta_dictionary(mc)
        population = np.zeros(
            [mc.zeus_parameters['nwalkers'], mc.ndim], dtype=np.double)

        for theta_name, theta_i in theta_dict.items():
            population[:, theta_i] = population_legacy[:,
                                                       theta_dict_legacy[theta_name]]
            mc.bounds[theta_i] = previous_boundaries[theta_dict_legacy[theta_name]]

        starting_point = np.median(population, axis=0)
        # print(starting_point)
        # print(population)

        print('Using previous population as starting point. ')
        print()
        sys.stdout.flush()

    elif mc.starting_point_flag or reloaded_optimize:

        if reloaded_optimize:
            print('Using the output from a previous optimize run as starting point')
            theta_dict_legacy = theta_dict.copy()
            starting_point_legacy = starting_point.copy()
            theta_dict = results_analysis.get_theta_dictionary(mc)
            for theta_name, theta_i in theta_dict.items():
                starting_point[theta_i] = starting_point_legacy[theta_dict_legacy[theta_name]]
        else:
            print('Using user-defined starting point from YAML file')
            mc.starting_points_setup()
            starting_point = mc.starting_point

        population = np.zeros(
            [mc.zeus_parameters['nwalkers'], mc.ndim], dtype=np.double)
        for ii in range(0, mc.zeus_parameters['nwalkers']):
            population[ii, :] = np.random.normal(starting_point, 0.0000001)

        print(
            'to write a synthetic population extremely close to the starting values.')
        print()

        sys.stdout.flush()

    else:
        """ None of the previous cases has been satisfied, we have to run PyDE """
        try:
            from pyde.de import DiffEvol
        except (ModuleNotFoundError,ImportError):
            print('ERROR! PyDE is not installed, run first with optimize instead of zeus')
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
                    mc.zeus_parameters['nwalkers'],
                    maximize=True,
                    pool=pool)

                de.optimize(int(mc.pyde_parameters['ngen']))
        else:
            de = DiffEvol(
                mc,
                mc.bounds,
                mc.zeus_parameters['nwalkers'],
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

    if reloaded_zeus:
        print('Original starting point of zeus:')
        print()

    results_analysis.results_summary(
        mc, starting_point, compute_lnprob=True, is_starting_point=True)
    sys.stdout.flush()

    print()
    print('*************************************************************')
    print()
    print('Running zeus')
    print()
    print('zeus version: ', zeus.__version__)
    print('Using threading pool for zeus:', mc.zeus_parameters.get('use_threading_pool', False))
    print()

    if reloaded_zeus:
        print('Using reloaded sampler')
        print()
    else:
        sampler = zeus.EnsembleSampler(
            mc.zeus_parameters['nwalkers'], mc.ndim, mc)

    if mc.zeus_parameters['nsave'] > 0:
        print('Saving temporary steps not supported')
        print()

    if mc.zeus_parameters['use_threading_pool']:
        with multiprocessing.Pool(num_threads) as pool:
            sampler.pool = pool
            sampler.run_mcmc(
                population,
                nsteps_todo,
                thin=mc.zeus_parameters['thin'],
                progress=True)

    else:
        sampler.run_mcmc(
            population,
            nsteps_todo,
            thin=mc.zeus_parameters['thin'],
            progress=True)

    population = sampler.get_chain()
    prob= sampler.get_log_prob()
    state = None

    sampled += nsteps_todo

    theta_dict = results_analysis.get_theta_dictionary(mc)
    zeus_save_to_cpickle(mc, starting_point, population,
                            prob, state, sampler, theta_dict,
                            samples=sampled)

    flatchain = zeus_flatchain(
        sampler.chain,
        mc.zeus_parameters['nburn'],
        mc.zeus_parameters['thin'])

    results_analysis.print_integrated_ACF(
        sampler.chain,
        theta_dict,
        mc.zeus_parameters['thin'])

    results_analysis.results_summary(mc, flatchain)

    print()
    print('zeus completed')
    print()

    """ A dummy file is written to let the cpulimit script to proceed with the next step"""
    zeus_write_dummy_file(mc)

    if return_output:
        return mc, sampler.chain,  sampler.lnprobability
