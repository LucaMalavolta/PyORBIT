from __future__ import print_function
from pyorbit.classes.model_container_ultranest_warmstart import ModelContainerUltranestWarmstart
from pyorbit.subroutines.input_parser import yaml_parser, pars_input
from pyorbit.subroutines.io_subroutines import nested_sampling_save_to_cpickle, \
    nested_sampling_load_from_cpickle, nested_sampling_write_dummy_file, \
    ultranest_sampler_save_to_cpickle, ultranest_sampler_load_from_cpickle
from pyorbit.subroutines.io_subroutines import *

import pyorbit.subroutines.results_analysis as results_analysis

import pyorbit.subroutines.results_analysis as results_analysis
import os
import sys
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from time import time

from pyorbit.subroutines.common import np, nested_sampling_prior_compute

__all__ = ["pyorbit_ultranest_warmstart"]

"""
def show(filepath):
    # open the output (pdf) file for the user
    if os.name == 'mac': subprocess.call(('open', filepath))
    elif os.name == 'nt': os.startfile(filepath)
"""


def pyorbit_ultranest_warmstart(config_in, input_datasets=None, return_output=None):


    mc = ModelContainerUltranestWarmstart()
    pars_input(config_in, mc, input_datasets)

    mc.output_directory = './' + config_in['output'] + '/ultranest_warmstart/'
    mc.pyde_dir_output = './' + config_in['output'] + '/pyde_warmup/'
    mc.emcee_dir_output = './' + config_in['output'] + '/emcee_warmup/'

    try:
        results = ultranest_sampler_load_from_cpickle(mc.output_directory)
        print('UltraNest results already saved in the respective directory, run PyORBIT_GetResults')
        if return_output:
            return mc
        else:
            return
    except FileNotFoundError:
        pass

    os.environ["OMP_NUM_THREADS"] = "1"
    try:
        num_threads = int(config_in['parameters'].get('cpu_threads',  multiprocessing.cpu_count()))
    except:
        print(" Something happened when trying to setup multiprocessing, switching back to 1 CPU")
        num_threads = 1


    reloaded_pyde = False
    reloaded_emcee = False

    try:
        mc, population, starting_point, theta_dict = pyde_load_from_cpickle(
            mc.pyde_dir_output, prefix='')
        reloaded_pyde = True
    except FileNotFoundError:
        pass

    try:
        mc, starting_point, population, prob, sampler_chain, \
            sampler_lnprobability, _, theta_dict = \
            emcee_load_from_cpickle(mc.emcee_dir_output)
        reloaded_emcee = True
    except FileNotFoundError:
        pass

    print('reloaded_pyde : ', reloaded_pyde)
    print('reloaded_emcee: ', reloaded_emcee)
    print()
    print('number of multiprocessing threads:', num_threads)
    print()

    safe_reload = config_in['parameters'].get('safe_reload', False)

    print('parameters/safe_reload flag (must be True for tinygp): ', safe_reload)
    print()


    if mc.nested_sampling_parameters['shutdown_jitter']:
        'Jitter term not included for evidence calculation'
        print()
        for dataset_name, dataset in mc.dataset_dict.items():
            dataset.shutdown_jitter()

    mc.model_setup()
    if not reloaded_pyde:
        mc.boundaries_setup()
    mc.initialize_logchi2()

    mc.starting_points_setup()

    results_analysis.print_bayesian_info(mc)
    theta_dictionary = results_analysis.get_theta_dictionary(mc)

    mc.emcee_parameters = mc.emcee_warmup_parameters
    mc.pyde_parameters = mc.pyde_warmup_parameters

    try:
        import ultranest
        import ultranest.stepsampler
    except (ModuleNotFoundError,ImportError):
        print("ERROR: ultranest not installed, this will not work")
        quit()

    try:
        import emcee
    except:
        print("ERROR: emcee not installed, this will not work")
        quit()

    try:
        from pyde.de import DiffEvol
    except (ModuleNotFoundError,ImportError):
        print('ERROR! PyDE is not installed, run first with optimize instead of emcee')
        quit()



    global ultranest_transform
    def ultranest_transform(cube):
        theta = np.zeros(len(cube), dtype=np.double)

        for i in range(0, len(cube)):
            theta[i] = nested_sampling_prior_compute(
                cube[i], mc.priors[i][0], mc.priors[i][2], mc.spaces[i])
        return theta

    global ultranest_call
    def ultranest_call(theta):
        _, log_likelihood = logprior_loglikelihood(theta)

        if not mc.check_bounds(theta):
            return -0.5e8

        if not np.isfinite(log_likelihood):
            return -0.5e8

        return log_likelihood

    global pyde_emcee_call
    def pyde_emcee_call(cube):
        theta = ultranest_transform(cube)

        """
        Constant term added either by dataset.model_logchi2() or gp.log_likelihood()
        """
        if not mc.check_bounds(theta):
            return -np.inf

        log_priors, log_likelihood = logprior_loglikelihood(theta, keep_finite=False)

        return log_priors + log_likelihood

    global logprior_loglikelihood
    def logprior_loglikelihood(theta, keep_finite=True):

        log_priors = 0.00
        log_likelihood = 0.00

        if mc.dynamical_model is not None:
            """ check if any keyword ahas get the output model from the dynamical tool
            we must do it here because all the planet are involved"""
            dynamical_output = mc.dynamical_model.compute(mc, theta)

        for model_name, model in mc.common_models.items():
            log_priors += model.return_priors(theta)

        delayed_lnlk_computation = []
        residuals_analysis = {}

        for dataset_name, dataset in mc.dataset_dict.items():

            logchi2_gp_model = None
            compute_gp_residuals = False

            dataset.model_reset()
            parameter_values = dataset.convert(theta)
            dataset.compute(parameter_values)

            log_priors += dataset.return_priors(theta)

            if 'none' in dataset.models or 'None' in dataset.models:
                continue
            if not dataset.models:
                continue

            skip_loglikelihood = False

            for model_name in dataset.models:

                log_priors += mc.models[model_name].return_priors(
                    theta, dataset_name)

                parameter_values = {}
                for common_ref in mc.models[model_name].common_ref:
                    parameter_values.update(
                        mc.common_models[common_ref].convert(theta))

                #TODO: remove try-except starting from version 11 !!
                try:
                    for planet_name in mc.models[model_name].multiple_planets:
                        parameter_values.update(
                            mc.common_models[planet_name].convert_with_name(theta, planet_name))
                except:
                    pass

                parameter_values.update(
                    mc.models[model_name].convert(theta, dataset_name))

                """ residuals will be computed following the definition in Dataset class
                """
                if getattr(mc.models[model_name], 'residuals_analysis', False):

                    skip_loglikelihood = True
                    if mc.models[model_name].gp_before_correlation:
                        compute_gp_residuals = True


                    try:
                        residuals_analysis[model_name]['parameter_values'] = parameter_values
                    except:
                        residuals_analysis[model_name] = {'parameter_values': parameter_values}

                    if dataset_name == mc.models[model_name].x_dataset:
                        residuals_dataset_label = 'x'
                        residuals_analysis[model_name]['x_gp_model'] = None
                        residuals_analysis[model_name]['x_gp_parameters'] = None
                    else:
                        residuals_dataset_label = 'y'
                        residuals_analysis[model_name]['y_gp_model'] = None
                        residuals_analysis[model_name]['y_gp_parameters'] = None
                    continue


                if getattr(mc.models[model_name], 'internal_likelihood', False):
                    logchi2_gp_model = model_name
                    continue

                # if getattr(mc.models[model_name], 'model_class', None) is 'common_jitter':
                if getattr(mc.models[model_name], 'jitter_model', False):
                    dataset.jitter += mc.models[model_name].compute(
                        parameter_values, dataset)
                    continue

                if getattr(dataset, 'dynamical', False):
                    dataset.external_model = dynamical_output[dataset_name]

                if dataset.normalization_model is None and (mc.models[model_name].unitary_model or mc.models[model_name].normalization_model):
                    dataset.normalization_model = np.ones(dataset.n, dtype=np.double)

                if mc.models[model_name].unitary_model:
                    dataset.unitary_model += mc.models[model_name].compute(
                    parameter_values, dataset)
                elif mc.models[model_name].normalization_model:
                    dataset.normalization_model *= mc.models[model_name].compute(
                        parameter_values, dataset)
                else:
                    dataset.additive_model += mc.models[model_name].compute(
                        parameter_values, dataset)

            dataset.compute_model()
            dataset.compute_residuals()

            """ Gaussian Process check MUST be the last one or the program will fail
             that's because for the GP to work we need to know the _deterministic_ part of the model
             (i.e. the theoretical values you get when you feed your model with the parameter values) """

            if logchi2_gp_model:

                parameter_values = {}
                for common_ref in mc.models[logchi2_gp_model].common_ref:
                    parameter_values.update(
                        mc.common_models[common_ref].convert(theta))

                parameter_values.update(
                    mc.models[logchi2_gp_model].convert(theta, dataset_name))

                """ GP Log-likelihood is not computed now because a single matrix must be
                    computed with the joint dataset"""
                if hasattr(mc.models[logchi2_gp_model], 'delayed_lnlk_computation'):

                    mc.models[logchi2_gp_model].add_internal_dataset(parameter_values, dataset)
                    if logchi2_gp_model not in delayed_lnlk_computation:
                        delayed_lnlk_computation.append(logchi2_gp_model)
                elif skip_loglikelihood:
                    residuals_analysis[model_name][residuals_dataset_label+'_gp_model'] = logchi2_gp_model
                    residuals_analysis[model_name][residuals_dataset_label+'_gp_parameters'] = parameter_values
                else:
                    log_likelihood += mc.models[logchi2_gp_model].lnlk_compute(
                        parameter_values, dataset)

                if compute_gp_residuals:
                    dataset.residuals_for_regression -= mc.models[logchi2_gp_model].sample_predict(parameter_values, dataset)
                    residuals_analysis[model_name][residuals_dataset_label+'_gp_model'] = None

            elif not skip_loglikelihood:
                log_likelihood += dataset.model_logchi2()

        """ Correlation between residuals of two datasets through orthogonal distance regression
            it must be computed after any other model has been removed from the
            independent variable, if not provided as ancillary dataset
        """
        for model_name in residuals_analysis:
            parameter_values =  residuals_analysis[model_name]['parameter_values']
            x_dataset = mc.dataset_dict[mc.models[model_name].x_dataset]
            y_dataset = mc.dataset_dict[mc.models[model_name].y_dataset]
            modelout_xx, modelout_yy = mc.models[model_name].compute(parameter_values, x_dataset, y_dataset)

            x_dataset.residuals -= modelout_xx
            if residuals_analysis[model_name]['x_gp_model']:
                logchi2_gp_model = residuals_analysis[model_name]['x_gp_model']
                parameter_values = residuals_analysis[model_name]['x_gp_parameters']
                log_likelihood += mc.models[logchi2_gp_model].lnlk_compute(parameter_values, x_dataset)
            else:
                log_likelihood += x_dataset.model_logchi2()

            y_dataset.residuals -= modelout_yy
            if residuals_analysis[model_name]['y_gp_model']:
                logchi2_gp_model = residuals_analysis[model_name]['y_gp_model']
                parameter_values = residuals_analysis[model_name]['y_gp_parameters']
                log_likelihood += mc.models[logchi2_gp_model].lnlk_compute(parameter_values, y_dataset)
            else:
                log_likelihood += y_dataset.model_logchi2()


        """ In case there is more than one GP model """
        for logchi2_gp_model in delayed_lnlk_computation:
            log_likelihood += mc.models[logchi2_gp_model].lnlk_compute()

        #""" check for finite log_priors and log_likelihood"""
        #if np.isnan(log_priors) or np.isnan(log_likelihood):
        #    if keep_finite:
        #        return -0.5e8, -0.5e8
        #    else:
        #        return -np.inf, -np.inf

        if not(np.isfinite(log_priors) and np.isfinite(log_likelihood)):
            if keep_finite:
                return -0.5e8, -0.5e8
            else:
                return -np.inf, -np.inf
            #log_likelihood = -np.inf
            #log_priors = -np.inf

        return log_priors, log_likelihood



    mc.emcee_parameters['nwalkers'] = mc.ndim * \
        mc.emcee_parameters['npop_mult']
    if mc.emcee_parameters['nwalkers'] % 2 == 1:
        mc.emcee_parameters['nwalkers'] += 1


    if reloaded_pyde is False:
        """ Run PyDE first"""

        time_start_pyde = time()

        if not os.path.exists(mc.pyde_dir_output):
            os.makedirs(mc.pyde_dir_output)


        print()
        print('*************************************************************')
        print()
        print('Running PyDE')
        print()

        print('Using threading pool for PyDE:',
            mc.pyde_parameters['use_threading_pool'])

        sys.stdout.flush()

        boundaries = np.empty([mc.ndim,2])
        boundaries[:, 0] = 0.
        boundaries[:, 1] = 1.

        if mc.pyde_parameters['use_threading_pool']:
            with multiprocessing.Pool(num_threads) as pool:

                de = DiffEvol(
                    pyde_emcee_call,
                    boundaries,
                    mc.emcee_parameters['nwalkers'],
                    maximize=True,
                    pool=pool)

                de.optimize(int(mc.pyde_parameters['ngen']))
        else:
            de = DiffEvol(
                pyde_emcee_call,
                boundaries,
                mc.emcee_parameters['nwalkers'],
                maximize=True)

            de.optimize(int(mc.pyde_parameters['ngen']))

        population = de.population
        starting_point = np.median(population, axis=0)

        theta_dict = results_analysis.get_theta_dictionary(mc)

        """ bounds redefinition and fix for PyDE anomalous results """
        mc.recenter_bounds_flag = 0.
        if mc.recenter_bounds_flag:
            pyde_save_to_pickle(
                mc, population, starting_point, theta_dict, prefix='orig')

            #mc.recenter_bounds(starting_point)
            population = mc.fix_population(starting_point, population)
            starting_point = np.median(population, axis=0)


            print('Boundaries redefined after PyDE output - CHECK THIS STEP!!')

        pyde_save_to_pickle(mc, population, starting_point, theta_dict)

        time_end_pyde=time()
        print('PyDE completed, it took {0:12.1f} seconds'.format(time_end_pyde-time_start_pyde))
        print()
        sys.stdout.flush()

        if safe_reload:
            print(' safe_reload flag on True, the program will now quit ')
            print(' You have to relaunch it again in order for emcee to work properly')
            print(' No worries, your PyDE results have been saved!')
            quit()


    elif reloaded_emcee is False:
        """ Run emcee second"""

        if not os.path.exists(mc.emcee_dir_output):
            os.makedirs(mc.emcee_dir_output)

        nsteps_todo = mc.emcee_parameters['nsteps']
        state = None
        emcee_skip_check = True

        sampler = emcee.EnsembleSampler(
            mc.emcee_parameters['nwalkers'], mc.ndim, pyde_emcee_call)


        if mc.emcee_parameters['use_threading_pool']:
            with multiprocessing.Pool(num_threads) as pool:
                sampler.pool = pool
                population, prob, state = sampler.run_mcmc(
                    population,
                    nsteps_todo,
                    thin=mc.emcee_parameters['thin'],
                    rstate0=state,
                    progress=True,
                    skip_initial_state_check=emcee_skip_check)

        else:
            print('Warning: NOT using threading, pool, performances will be slower')
            population, prob, state = sampler.run_mcmc(
                population,
                nsteps_todo,
                thin=mc.emcee_parameters['thin'],
                rstate0=state,
                progress=True,
                    skip_initial_state_check=emcee_skip_check)



        theta_dict = results_analysis.get_theta_dictionary(mc)
        emcee_save_to_cpickle(mc, starting_point, population,
                              prob, state, sampler, theta_dict,
                              samples=nsteps_todo)

        flatchain = emcee_flatchain(
            sampler.chain,
            mc.emcee_parameters['nburn'],
            mc.emcee_parameters['thin'])

        flat_lnprob, sampler_lnprob = emcee_flatlnprob(
            sampler.lnprobability, mc.emcee_parameters['nburn'], mc.emcee_parameters['thin'], population, mc.emcee_parameters['nwalkers'])

        results_analysis.print_integrated_ACF(
            sampler.chain,
            theta_dict,
            mc.emcee_parameters['thin'])

        nsample, npams = np.shape(flatchain)
        converted_flatchain = np.zeros_like(flatchain)
        for ii in range(0, nsample):
            converted_flatchain[ii,:] = ultranest_transform(flatchain[ii,:])

        results_analysis.results_summary(mc, converted_flatchain)

        print()

        print()
        print('emcee completed')
        print()
        print('*** writing up warmup file for ultranest ***')

        weights = np.ones((len(flatchain), 1)) / len(flatchain)
        logl = np.zeros(len(flatchain)).reshape((-1, 1))

        labels_array = [None] * len(theta_dictionary)
        for key_name, key_value in theta_dictionary.items():
            labels_array[key_value] = re.sub('_', '-', key_name)

        flatchain_mask00 = (flatchain < 0.001)
        flatchain_mask10 = (flatchain > 0.999)
        flatchain[flatchain_mask00] = 0.001
        flatchain[flatchain_mask10] = 0.999

        np.savetxt(
            mc.emcee_dir_output + '/custom-weighted_post_untransformed.txt',
            np.hstack((weights, logl, flatchain)),
            header=' '.join(['weight', 'logl'] + labels_array),
            fmt='%f'
        )

        print(' You have to relaunch the code if you want to use ultranest with MPI')
        print()
        quit()

    else:

        """ Run ultranest"""
        if not os.path.exists(mc.output_directory):
            os.makedirs(mc.output_directory)

        print('*** emcee results ***')
        print()
        flatchain = emcee_flatchain(
            sampler_chain,
            mc.emcee_parameters['nburn'],
            mc.emcee_parameters['thin'])

        nsample, npams = np.shape(flatchain)
        converted_flatchain = np.zeros_like(flatchain)
        for ii in range(0, nsample):
            converted_flatchain[ii,:] = ultranest_transform(flatchain[ii,:])

        results_analysis.results_summary(mc, converted_flatchain)
        print()
        print()

        posterior_upoints_file = mc.emcee_dir_output + '/custom-weighted_post_untransformed.txt'

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

        from ultranest.integrator import warmstart_from_similar_file

        all_points = np.genfromtxt(posterior_upoints_file)
        upoints = all_points[:,2:]
        sel_min = (upoints <= 0)
        sel_max = (upoints >= 1)
        print('******* ', np.sum(sel_min), np.sum(sel_max))

        mask = np.logical_and(upoints > 0, upoints < 1).all(axis=1)
        assert np.all(mask), (
            'upoints must be between 0 and 1, have:', upoints[~mask,:])

        aux_paramnames, aux_log_likelihood, aux_prior_transform, vectorized = warmstart_from_similar_file(
            posterior_upoints_file, labels_array, ultranest_call, ultranest_transform)

        sampler = ultranest.ReactiveNestedSampler(
            aux_paramnames,
            aux_log_likelihood,
            transform=aux_prior_transform,
            log_dir=mc.output_directory, # folder where to store files
            resume=True, # whether to resume from there (otherwise start from scratch)
        )
        #sampler.stepsampler = ultranest.stepsampler.SliceSampler(
        #    nsteps=nsteps,
        #    generate_direction=ultranest.stepsampler.generate_mixture_random_direction,
        #    # adaptive_nsteps=False,
        #    # max_nsteps=400
        #)
        sampler.run(
            min_num_live_points=nlive,
            dlogz=mc.nested_sampling_parameters['desired_accuracy'], # desired accuracy on logz
            min_ess=mc.nested_sampling_parameters['min_ess'], # number of effective samples
            max_num_improvement_loops=mc.nested_sampling_parameters['improvement_loops'], # how many times to go back and improve
            #frac_remain=0.5
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
