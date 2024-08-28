from __future__ import print_function
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
from pyorbit.subroutines.common import np, nested_sampling_prior_compute

__all__ = ["pyorbit_dynesty"]

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
    save_checkpoint = mc.output_directory + 'dynesty.save'

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

    os.environ["OMP_NUM_THREADS"] = "1"
    try:
        num_threads = int(config_in['parameters'].get('cpu_threads',  multiprocessing.cpu_count()))
    except:
        print(" Something happened when trying to setup multiprocessing, switching back to 1 CPU")
        num_threads = 1

    mc.model_setup()
    mc.boundaries_setup()
    mc.initialize_logchi2()

    mc.starting_points_setup()

    results_analysis.print_bayesian_info(mc)


    if 'nlive_mult' in mc.nested_sampling_parameters:
        nlive = mc.ndim * mc.nested_sampling_parameters['nlive_mult']
    else:
        nlive = mc.nested_sampling_parameters['nlive']


    global dynesty_priors
    def dynesty_priors(cube):
        theta = np.zeros(len(cube), dtype=np.double)

        for i in range(0, len(cube)):
            theta[i] = nested_sampling_prior_compute(
                cube[i], mc.priors[i][0], mc.priors[i][2], mc.spaces[i])
        return theta

    global dynesty_loglikelihood
    def dynesty_loglikelihood(theta):

        log_likelihood = 0.00

        """
        Constant term added either by dataset.model_logchi2() or gp.log_likelihood()
        """
        if not mc.check_bounds(theta):
            return -np.inf

        if mc.dynamical_model is not None:
            """ check if any keyword ahas get the output model from the dynamical tool
            we must do it here because all the planet are involved"""
            dynamical_output = mc.dynamical_model.compute(mc, theta)

        delayed_lnlk_computation = []
        residuals_analysis = {}

        for dataset_name, dataset in mc.dataset_dict.items():

            logchi2_gp_model = None
            compute_gp_residuals = False

            dataset.model_reset()
            parameter_values = dataset.convert(theta)
            dataset.compute(parameter_values)

            if 'none' in dataset.models or 'None' in dataset.models:
                continue
            if not dataset.models:
                continue

            skip_loglikelihood = False

            for model_name in dataset.models:

                parameter_values = {}
                for common_ref in mc.models[model_name].common_ref:
                    parameter_values.update(
                        mc.common_models[common_ref].convert(theta))

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
            it must be coomputed after any other model has been removed from the
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

        """ check for finite log_likelihood"""
        if  np.isnan(log_likelihood):
            log_likelihood = -np.inf

        return log_likelihood




    print('Number of live points:', nlive)
    print()
    print('number of multiprocessing threads:', num_threads)

    print()
    print('Reference Time Tref: ', mc.Tref)
    print()
    print('*************************************************************')
    print()

    try:
        import dynesty
    except (ModuleNotFoundError,ImportError):
        print("ERROR: dynesty not installed, this will not work")
        quit()


    dlogz = mc.nested_sampling_parameters.get('dlogz', 0.01)
    use_threading_pool = mc.nested_sampling_parameters.get('use_threading_pool', True)
    pfrac_value = mc.nested_sampling_parameters.get('pfrac', 0.000)
    use_default = (mc.nested_sampling_parameters.get('default', False) | mc.nested_sampling_parameters.get('use_default', False))

    if use_default:
        pfrac_value = 0.80

    print('Using threading pool for dynesty:', use_threading_pool)


    if pfrac_value > 1.0 or  pfrac_value < 0.0:
        print('Double run with 0/100 posterior/evidence split; then 100/0 posterior/evidence split')
        print()
        pfrac_value = 0.000

        print('Setting up the Dynamic Nested Sampling, posterior/evidence split = {0:4.3f}'.format(pfrac))
        print()

        if use_threading_pool:


            with multiprocessing.Pool(num_threads) as pool:

                # "Dynamic" nested sampling.

                dsampler_maxevidence = dynesty.DynamicNestedSampler(dynesty_loglikelihood,
                                                        dynesty_priors,
                                                        mc.ndim,
                                                        nlive=nlive,
                                                        pool=pool,
                                                        queue_size=num_threads,
                                                        bound= mc.nested_sampling_parameters['bound'],
                                                        sample= mc.nested_sampling_parameters['sample'],
                                                        use_pool={
                                                            'prior_transform': False},
                                                        )
                print('Running Dynamic Nested Sampling')
                dsampler_maxevidence.run_nested(dlogz=dlogz, wt_kwargs={'pfrac': pfrac_value})
                print()
        else:

            dsampler_maxevidence = dynesty.DynamicNestedSampler(dynesty_loglikelihood,
                                                    dynesty_priors,
                                                    mc.ndim,
                                                    nlive=nlive,
                                                    bound= mc.nested_sampling_parameters['bound'],
                                                    sample= mc.nested_sampling_parameters['sample'],
                                                    )
            print('Running Dynamic Nested Sampling')
            dsampler_maxevidence.run_nested(dlogz=dlogz, wt_kwargs={'pfrac': pfrac_value})
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

        pfrac_value = 1.000

    if use_default:
        print('Setting up the Dynamic Nested Sampling, number of live points = {0:6.0f}'.format(nlive))
        print('Using default Dynesty hyperparameters')
        print()
    else:
        print('Setting up the Dynamic Nested Sampling, number of live points = {0:6.0f}'.format(nlive))
        print('                                        posterior/evidence split = {0:4.3f}'.format(pfrac_value))
        print('                                        inizial stopping criterion = {0:5.4f}'.format(dlogz))
        print()

    if use_threading_pool:

        with multiprocessing.Pool(num_threads) as pool:
            dsampler = dynesty.DynamicNestedSampler(dynesty_loglikelihood,
                                                    dynesty_priors,
                                                    mc.ndim,
                                                    nlive=nlive,
                                                    pool=pool,
                                                    queue_size=num_threads,
                                                    bound= mc.nested_sampling_parameters['bound'],
                                                    sample= mc.nested_sampling_parameters['sample'],
                                                    use_pool={
                                                        'prior_transform': False},
                                                    )

            print('Running Dynamic Nested Sampling')
            if use_default:
                dsampler.run_nested()
            else:
                dsampler.run_nested(dlogz_init=dlogz, wt_kwargs={'pfrac': pfrac_value})
            print()
    else:

        dsampler = dynesty.DynamicNestedSampler(dynesty_loglikelihood,
                                                dynesty_priors,
                                                mc.ndim,
                                                nlive=nlive,
                                                bound= mc.nested_sampling_parameters['bound'],
                                                sample= mc.nested_sampling_parameters['sample'])

        print('Running Dynamic Nested Sampling without parallelization')
        if use_default:
            dsampler.run_nested()
        else:
            dsampler.run_nested(dlogz_init=dlogz, wt_kwargs={'pfrac': pfrac_value})
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
