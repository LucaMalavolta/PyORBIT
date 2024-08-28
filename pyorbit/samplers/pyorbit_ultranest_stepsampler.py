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
import matplotlib.pyplot as plt
from pyorbit.subroutines.common import np, nested_sampling_prior_compute

__all__ = ["pyorbit_ultranest_stepsampler"]

"""
def show(filepath):
    # open the output (pdf) file for the user
    if os.name == 'mac': subprocess.call(('open', filepath))
    elif os.name == 'nt': os.startfile(filepath)
"""


def pyorbit_ultranest_stepsampler(config_in, input_datasets=None, return_output=None):

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
        import ultranest
        import ultranest.stepsampler
    except (ModuleNotFoundError,ImportError):
        print("ERROR: ultranest not installed, this will not work")
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

        log_likelihood = 0.00

        """
        Constant term added either by dataset.model_logchi2() or gp.log_likelihood()
        """
        if not mc.check_bounds(theta):
            return -0.5e10

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
        #if  np.isnan(log_likelihood):
        #    log_likelihood = -np.inf

        if  np.isnan(log_likelihood):
            log_likelihood = -0.5e10


        return log_likelihood

    nsteps = 2 * mc.ndim

    sampler = ultranest.ReactiveNestedSampler(
        labels_array,
        mc.ultranest_call,
        transform=mc.ultranest_transform,
        log_dir=mc.output_directory, # folder where to store files
        resume=True, # whether to resume from there (otherwise start from scratch)
    )
    sampler.stepsampler = ultranest.stepsampler.SliceSampler(
        nsteps=nsteps,
        generate_direction=ultranest.stepsampler.generate_mixture_random_direction,
        # adaptive_nsteps=False,
        # max_nsteps=400
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
