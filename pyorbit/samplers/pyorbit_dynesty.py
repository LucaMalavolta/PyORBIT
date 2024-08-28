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


def pyorbit_dynesty(config_in, input_datasets=None, return_output=None, run_nested=True):

    multiprocessing.set_start_method('fork')

    reloaded_dynesty = False
    output_directory = './' + config_in['output'] + '/dynesty/'
    save_checkpoint = output_directory + 'dynesty.save'
    save_checkpoint_maxevidence = output_directory + 'dynesty_maxevidence.save'

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    try:
        mc = nested_sampling_load_from_cpickle(output_directory)
        reloaded_dynesty = True

        pars_input(config_in, mc, input_datasets, reload_emcee=True)
        mc.output_directory = output_directory

    except:
        mc = ModelContainerDynesty()

        pars_input(config_in, mc, input_datasets)
        mc.output_directory = output_directory

    try:
        results = dynesty_results_load_from_cpickle(mc.output_directory)
        print('Dynesty results already saved in the respective directory, run pyorbit_results')
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

    if reloaded_dynesty:
        mc.model_setup()
        #mc.boundaries_setup()
        mc.initialize_logchi2()
    else:
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
        return mc.dynesty_priors(cube)

    global dynesty_loglikelihood
    def dynesty_loglikelihood(theta):
        return mc.dynesty_call(theta)



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

    if not reloaded_dynesty:
        """ A dummy file is created to let the cpulimit script to proceed with the next step"""
        nested_sampling_write_dummy_file(mc)
        nested_sampling_save_to_cpickle(mc)

    if pfrac_value > 1.0 or  pfrac_value < 0.0:
        print('Double run with 0/100 posterior/evidence split; then 100/0 posterior/evidence split')
        print()
        pfrac_value = 0.000

        print('Setting up the Dynamic Nested Sampling, posterior/evidence split = {0:4.3f}'.format(pfrac))
        print()

        if use_threading_pool:


            with multiprocessing.Pool(num_threads) as pool:

                # "Dynamic" nested sampling.

                try:
                    dsampler_maxevidence = dynesty.DynamicNestedSampler.restore(save_checkpoint_maxevidence, pool=pool)
                    print('Restoring Dynamic Nested Sampling')
                    if run_nested:
                        dsampler_maxevidence.run_nested(resume=True)
                except:

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
                    dsampler_maxevidence.run_nested(checkpoint_file=save_checkpoint_maxevidence, dlogz=dlogz, wt_kwargs={'pfrac': pfrac_value})
                    print()



        else:

            try:
                dsampler_maxevidence = dynesty.DynamicNestedSampler.restore(save_checkpoint_maxevidence)
                print('Restoring Dynamic Nested Sampling for MaxEvidence')
                results_maxevidence = dsampler_maxevidence.results
                dynesty_results_maxevidence_save_to_cpickle(mc.output_directory, results_maxevidence)
                if run_nested:
                    dsampler_maxevidence.run_nested(resume=True)
            except:


                dsampler_maxevidence = dynesty.DynamicNestedSampler(dynesty_loglikelihood,
                                                        dynesty_priors,
                                                        mc.ndim,
                                                        nlive=nlive,
                                                        bound= mc.nested_sampling_parameters['bound'],
                                                        sample= mc.nested_sampling_parameters['sample'],
                                                        )
                print('Running Dynamic Nested Sampling')
                dsampler_maxevidence.run_nested(checkpoint_file=save_checkpoint_maxevidence, dlogz=dlogz, wt_kwargs={'pfrac': pfrac_value})
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
        print('                                        initial stopping criterion = {0:5.4f}'.format(dlogz))
        print('                                        bound = ', mc.nested_sampling_parameters['bound'])
        print('                                        sample = ', mc.nested_sampling_parameters['sample'])
        print()

    if use_threading_pool:

        with multiprocessing.Pool(num_threads) as pool:

            try:
                dsampler = dynesty.DynamicNestedSampler.restore(save_checkpoint, pool=pool)
                print('Restoring Dynamic Nested Sampling')
                results = dsampler.results
                dynesty_results_save_to_cpickle(output_directory, results)
                if run_nested:
                    dsampler.run_nested(resume=True)
            except:

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
                    dsampler.run_nested(checkpoint_file=save_checkpoint)
                else:
                    dsampler.run_nested(checkpoint_file=save_checkpoint, dlogz_init=dlogz, wt_kwargs={'pfrac': pfrac_value})
                print()
    else:
        try:
            dsampler = dynesty.DynamicNestedSampler.restore(save_checkpoint)
            print('Restoring Dynamic Nested Sampling')
            if run_nested:
                dsampler.run_nested(resume=True)
        except:
            dsampler = dynesty.DynamicNestedSampler(dynesty_loglikelihood,
                                                    dynesty_priors,
                                                    mc.ndim,
                                                    nlive=nlive,
                                                    bound= mc.nested_sampling_parameters['bound'],
                                                    sample= mc.nested_sampling_parameters['sample'])

            print('Running Dynamic Nested Sampling without parallelization')
            if use_default:
                dsampler.run_nested(checkpoint_file=save_checkpoint)
            else:
                dsampler.run_nested(checkpoint_file=save_checkpoint, dlogz_init=dlogz, wt_kwargs={'pfrac': pfrac_value})
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
    #nested_sampling_write_dummy_file(mc)
    #nested_sampling_save_to_cpickle(mc)
    dynesty_results_save_to_cpickle(mc.output_directory, results)

    if return_output:
        return mc
    else:
        return
