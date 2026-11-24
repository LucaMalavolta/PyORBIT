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

import pyorbit.common.gpu_backend as gpu_backend
import pyorbit.samplers.dynesty_pool as dynesty_pool
from contextlib import contextmanager

__all__ = ["pyorbit_dynesty"]


@contextmanager
def _gpu_pool_ctx(gpu_workers, config_in, input_datasets, mc, vram_fraction=None):
    """Spawn K dynesty workers that share the GPU (each VRAM-capped) so the
    sampler's independent live points are evaluated concurrently.

    dynesty's pool uses ``mp.get_context('spawn')`` -- each worker is a fresh
    interpreter that builds its own CUDA context (no fork-after-CUDA). The
    per-process VRAM cap is set *before* the workers spawn (they inherit the
    environment), and the master is told to reuse its already-built container
    for the in-master prior transform so it does not open a second JAX context.
    """
    os.environ['XLA_PYTHON_CLIENT_PREALLOCATE'] = 'false'
    if vram_fraction is None:
        vram_fraction = min(0.30, max(0.05, 0.90 / (gpu_workers + 1)))
    os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = str(vram_fraction)
    dynesty_pool.set_master_mc(mc)
    import dynesty
    pool = dynesty.pool.Pool(gpu_workers,
                             dynesty_pool.pool_loglike,
                             dynesty_pool.pool_priors,
                             logl_args=(config_in, input_datasets),
                             ptform_args=(config_in, input_datasets))
    with pool:
        yield pool

"""
def show(filepath):
    # open the output (pdf) file for the user
    if os.name == 'mac': subprocess.call(('open', filepath))
    elif os.name == 'nt': os.startfile(filepath)
"""


def pyorbit_dynesty(config_in, input_datasets=None, return_output=None, run_nested=True):

    mp_method = config_in['parameters'].get('mp_method', 'fork')
    multiprocessing.set_start_method(mp_method)

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
        mc.print_model_info()
        #mc.boundaries_setup()
        mc.initialize_logchi2()
    else:
        mc.model_setup()
        mc.print_model_info()
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
    maxcall = mc.nested_sampling_parameters.get('maxcall', None)
    use_threading_pool = mc.nested_sampling_parameters.get('use_threading_pool', True)

    # Step 2: fill the GPU with K concurrent workers. ``gpu_workers > 1``
    # enables dynesty's spawn-based pool (CUDA-safe: fresh interpreter per
    # worker, each builds its own CUDA context -- no fork-after-CUDA), so the
    # sampler's independent live points are evaluated concurrently and the GPU
    # / VRAM stay busy instead of idling between per-point round-trips.
    # Default 1 = unchanged single-process behaviour.
    gpu_workers = int(config_in['parameters'].get('gpu_workers', 1))
    gpu_vram_fraction = config_in['parameters'].get('gpu_vram_fraction', None)
    use_worker_pool = gpu_workers > 1

    if gpu_backend.gpu_enabled() and use_threading_pool and not use_worker_pool:
        # Forking a worker after the parent has used CUDA is unsafe, and the
        # spawn pool is not requested here, so fall back to the single-process
        # path (both 'use_threading_pool' branches below take the no-pool route).
        print('GPU enabled: single-process (use_threading_pool is ignored).')
        print()
        use_threading_pool = False

    # Precedence: spawn worker pool (K processes) > CPU fork pool > single-process.
    if use_worker_pool:
        use_threading_pool = False

    pfrac_value = mc.nested_sampling_parameters.get('pfrac', 0.000)
    use_default = (mc.nested_sampling_parameters.get('default', False) | mc.nested_sampling_parameters.get('use_default', False))

    if use_default:
        pfrac_value = 0.80

    if use_worker_pool:
        print('Using {0}-worker dynesty pool (per-worker VRAM fraction {1})'.format(
            gpu_workers, gpu_vram_fraction if gpu_vram_fraction is not None else 'auto'))
    else:
        print('Using threading pool for dynesty:', use_threading_pool)

    if not run_nested:
        print()
        print('*************************************************************')
        print()
        print('Dry run requested: validating configuration without sampling')
        print()
        print('*************************************************************')
        print()

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

        """ Constructing the sampler (without a pool, without a checkpoint file) validates
        nlive/bound/sample eagerly and cheaply: dynesty raises immediately on an unsupported
        bounding method or sampling method, and no likelihood evaluation happens until
        run_nested() is actually called. """
        dynesty.DynamicNestedSampler(dynesty_loglikelihood,
                                      dynesty_priors,
                                      mc.ndim,
                                      nlive=nlive,
                                      bound=mc.nested_sampling_parameters['bound'],
                                      sample=mc.nested_sampling_parameters['sample'])

        theta = dynesty_priors(0.5 * np.ones(mc.ndim))
        results_analysis.results_summary(
            mc, theta, compute_lnprob=True, is_starting_point=True)

        print('Dry run completed: setup validated, sampling not performed')
        print()

        if return_output:
            return mc
        else:
            return

    if not reloaded_dynesty:
        """ A dummy file is created to let the cpulimit script to proceed with the next step"""
        nested_sampling_write_dummy_file(mc)
        nested_sampling_save_to_cpickle(mc)

    if pfrac_value > 1.0 or  pfrac_value < 0.0:
        print('Double run with 0/100 posterior/evidence split; then 100/0 posterior/evidence split')
        print()
        pfrac_value = 0.000

        print('Setting up the Dynamic Nested Sampling, posterior/evidence split = {0:4.3f}'.format(pfrac_value))
        print()

        if use_worker_pool:

            with _gpu_pool_ctx(gpu_workers, config_in, input_datasets, mc, gpu_vram_fraction) as pool:

                # "Dynamic" nested sampling.

                try:
                    dsampler_maxevidence = dynesty.DynamicNestedSampler.restore(save_checkpoint_maxevidence, pool=pool)
                    print('Restoring Dynamic Nested Sampling')
                except Exception as err:
                    if os.path.exists(save_checkpoint_maxevidence):
                        corrupt_checkpoint = save_checkpoint_maxevidence + '.corrupt'
                        os.replace(save_checkpoint_maxevidence, corrupt_checkpoint)
                        print('WARNING: could not restore checkpoint {0} ({1}); moved aside to {2} and starting a fresh run'.format(save_checkpoint_maxevidence, err, corrupt_checkpoint))

                    dsampler_maxevidence = dynesty.DynamicNestedSampler(pool.loglike,
                                                            pool.prior_transform,
                                                            mc.ndim,
                                                            nlive=nlive,
                                                            pool=pool,
                                                            queue_size=gpu_workers,
                                                            bound= mc.nested_sampling_parameters['bound'],
                                                            sample= mc.nested_sampling_parameters['sample'],
                                                            use_pool={
                                                                'prior_transform': False},
                                                            )
                    print('Running Dynamic Nested Sampling')
                    dsampler_maxevidence.run_nested(checkpoint_file=save_checkpoint_maxevidence, dlogz_init=dlogz, maxcall=maxcall, wt_kwargs={'pfrac': pfrac_value})
                    print()
                else:
                    results_maxevidence = dsampler_maxevidence.results
                    dynesty_results_maxevidence_save_to_cpickle(mc.output_directory, results_maxevidence)
                    dsampler_maxevidence.run_nested(resume=True)



        elif use_threading_pool:


            with multiprocessing.Pool(num_threads) as pool:

                # "Dynamic" nested sampling.

                try:
                    dsampler_maxevidence = dynesty.DynamicNestedSampler.restore(save_checkpoint_maxevidence, pool=pool)
                    print('Restoring Dynamic Nested Sampling')
                except Exception as err:
                    if os.path.exists(save_checkpoint_maxevidence):
                        corrupt_checkpoint = save_checkpoint_maxevidence + '.corrupt'
                        os.replace(save_checkpoint_maxevidence, corrupt_checkpoint)
                        print('WARNING: could not restore checkpoint {0} ({1}); moved aside to {2} and starting a fresh run'.format(save_checkpoint_maxevidence, err, corrupt_checkpoint))

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
                    dsampler_maxevidence.run_nested(checkpoint_file=save_checkpoint_maxevidence, dlogz_init=dlogz, maxcall=maxcall, wt_kwargs={'pfrac': pfrac_value})
                    print()
                else:
                    results_maxevidence = dsampler_maxevidence.results
                    dynesty_results_maxevidence_save_to_cpickle(mc.output_directory, results_maxevidence)
                    dsampler_maxevidence.run_nested(resume=True)



        else:

            try:
                dsampler_maxevidence = dynesty.DynamicNestedSampler.restore(save_checkpoint_maxevidence)
                print('Restoring Dynamic Nested Sampling for MaxEvidence')
            except Exception as err:
                if os.path.exists(save_checkpoint_maxevidence):
                    corrupt_checkpoint = save_checkpoint_maxevidence + '.corrupt'
                    os.replace(save_checkpoint_maxevidence, corrupt_checkpoint)
                    print('WARNING: could not restore checkpoint {0} ({1}); moved aside to {2} and starting a fresh run'.format(save_checkpoint_maxevidence, err, corrupt_checkpoint))

                dsampler_maxevidence = dynesty.DynamicNestedSampler(dynesty_loglikelihood,
                                                        dynesty_priors,
                                                        mc.ndim,
                                                        nlive=nlive,
                                                        bound= mc.nested_sampling_parameters['bound'],
                                                        sample= mc.nested_sampling_parameters['sample'],
                                                        )
                print('Running Dynamic Nested Sampling')
                dsampler_maxevidence.run_nested(checkpoint_file=save_checkpoint_maxevidence, dlogz_init=dlogz, maxcall=maxcall, wt_kwargs={'pfrac': pfrac_value})
                print()
            else:
                results_maxevidence = dsampler_maxevidence.results
                dynesty_results_maxevidence_save_to_cpickle(mc.output_directory, results_maxevidence)
                dsampler_maxevidence.run_nested(resume=True)

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

    if use_worker_pool:

        with _gpu_pool_ctx(gpu_workers, config_in, input_datasets, mc, gpu_vram_fraction) as pool:

            try:
                dsampler = dynesty.DynamicNestedSampler.restore(save_checkpoint, pool=pool)
                print('Restoring Dynamic Nested Sampling')
            except Exception as err:
                if os.path.exists(save_checkpoint):
                    corrupt_checkpoint = save_checkpoint + '.corrupt'
                    os.replace(save_checkpoint, corrupt_checkpoint)
                    print('WARNING: could not restore checkpoint {0} ({1}); moved aside to {2} and starting a fresh run'.format(save_checkpoint, err, corrupt_checkpoint))

                dsampler = dynesty.DynamicNestedSampler(pool.loglike,
                                                        pool.prior_transform,
                                                        mc.ndim,
                                                        nlive=nlive,
                                                        pool=pool,
                                                        queue_size=gpu_workers,
                                                        bound= mc.nested_sampling_parameters['bound'],
                                                        sample= mc.nested_sampling_parameters['sample'],
                                                        use_pool={
                                                            'prior_transform': False},
                                                        )

                print('Running Dynamic Nested Sampling')
                if use_default:
                    dsampler.run_nested(checkpoint_file=save_checkpoint, maxcall=maxcall)
                else:
                    dsampler.run_nested(checkpoint_file=save_checkpoint, dlogz_init=dlogz, maxcall=maxcall, wt_kwargs={'pfrac': pfrac_value})
                print()
            else:
                results = dsampler.results
                dynesty_results_save_to_cpickle(mc.output_directory, results)
                dsampler.run_nested(resume=True)

    elif use_threading_pool:

        with multiprocessing.Pool(num_threads) as pool:

            try:
                dsampler = dynesty.DynamicNestedSampler.restore(save_checkpoint, pool=pool)
                print('Restoring Dynamic Nested Sampling')
            except Exception as err:
                if os.path.exists(save_checkpoint):
                    corrupt_checkpoint = save_checkpoint + '.corrupt'
                    os.replace(save_checkpoint, corrupt_checkpoint)
                    print('WARNING: could not restore checkpoint {0} ({1}); moved aside to {2} and starting a fresh run'.format(save_checkpoint, err, corrupt_checkpoint))

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
                    dsampler.run_nested(checkpoint_file=save_checkpoint, maxcall=maxcall)
                else:
                    dsampler.run_nested(checkpoint_file=save_checkpoint, dlogz_init=dlogz, maxcall=maxcall, wt_kwargs={'pfrac': pfrac_value})
                print()
            else:
                results = dsampler.results
                dynesty_results_save_to_cpickle(mc.output_directory, results)
                dsampler.run_nested(resume=True)
    else:
        try:
            dsampler = dynesty.DynamicNestedSampler.restore(save_checkpoint)
            print('Restoring Dynamic Nested Sampling')
        except Exception as err:
            if os.path.exists(save_checkpoint):
                corrupt_checkpoint = save_checkpoint + '.corrupt'
                os.replace(save_checkpoint, corrupt_checkpoint)
                print('WARNING: could not restore checkpoint {0} ({1}); moved aside to {2} and starting a fresh run'.format(save_checkpoint, err, corrupt_checkpoint))

            dsampler = dynesty.DynamicNestedSampler(dynesty_loglikelihood,
                                                    dynesty_priors,
                                                    mc.ndim,
                                                    nlive=nlive,
                                                    bound= mc.nested_sampling_parameters['bound'],
                                                    sample= mc.nested_sampling_parameters['sample'])

            print('Running Dynamic Nested Sampling without parallelization')
            if use_default:
                dsampler.run_nested(checkpoint_file=save_checkpoint, maxcall=maxcall)
            else:
                dsampler.run_nested(checkpoint_file=save_checkpoint, dlogz_init=dlogz, maxcall=maxcall, wt_kwargs={'pfrac': pfrac_value})
            print()
        else:
            results = dsampler.results
            dynesty_results_save_to_cpickle(mc.output_directory, results)
            dsampler.run_nested(resume=True)

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

    try:
        os.remove(save_checkpoint)
        os.remove(save_checkpoint_maxevidence)
    except:
        pass
    
    if return_output:
        return mc
    else:
        return
