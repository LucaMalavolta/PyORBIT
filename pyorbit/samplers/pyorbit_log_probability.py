from __future__ import print_function

from pyorbit.subroutines.common import np
from pyorbit.classes.model_container_emcee import ModelContainerEmcee
from pyorbit.subroutines.input_parser import pars_input

from pyorbit.subroutines.io_subroutines import pyde_save_to_pickle,\
    pyde_load_from_cpickle,\
    emcee_save_to_cpickle, emcee_load_from_cpickle, emcee_flatchain,\
    emcee_write_dummy_file, starting_point_load_from_cpickle, emcee_simpler_load_from_cpickle

import pyorbit.subroutines.results_analysis as results_analysis
import os
import multiprocessing
import emcee

__all__ = [
    "coerce_theta",
    "prepare_log_probability_model",
    "pyorbit_log_probability",
]


def _output_directories(config_in):
    return {
        "pyde": os.path.join(".", config_in["output"], "pyde") + os.sep,
        "emcee": os.path.join(".", config_in["output"], "emcee") + os.sep,
    }


def _try_load_emcee(emcee_dir_output, explicit=False):
    try:
        mc, _, _, _, _, _, _, _ = emcee_load_from_cpickle(emcee_dir_output)
        return mc
    except FileNotFoundError:
        if explicit:
            raise
        return None


def _try_load_pyde_boundaries(pyde_dir_output, explicit=False):
    try:
        mc, _, _, theta_dict = pyde_load_from_cpickle(pyde_dir_output, prefix="")
        return mc.bounds.copy(), theta_dict
    except FileNotFoundError:
        if explicit:
            raise
        return None, None


def _apply_reloaded_boundaries(mc, previous_boundaries, previous_theta_dict):
    theta_dict = results_analysis.get_theta_dictionary(mc)
    for theta_name, theta_i in theta_dict.items():
        if theta_name in previous_theta_dict:
            mc.bounds[theta_i] = previous_boundaries[previous_theta_dict[theta_name]]


def coerce_theta(mc, theta):
    """
    Convert a vector or a name/value mapping to the theta array used by PyORBIT.
    """

    if isinstance(theta, dict):
        if "theta" in theta and len(theta) == 1:
            theta = theta["theta"]
        else:
            theta_dict = results_analysis.get_theta_dictionary(mc)
            missing = [
                theta_name for theta_name in theta_dict.keys()
                if theta_name not in theta
            ]
            if missing:
                raise ValueError(
                    "Missing theta values for: {0}".format(", ".join(missing))
                )

            unknown = [
                theta_name for theta_name in theta.keys()
                if theta_name not in theta_dict
            ]
            if unknown:
                raise ValueError(
                    "Unknown theta names: {0}".format(", ".join(unknown))
                )

            theta_array = np.zeros(mc.ndim, dtype=np.double)
            for theta_name, theta_i in theta_dict.items():
                theta_array[theta_i] = theta[theta_name]
            theta = theta_array

    theta = np.asarray(theta, dtype=np.double)

    if theta.ndim == 0:
        raise ValueError("theta must be a one-dimensional vector or a matrix")

    if theta.ndim == 1 and theta.size != mc.ndim:
        raise ValueError(
            "theta has {0} values, but this model expects {1}".format(
                theta.size, mc.ndim
            )
        )

    if theta.ndim == 2 and theta.shape[1] != mc.ndim:
        raise ValueError(
            "theta has {0} columns, but this model expects {1}".format(
                theta.shape[1], mc.ndim
            )
        )

    if theta.ndim > 2:
        raise ValueError("theta must be one- or two-dimensional")

    return theta


def prepare_log_probability_model(config_in, input_datasets=None, reload_from="auto"):



    """
    Compute log_probability = log_prior + log_likelihood for theta.
    """


    os.environ["OMP_NUM_THREADS"] = "1"
    try:
        num_threads = int(config_in['parameters'].get('cpu_threads',  multiprocessing.cpu_count()))
    except:
        print(" Something happened when trying to setup multiprocessing, switching back to 1 CPU")
        num_threads = 1

    mp_method = config_in['parameters'].get('mp_method', 'fork')
    multiprocessing.set_start_method(mp_method)

    optimize_dir_output = './' + config_in['output'] + '/optimize/'
    pyde_dir_output = './' + config_in['output'] + '/pyde/'
    emcee_dir_output = './' + config_in['output'] + '/emcee/'

    reloaded_optimize = False
    reloaded_pyde = False
    reloaded_emcee = False

    emcee_skip_check = False

    try:
        mc, population, starting_point, theta_dict = pyde_load_from_cpickle(
            pyde_dir_output, prefix='')
        reloaded_pyde = True
    except FileNotFoundError:
        pass

    try:
        mc, starting_point, population, prob, sampler_chain, \
            sampler_lnprobability, _, theta_dict = \
            emcee_load_from_cpickle(emcee_dir_output)
        reloaded_emcee = True
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
    print('reloaded_emcee: ', reloaded_emcee)
    print()

    if reloaded_pyde or reloaded_emcee:
        previous_boundaries = mc.bounds

    if reloaded_emcee:
        print('Requested steps:', mc.emcee_parameters['nsteps'])

        mc.emcee_parameters['completed_nsteps'] = \
            int(sampler_chain.shape[1] * mc.emcee_parameters['thin'])

        pars_input(config_in, mc, input_datasets, reload_emcee=True)

        """ There's no need to do anything"""
        flatchain = emcee_flatchain(
            sampler_chain, mc.emcee_parameters['nburn'], mc.emcee_parameters['thin'])
        mc.model_setup()
        mc.initialize_logchi2()

        results_analysis.print_bayesian_info(mc)

        results_analysis.print_integrated_ACF(
            sampler_chain, theta_dict, mc.emcee_parameters['thin'])

        """ In case the current startin point comes from a previous analysis """
        mc.emcee_dir_output = emcee_dir_output

        if mc.emcee_parameters['nsteps'] <= mc.emcee_parameters['completed_nsteps']:

            print('Reference Time Tref: ', mc.Tref)
            print()
            print('Dimensions = ', mc.ndim)
            print('Nwalkers = ', mc.emcee_parameters['nwalkers'])
            print('Dimensions = ', mc.ndim)
            print('Steps = ', mc.emcee_parameters['nsteps'])
            print()
            print('Original starting point of emcee:')
            print()

            results_analysis.results_summary(
                mc, starting_point, compute_lnprob=True, is_starting_point=True)

            print('emcee results:')
            print()

            results_analysis.results_summary(mc, flatchain)

            print()
            print('emcee completed')
            print()

            if return_output:
                return mc, sampler_chain, sampler_lnprobability
            else:
                return

    if reloaded_emcee:
        sampled = mc.emcee_parameters['completed_nsteps']
        nsteps_todo = mc.emcee_parameters['nsteps'] \
            - mc.emcee_parameters['completed_nsteps']
    else:

        mc = ModelContainerEmcee()
        pars_input(config_in, mc, input_datasets)

        if mc.pyde_parameters['shutdown_jitter'] or mc.emcee_parameters['shutdown_jitter']:
            for dataset_name, dataset in mc.dataset_dict.items():
                dataset.shutdown_jitter()

        # keep track of which version has been used to perform emcee computations
        mc.emcee_parameters['version'] = emcee.__version__[0]

        mc.model_setup()
        mc.boundaries_setup()
        mc.initialize_logchi2()

        results_analysis.print_bayesian_info(mc)

        mc.pyde_dir_output = pyde_dir_output
        mc.emcee_dir_output = emcee_dir_output

        mc.emcee_parameters['nwalkers'] = mc.ndim * \
            mc.emcee_parameters['npop_mult']
        if mc.emcee_parameters['nwalkers'] % 2 == 1:
            mc.emcee_parameters['nwalkers'] += 1

        if not os.path.exists(mc.emcee_dir_output):
            os.makedirs(mc.emcee_dir_output)

        state = None
        sampled = 0
        nsteps_todo = mc.emcee_parameters['nsteps']

        return mc

def pyorbit_log_probability(
    config_in,
    theta,
    input_datasets=None,
    reload_from="auto",
    return_output=False,
):
    """
    Compute log_probability = log_prior + log_likelihood for theta.
    """

    mc = prepare_log_probability_model(
        config_in, input_datasets=input_datasets, reload_from=reload_from
    )
    theta = coerce_theta(mc, theta)

    global log_priors_likelihood
    def log_priors_likelihood(theta):

        #start = time.time()

        log_priors, log_likelihood = mc.log_priors_likelihood(theta)
        #end = time.time()
        #print("--------Computation took {0:.16f} seconds".format(end - start))

        return log_priors, log_likelihood, log_priors + log_likelihood



    if theta.ndim == 1:
        log_prior, log_likelihood, log_probability = log_priors_likelihood(theta)
    else:
        n_samples = theta.shape[0]
        log_prior = np.zeros(n_samples, dtype=np.double)
        log_likelihood = np.zeros(n_samples, dtype=np.double)
        log_probability = np.zeros(n_samples, dtype=np.double)

        for ii in range(0, n_samples):
            log_prior[ii], log_likelihood[ii], log_probability[ii] = (
                log_priors_likelihood(theta[ii, :])
            )

    output = {
        "log_prior": log_prior,
        "log_likelihood": log_likelihood,
        "log_probability": log_probability,
        "theta": theta,
        "theta_dictionary": results_analysis.get_theta_dictionary(mc),
    }

    if return_output:
        return mc, output

    return output
