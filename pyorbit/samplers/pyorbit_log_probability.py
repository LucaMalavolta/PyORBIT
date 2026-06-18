from __future__ import print_function

from pyorbit.subroutines.common import np
from pyorbit.classes.model_container_emcee import ModelContainerEmcee
from pyorbit.subroutines.input_parser import pars_input
from pyorbit.subroutines.io_subroutines import (
    emcee_load_from_cpickle,
    pyde_load_from_cpickle,
)
import pyorbit.subroutines.results_analysis as results_analysis
import os

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


def prepare_log_probability_model(config_in, input_datasets=None, reload_from="auto"):
    """
    Build the model container needed to evaluate log_prior + log_likelihood.

    Parameters
    ----------
    config_in : dict
        Parsed PyORBIT configuration.
    input_datasets : optional
        Optional in-memory datasets, forwarded to ``pars_input``.
    reload_from : {"auto", "none", "emcee", "pyde"}
        Source used to recover saved state. ``auto`` prefers a saved emcee
        container, then pyde boundaries, then a fresh YAML initialization.
    """

    if reload_from not in ["auto", "none", "emcee", "pyde"]:
        raise ValueError(
            "reload_from must be one of: auto, none, emcee, pyde"
        )

    os.environ["OMP_NUM_THREADS"] = "1"

    output_directories = _output_directories(config_in)
    source = "config"

    if reload_from in ["auto", "emcee"]:
        mc = _try_load_emcee(
            output_directories["emcee"], explicit=(reload_from == "emcee")
        )
        if mc is not None:
            pars_input(config_in, mc, input_datasets, reload_emcee=True)
            mc.model_setup()
            mc.initialize_logchi2()
            mc.emcee_dir_output = output_directories["emcee"]
            return mc, "emcee"

    previous_boundaries = None
    previous_theta_dict = None
    if reload_from in ["auto", "pyde"]:
        previous_boundaries, previous_theta_dict = _try_load_pyde_boundaries(
            output_directories["pyde"], explicit=(reload_from == "pyde")
        )
        if previous_boundaries is not None:
            source = "pyde"

    mc = ModelContainerEmcee()
    pars_input(config_in, mc, input_datasets)

    if mc.pyde_parameters["shutdown_jitter"] or mc.emcee_parameters["shutdown_jitter"]:
        for dataset_name, dataset in mc.dataset_dict.items():
            dataset.shutdown_jitter()

    mc.model_setup()
    mc.boundaries_setup()

    if previous_boundaries is not None:
        _apply_reloaded_boundaries(mc, previous_boundaries, previous_theta_dict)

    mc.initialize_logchi2()
    mc.pyde_dir_output = output_directories["pyde"]
    mc.emcee_dir_output = output_directories["emcee"]

    return mc, source


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


def _evaluate_log_probability(mc, theta):
    log_prior, log_likelihood = mc.log_priors_likelihood(theta)
    return log_prior, log_likelihood, log_prior + log_likelihood


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

    mc, source = prepare_log_probability_model(
        config_in, input_datasets=input_datasets, reload_from=reload_from
    )
    theta = coerce_theta(mc, theta)

    if theta.ndim == 1:
        log_prior, log_likelihood, log_probability = _evaluate_log_probability(
            mc, theta
        )
    else:
        n_samples = theta.shape[0]
        log_prior = np.zeros(n_samples, dtype=np.double)
        log_likelihood = np.zeros(n_samples, dtype=np.double)
        log_probability = np.zeros(n_samples, dtype=np.double)

        for ii in range(0, n_samples):
            log_prior[ii], log_likelihood[ii], log_probability[ii] = (
                _evaluate_log_probability(mc, theta[ii, :])
            )

    output = {
        "log_prior": log_prior,
        "log_likelihood": log_likelihood,
        "log_probability": log_probability,
        "model_source": source,
        "theta": theta,
        "theta_dictionary": results_analysis.get_theta_dictionary(mc),
    }

    if return_output:
        return mc, output

    return output
