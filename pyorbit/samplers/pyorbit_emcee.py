from __future__ import print_function
from pyorbit.subroutines.common import np
from pyorbit.classes.model_container_emcee import ModelContainerEmcee
from pyorbit.subroutines.input_parser import yaml_parser, pars_input
from pyorbit.subroutines.io_subroutines import (
    pyde_save_to_pickle,
    pyde_load_from_cpickle,
    emcee_save_to_cpickle,
    emcee_load_from_cpickle,
    emcee_flatchain,
    emcee_write_dummy_file,
    starting_point_load_from_cpickle,
    emcee_simpler_load_from_cpickle,
)
import pyorbit.subroutines.results_analysis as results_analysis
import os
import sys
import multiprocessing
from time import time

import pyorbit.common.gpu_backend as gpu_backend
# __all__ = ["pyorbit_emcee"]


def _fbatch_2d_results(mc, theta):
    """F-batch: per-walker (lnprob, log_priors) results for the emcee
    vectorize=True 2-D (nwalkers, ndim) ``log_prob_fn`` call.

    Cardinal rule (the batch must never change a fit, and the recorded blob
    must be bit-equal to the per-theta wrapper on EVERY row):

    * finite rows keep the batched (lnprob, log_priors);
    * every non-finite row resolves through the EXACT per-theta path for
      that row, and its blob (the per-theta log_priors) is what gets
      recorded. This is the F1 fix: for a finite-prior row whose batched
      likelihood is -inf, per-theta mode records the FINITE log_priors --
      the old ``blobs[bad] = -inf`` masking diverged from per-theta mode on
      exactly that class of rows;
    * the lnprob for a non-finite row is the batched value (non-finite) in
      both modes -- it is never flipped to finite here, so emcee's
      accept/reject stream is untouched.

    Returns the list of ``(lnprob, log_priors)`` tuples emcee 3's
    ``compute_log_prob`` expects (one per walker).
    """
    lnprobs = np.asarray(mc.lnprob_batch(theta), dtype=np.double)
    blobs = np.asarray(mc._lnprob_batch_priors(theta), dtype=np.double)
    bad = ~np.isfinite(lnprobs)
    if bad.any():
        for i in np.nonzero(bad)[0]:
            log_priors_i, _log_likelihood_i = mc.log_priors_likelihood(theta[i])
            blobs[i] = float(log_priors_i)
    return list(zip(lnprobs.tolist(), blobs.tolist()))


def _emcee_vectorize_gate(version_str):
    """F-batch emcee major-version gate (FIX 2).

    The 2-D vectorized contract used by :func:`_fbatch_2d_results` is only
    validated against emcee 3.x (3.1.6). The decision is a pure function of
    the version string so it can be unit-checked with any fake version.

    Returns ``(ok, message)``: ``ok`` is True only for the validated major;
    ``message`` is the single log line emitted (one line) when gating off.
    Never raises on a malformed version -- an unparseable version gates off.
    """
    try:
        major = int(str(version_str).strip().split(".")[0])
    except Exception:
        return False, (
            "emcee version '{0}' not validated for the 2-D vectorized "
            "path; falling back to per-theta".format(version_str or "unknown")
        )
    if major == 3:
        return True, ""
    return False, (
        "emcee {0} not validated for the 2-D vectorized path; "
        "falling back to per-theta".format(version_str)
    )


def pyorbit_emcee(config_in, input_datasets=None, return_output=None):

    try:
        import emcee
    except:
        print("ERROR: emcee not installed, this will not work")
        quit()

    moves_definition = {
        "stretch": emcee.moves.StretchMove(),
        "differential_evolution": emcee.moves.DEMove(),
        "de_snooker": emcee.moves.DESnookerMove(),
        "de": emcee.moves.DEMove(),
        "redblue": emcee.moves.RedBlueMove(),
        "walk": emcee.moves.WalkMove(),
    }

    os.environ["OMP_NUM_THREADS"] = "1"
    try:
        num_threads = int(
            config_in["parameters"].get("cpu_threads", multiprocessing.cpu_count())
        )
    except:
        print(
            " Something happened when trying to setup multiprocessing, switching back to 1 CPU"
        )
        num_threads = 1

    mp_method = config_in["parameters"].get("mp_method", "fork")
    multiprocessing.set_start_method(mp_method)

    optimize_dir_output = "./" + config_in["output"] + "/optimize/"
    pyde_dir_output = "./" + config_in["output"] + "/pyde/"
    emcee_dir_output = "./" + config_in["output"] + "/emcee/"

    reloaded_optimize = False
    reloaded_pyde = False
    reloaded_emcee = False

    emcee_skip_check = False

    try:
        mc, population, starting_point, theta_dict = pyde_load_from_cpickle(
            pyde_dir_output, prefix=""
        )
        reloaded_pyde = True
    except FileNotFoundError:
        pass

    try:
        (
            mc,
            starting_point,
            population,
            prob,
            sampler_chain,
            sampler_lnprobability,
            _,
            theta_dict,
            prior,
        ) = emcee_load_from_cpickle(emcee_dir_output)
        reloaded_emcee = True
    except FileNotFoundError:
        pass

    try:
        starting_point, previous_boundaries, theta_dict = (
            starting_point_load_from_cpickle(optimize_dir_output)
        )
        reloaded_optimize = True
    except FileNotFoundError:
        pass

    print()
    print("reloaded_optimize: ", reloaded_optimize)
    print("reloaded_pyde: ", reloaded_pyde)
    print("reloaded_emcee: ", reloaded_emcee)
    print()
    print("number of multiprocessing threads:", num_threads)
    print(
        "multiprocessing method (should be fork): ", multiprocessing.get_start_method()
    )
    print()

    safe_reload = config_in["parameters"].get("safe_reload", False)

    print("parameters/safe_reload flag (must be True for tinygp): ", safe_reload)
    print()

    if gpu_backend.gpu_enabled():
        # Forking a worker after the parent has used CUDA is unsafe, and JAX
        # already parallelizes on the GPU, so the multiprocessing pool adds
        # hazard with little gain. Force the single-process path (the
        # use_threading_pool guards below are disabled while the GPU is on).
        print("GPU enabled: single-process (use_threading_pool is ignored).")
        print()

    if reloaded_pyde or reloaded_emcee:
        previous_boundaries = mc.bounds

    if reloaded_emcee:
        print("Requested steps:", mc.emcee_parameters["nsteps"])

        mc.emcee_parameters["completed_nsteps"] = int(
            sampler_chain.shape[1] * mc.emcee_parameters["thin"]
        )

        print("Completed:", mc.emcee_parameters["completed_nsteps"])
        pars_input(config_in, mc, input_datasets, reload_emcee=True)
        print("Total:", mc.emcee_parameters["nsteps"])

        """ There's no need to do anything"""
        flatchain = emcee_flatchain(
            sampler_chain, mc.emcee_parameters["nburn"], mc.emcee_parameters["thin"]
        )
        mc.model_setup()
        mc.print_model_info()
        mc.initialize_logchi2()

        results_analysis.print_bayesian_info(mc)

        results_analysis.print_integrated_ACF(
            sampler_chain, theta_dict, mc.emcee_parameters["thin"]
        )

        """ In case the current startin point comes from a previous analysis """
        mc.emcee_dir_output = emcee_dir_output

        if mc.emcee_parameters["nsteps"] <= mc.emcee_parameters["completed_nsteps"]:
            print("Reference Time Tref: ", mc.Tref)
            print()
            print("Dimensions = ", mc.ndim)
            print("Nwalkers = ", mc.emcee_parameters["nwalkers"])
            print("Dimensions = ", mc.ndim)
            print("Steps = ", mc.emcee_parameters["nsteps"])
            print()
            print("Original starting point of emcee:")
            print()

            results_analysis.results_summary(
                mc, starting_point, compute_lnprob=True, is_starting_point=True
            )

            print("emcee results:")
            print()

            results_analysis.results_summary(mc, flatchain)

            print()
            print("emcee completed")
            print()

            if return_output:
                return mc, sampler_chain, sampler_lnprobability
            else:
                return

    if reloaded_emcee:
        sampled = mc.emcee_parameters["completed_nsteps"]
        nsteps_todo = (
            mc.emcee_parameters["nsteps"] - mc.emcee_parameters["completed_nsteps"]
        )

        print("Resuming from a previous run:")
        print("Performed steps = ", mc.emcee_parameters["completed_nsteps"])
        print("Final # of steps = ", mc.emcee_parameters["nsteps"])
        print("Steps to be performed = ", nsteps_todo)
        print()

    else:
        mc = ModelContainerEmcee()
        pars_input(config_in, mc, input_datasets)

        if (
            mc.pyde_parameters["shutdown_jitter"]
            or mc.emcee_parameters["shutdown_jitter"]
        ):
            for dataset_name, dataset in mc.dataset_dict.items():
                dataset.shutdown_jitter()

        # keep track of which version has been used to perform emcee computations
        mc.emcee_parameters["version"] = emcee.__version__[0]

        mc.model_setup()
        mc.print_model_info()
        mc.boundaries_setup()
        mc.initialize_logchi2()

        results_analysis.print_bayesian_info(mc)

        mc.pyde_dir_output = pyde_dir_output
        mc.emcee_dir_output = emcee_dir_output

        mc.emcee_parameters["nwalkers"] = mc.ndim * mc.emcee_parameters["npop_mult"]
        if mc.emcee_parameters["nwalkers"] % 2 == 1:
            mc.emcee_parameters["nwalkers"] += 1

        if not os.path.exists(mc.emcee_dir_output):
            os.makedirs(mc.emcee_dir_output)

        state = None
        sampled = 0
        nsteps_todo = mc.emcee_parameters["nsteps"]

    emcee_input_moves = mc.emcee_parameters.get("moves", "stretch")

    if type(emcee_input_moves) == dict:
        emcee_moves = []
        for key, val in emcee_input_moves.items():
            if key not in moves_definition:
                print(
                    "ERROR: move ",
                    key,
                    " not recognized, available moves are: ",
                    list(moves_definition.keys()),
                )
                quit()
            else:
                emcee_moves.append((moves_definition[key], val))
    else:
        try:
            emcee_moves.append((moves_definition[emcee_input_moves], 1.0))
        except KeyError:
            print(
                "ERROR: move ",
                emcee_input_moves,
                " not recognized, available moves are: ",
                list(moves_definition.keys()),
            )
            quit()

    global log_priors_likelihood

    def log_priors_likelihood(theta):

        # start = time.time()

        log_priors, log_likelihood = mc.log_priors_likelihood(theta)
        # end = time.time()
        # print("--------Computation took {0:.16f} seconds".format(end - start))

        return log_priors + log_likelihood

    global log_priors_likelihood_blob

    def log_priors_likelihood_blob(theta):

        # start = time.time()

        log_priors, log_likelihood = mc.log_priors_likelihood(theta)
        # end = time.time()
        # print("--------Computation took {0:.16f} seconds".format(end - start))

        return log_priors + log_likelihood, log_priors

    """ F-batch (GPU-only, additive): when the container is structurally
        batchable on GPU, route the full (nwalkers, ndim) batch through the
        proven mc.lnprob_batch() seam -- ONE jitted dispatch for the whole
        ensemble instead of one per walker -- via emcee 3 vectorize=True.

        Cardinal rule: strictly additive. The per-theta blob wrapper above
        stays the exact fallback: every container that is not
        GPU-batchable (gpu disabled, priors excluded from lnprob, or a
        structure without a proven batch counterpart) keeps today's
        vectorize=False behaviour, byte-identical. lnprob_batch is never
        allowed to raise: any un-batchable term falls back to the exact
        per-theta evaluations internally (logged with the triggering term).

        Gating: include_priors must be True because lnprob_batch returns
        (priors + likelihood) exactly like the per-theta blob wrapper only
        in that case (its parity proof covers that setting).

        Blob parity (F1 fixed): emcee State iteration (and the 4-way unpack
        in the run_mcmc call below, plus pyde/getresults "sampler_lnprior")
        expects the legacy per-walker (lnprob, log_priors) tuples. The
        recorded blob is now identical in both modes on EVERY row: finite
        rows use the batched log_priors (same closed-form return_priors
        calls, same order -- bit-equal), and every non-finite row resolves
        through the exact per-theta path, so its blob is the per-theta
        log_priors (finite for finite-prior/-inf-likelihood positions,
        -inf for out-of-bounds). Fit, chain, positions, lnprob and
        acceptance are identical in both modes.

        emcee version gate: the 2-D vectorized contract is only validated
        on emcee 3.x; any other emcee major logs one line and falls back
        to the per-theta path (no crash, no dependency change). """
    emcee_vectorize = False
    if gpu_backend.gpu_enabled() and bool(getattr(mc, "include_priors", True)):
        try:
            _batch_ok, _batch_why = mc._lnprob_batch_supported()
        except Exception:
            _batch_ok, _batch_why = False, "gate check failed"
        emcee_vectorize = bool(_batch_ok)
        if emcee_vectorize:
            _gate_ok, _gate_msg = _emcee_vectorize_gate(
                getattr(emcee, "__version__", "")
            )
            if not _gate_ok:
                emcee_vectorize = False
                _batch_why = _gate_msg
        if emcee_vectorize:

            def log_priors_likelihood_blob(theta):
                theta = np.asarray(theta)
                if theta.ndim == 2:
                    # emcee vectorize=True passes the full (nwalkers, ndim)
                    # batch -> per-walker (lnprob, log_priors) in one dispatch
                    try:
                        return _fbatch_2d_results(mc, theta)
                    except Exception:
                        # cardinal rule: exact per-theta semantics, never a raise
                        out = []
                        for row in theta:
                            lp, ll = mc.log_priors_likelihood(row)
                            out.append((lp + ll, lp))
                        return out
                # 1-D (defensive; emcee vectorize=True never sends it here)
                log_priors, log_likelihood = mc.log_priors_likelihood(theta)
                return log_priors + log_likelihood, log_priors
        else:
            print(
                "emcee: F-batch vectorized lnprob not available "
                "({0}) -- using the per-theta path".format(
                    "gpu disabled or priors excluded"
                    if not (
                        gpu_backend.gpu_enabled()
                        and bool(getattr(mc, "include_priors", True))
                    )
                    else _batch_why
                )
            )

    if reloaded_emcee:
        state, sampler = emcee_simpler_load_from_cpickle(emcee_dir_output)

    """
      emcee:
    npop_mult: 4
    nsteps: 20000
    nburn: 7000
    nsave: 5000
    thin: 100
    moves: 
      de: 0.5
      de_snooker: 0.5
    """

    print("Include priors: ", mc.include_priors)
    print()
    print("Reference Time Tref: ", mc.Tref)
    print()
    print("Dimensions = ", mc.ndim)
    print("Nwalkers = ", mc.emcee_parameters["nwalkers"])
    print("Moves = ", emcee_input_moves)
    print()

    if mc.emcee_parameters["version"] == "2":
        mc.emcee_parameters["use_threading_pool"] = False

    try:
        progress_bar = config_in["solver"]["emcee"].get(
            "progress_bar", True
        ) and config_in["parameters"].get("progress_bar", True)
    except:
        progress_bar = True

    if reloaded_emcee:
        sys.stdout.flush()
        pass
    elif reloaded_pyde:
        theta_dict_legacy = theta_dict.copy()
        population_legacy = population.copy()

        theta_dict = results_analysis.get_theta_dictionary(mc)
        population = np.zeros(
            [mc.emcee_parameters["nwalkers"], mc.ndim], dtype=np.double
        )

        for theta_name, theta_i in theta_dict.items():
            population[:, theta_i] = population_legacy[:, theta_dict_legacy[theta_name]]
            mc.bounds[theta_i] = previous_boundaries[theta_dict_legacy[theta_name]]

        starting_point = np.median(population, axis=0)
        # print(starting_point)
        # print(population)

        print("Using previous population as starting point. ")
        print()
        sys.stdout.flush()

    elif mc.starting_point_flag or reloaded_optimize:
        if reloaded_optimize:
            print("Using the output from a previous optimize run as starting point")
            theta_dict_legacy = theta_dict.copy()
            starting_point_legacy = starting_point.copy()
            theta_dict = results_analysis.get_theta_dictionary(mc)
            for theta_name, theta_i in theta_dict.items():
                starting_point[theta_i] = starting_point_legacy[
                    theta_dict_legacy[theta_name]
                ]
        else:
            print("Using user-defined starting point from YAML file")
            mc.starting_points_setup()
            starting_point = mc.starting_point

        population = np.zeros(
            [mc.emcee_parameters["nwalkers"], mc.ndim], dtype=np.double
        )

        for jj, val in enumerate(starting_point):
            if np.isfinite(val):
                if mc.emcee_parameters.get("starts_relative_dispersion", True):
                    if val == 0.00:
                        population[:, jj] = np.random.normal(
                            val, 0.0001, size=mc.emcee_parameters["nwalkers"]
                        )
                    else:
                        population[:, jj] = np.random.normal(
                            val,
                            np.abs(val) * 0.0001,
                            size=mc.emcee_parameters["nwalkers"],
                        )
                else:
                    population[:, jj] = np.random.normal(
                        val, 0.0001, size=mc.emcee_parameters["nwalkers"]
                    )

            else:
                population[:, jj] = np.random.uniform(
                    mc.bounds[jj, 0],
                    mc.bounds[jj, 1],
                    size=mc.emcee_parameters["nwalkers"],
                )

        # for ii in range(0, mc.emcee_parameters['nwalkers']):
        #    population[ii, :] = np.random.normal(starting_point, 0.0000001)
        # for ii in range(0, mc.emcee_parameters['nwalkers']):
        #     print(population[ii, :])

        print("to write a synthetic population extremely close to the starting values.")
        print("Undefined values have a uniform random value withing the boundaries.")
        print("WARNING: Initial state check of emcee will ne skipped")
        print()
        emcee_skip_check = True
        sys.stdout.flush()

    else:
        """ None of the previous cases has been satisfied, we have to run PyDE """
        try:
            from pyde.de import DiffEvol
        except (ModuleNotFoundError, ImportError):
            print(
                "ERROR! PyDE is not installed, run first with optimize instead of emcee"
            )
            quit()

        time_start_pyde = time()

        if not os.path.exists(mc.pyde_dir_output):
            os.makedirs(mc.pyde_dir_output)

        print(
            "Using threading pool for PyDE:", mc.pyde_parameters["use_threading_pool"]
        )

        print("PyDE running")
        sys.stdout.flush()

        if mc.pyde_parameters["use_threading_pool"] and not gpu_backend.gpu_enabled():
            with multiprocessing.Pool(num_threads) as pool:
                de = DiffEvol(
                    log_priors_likelihood,
                    mc.bounds,
                    mc.emcee_parameters["nwalkers"],
                    maximize=True,
                    pool=pool,
                )

                de.optimize(int(mc.pyde_parameters["ngen"]))
        else:
            de = DiffEvol(
                log_priors_likelihood,
                mc.bounds,
                mc.emcee_parameters["nwalkers"],
                maximize=True,
            )

            de.optimize(int(mc.pyde_parameters["ngen"]))

        population = de.population
        starting_point = np.median(population, axis=0)

        theta_dict = results_analysis.get_theta_dictionary(mc)

        """ bounds redefinition and fix for PyDE anomalous results """
        if mc.recenter_bounds_flag:
            pyde_save_to_pickle(
                mc, population, starting_point, theta_dict, prefix="orig"
            )

            mc.recenter_bounds(starting_point)
            population = mc.fix_population(starting_point, population)
            starting_point = np.median(population, axis=0)

            print("Boundaries redefined after PyDE output")

        pyde_save_to_pickle(mc, population, starting_point, theta_dict)

        time_end_pyde = time()
        print(
            "PyDE completed, it took {0:12.1f} seconds".format(
                time_end_pyde - time_start_pyde
            )
        )
        print()
        sys.stdout.flush()

        if safe_reload:
            print(" safe_reload flag on True, the program will now quit ")
            print(" You have to relaunch it again in order for emcee to work properly")
            print(" No worries, your PyDE results have been saved!")
            quit()

    if reloaded_emcee:
        print("Original starting point of emcee:")
        print()

    results_analysis.results_summary(
        mc, starting_point, compute_lnprob=not (safe_reload), is_starting_point=True
    )
    sys.stdout.flush()

    print()
    print("*************************************************************")
    print()
    print("Running emcee")
    print()
    print("emcee version: ", emcee.__version__)
    if mc.emcee_parameters["version"] == "2":
        print("WARNING: upgrading to version 3 is strongly advised")
    print(
        "Using threading pool for emcee:",
        mc.emcee_parameters.get("use_threading_pool", False),
    )
    print()

    if reloaded_emcee:
        print("Using reloaded sampler")
        print()
    else:
        emcee_vec_kwargs = {"vectorize": True} if emcee_vectorize else {}
        if emcee_vectorize:
            print(
                "emcee: vectorize=True -- F-batch lnprob (whole ensemble, "
                "one dispatch per step)"
            )
            print()
        sampler = emcee.EnsembleSampler(
            mc.emcee_parameters["nwalkers"],
            mc.ndim,
            log_priors_likelihood_blob,
            moves=emcee_moves,
            **emcee_vec_kwargs,
        )

    if mc.emcee_parameters["nsave"] > 0:
        print("Saving temporary steps")
        print()
        niter = int(np.ceil(nsteps_todo / mc.emcee_parameters["nsave"]))

        for i in range(0, niter):
            if (
                mc.emcee_parameters["use_threading_pool"]
                and not gpu_backend.gpu_enabled()
            ):
                with multiprocessing.Pool(num_threads) as pool:
                    sampler.pool = pool
                    population, prob, state, prior = sampler.run_mcmc(
                        population,
                        int(mc.emcee_parameters["nsave"]),
                        thin=mc.emcee_parameters["thin"],
                        rstate0=state,
                        progress=progress_bar,
                        skip_initial_state_check=emcee_skip_check,
                    )
            else:
                population, prob, state, prior = sampler.run_mcmc(
                    population,
                    int(mc.emcee_parameters["nsave"]),
                    thin=mc.emcee_parameters["thin"],
                    rstate0=state,
                    progress=progress_bar,
                    skip_initial_state_check=emcee_skip_check,
                )

            sampled += mc.emcee_parameters["nsave"]
            theta_dict = results_analysis.get_theta_dictionary(mc)
            emcee_save_to_cpickle(
                mc,
                starting_point,
                population,
                prob,
                state,
                sampler,
                theta_dict,
                samples=sampled,
            )

            flatchain = emcee_flatchain(
                sampler.chain, mc.emcee_parameters["nburn"], mc.emcee_parameters["thin"]
            )

            results_analysis.print_integrated_ACF(
                sampler.chain, theta_dict, mc.emcee_parameters["thin"]
            )

            results_analysis.results_summary(mc, flatchain)

            print()
            print(sampled, "  steps completed, average lnprob:, ", np.median(prob))
            print()
            sys.stdout.flush()

            if safe_reload:
                print(" safe_reload flag on True, the program will now quit ")
                print(
                    " You have to relaunch it again in order for emcee to work properly"
                )
                print(" No worries, your PyDE results have been saved!")
                quit()

            # It turns out that reloading the sampler from the file will
            # result in faster parallelization...
            state, sampler = emcee_simpler_load_from_cpickle(emcee_dir_output)

    else:
        if mc.emcee_parameters["use_threading_pool"] and not gpu_backend.gpu_enabled():
            with multiprocessing.Pool(num_threads) as pool:
                sampler.pool = pool
                population, prob, state, prior = sampler.run_mcmc(
                    population,
                    nsteps_todo,
                    thin=mc.emcee_parameters["thin"],
                    rstate0=state,
                    progress=progress_bar,
                    skip_initial_state_check=emcee_skip_check,
                )

        else:
            print("Warning: NOT using threading, pool, performances will be slower")
            population, prob, state, prior = sampler.run_mcmc(
                population,
                nsteps_todo,
                thin=mc.emcee_parameters["thin"],
                rstate0=state,
                progress=progress_bar,
                skip_initial_state_check=emcee_skip_check,
            )

        sampled += nsteps_todo

        theta_dict = results_analysis.get_theta_dictionary(mc)
        emcee_save_to_cpickle(
            mc,
            starting_point,
            population,
            prob,
            state,
            sampler,
            theta_dict,
            samples=sampled,
        )

        flatchain = emcee_flatchain(
            sampler.chain, mc.emcee_parameters["nburn"], mc.emcee_parameters["thin"]
        )

        results_analysis.print_integrated_ACF(
            sampler.chain, theta_dict, mc.emcee_parameters["thin"]
        )

        results_analysis.results_summary(mc, flatchain)

    print()
    print("emcee completed")
    print()

    """ A dummy file is written to let the cpulimit script to proceed with the next step
        Nopt sure if this step is needed anymore
    """
    emcee_write_dummy_file(mc)

    if return_output:
        return mc, sampler.chain, sampler.lnprobability
