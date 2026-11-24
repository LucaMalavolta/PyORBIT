#!/usr/bin/env python
"""F-batch emcee A/B: the SAME emcee run shape (vectorize=True, identical
seed/initial population/move) with the log_prob_fn routed two ways, so the
A/B isolates exactly the jitted-dispatch collapse:

  BEFORE  per-theta : log_prob(theta2d) -> [ per-theta lnprob(theta2d[i])
                                             for i in range(nwalkers) ]
                        => nwalkers dispatches per emcee step.
  AFTER   batched   : log_prob(theta2d) -> mc.lnprob_batch(theta2d)
                        => ONE dispatch per step (or lnprob_batch's exact
                           internal per-theta fallback, which we DETECT).

HARD GATE (the batch must never change a single number):
  * parity on the same theta2d: AFTER vs BEFORE maxdiff < 1e-6 (expect 0.0)
    and the -inf (out-of-bounds) row sets identical;
  * both legs (vec=True) must produce the IDENTICAL final chain
    (same seed, same accept/reject stream) -> final-population maxdiff 0.0.

Cardinal rule evidence: AFTER's lnprob values are the production
per-theta values (lnprob_batch proven bit-exact; any un-batchable term
falls back internally and we REPORT it, we do not force it).

Run:
  source ~/anaconda3/etc/profile.d/conda.sh && conda activate pyorbit-gpu
  python development/emcee_v2_ab.py
Output: /tmp/opencode/pyorbit_emcee_v2/ (results.json, stdout.log, progress.log)
"""

import json
import os
import statistics
import sys
import time

OUT = "/tmp/opencode/pyorbit_emcee_v2"
N_POINTS = 600
N_STEPS = 60
WARMUP_STEPS = 5
TIMED_REPS = 2
SEED = 42

os.makedirs(OUT, exist_ok=True)
LOG = open(os.path.join(OUT, "stdout.log"), "w")


def say(*args):
    msg = " ".join(str(a) for a in args)
    print(msg)
    LOG.write(msg + "\n")
    LOG.flush()


def build_mc(n_points):
    """Same synthetic RV container as development/gpu_perf_surfaces.py
    _build_mc_leg (identical data generator, ordering, single-planet
    Kepler RV) with gpu: true."""
    import numpy as np
    import yaml

    import pyorbit.subroutines.input_parser as ip
    import pyorbit.common.gpu_backend as gpu_backend
    from pyorbit.classes.model_container_emcee import ModelContainerEmcee

    workdir = os.path.join(OUT, "mc_leg_gpu")
    os.makedirs(workdir, exist_ok=True)

    tref = 2457000.0
    rng = np.random.default_rng(1234)
    epoch = tref - 150.0 + np.sort(rng.uniform(0.0, 300.0, n_points))
    P_true, K_true = 21.0, 18.0
    phase = 2.0 * np.pi * ((epoch - tref) / P_true)
    value = K_true * np.cos(phase) + 0.5 * rng.normal(size=n_points)
    error = np.full(n_points, 1.2)
    data_path = os.path.join(workdir, "synthetic_rv.dat")
    with open(data_path, "w") as fh:
        for tv, vv, ee in zip(epoch, value, error):
            fh.write("{0:.6f}\t{1:.5f}\t{2:.3f}\t0\n".format(tv, vv, ee))

    conf = {
        "output": os.path.join(workdir, "out"),
        "inputs": {
            "synthetic_rv": {
                "file": data_path,
                "kind": "RV",
                "models": ["radial_velocities"],
            }
        },
        "common": {
            "planets": {
                "p": {
                    "orbit": "keplerian",
                    "boundaries": {
                        "P": [15.0, 30.0],
                        "K": [0.1, 100.0],
                        "e": [0.0, 0.6],
                    },
                }
            },
            "star": {"star_parameters": {"priors": {}}},
        },
        "models": {"radial_velocities": {"planets": ["p"]}},
        "parameters": {"Tref": tref, "gpu": True},
        "solver": {"emcee": {"npop_mult": 4, "nsteps": 10, "nburn": 5, "thin": 1}},
    }
    conf_path = os.path.join(workdir, "config.yaml")
    with open(conf_path, "w") as fh:
        yaml.safe_dump(conf, fh)

    config_in = ip.yaml_parser(conf_path)
    gpu_backend.initialize(config_in["parameters"])

    mc = ModelContainerEmcee()
    ip.pars_input(config_in, mc, input_datasets=None)
    mc.model_setup()
    mc.boundaries_setup()
    mc.initialize_logchi2()
    return mc


def parity_check(mc, nwalkers):
    """HARD GATE: same theta2d -> BEFORE (per-theta) vs AFTER (lnprob_batch)."""
    import numpy as np

    # deterministic theta set: in-bounds perturbations of the midpoint +
    # out-of-bounds rows (one per violated bound), shaped (nwalkers, ndim)
    theta2d = np.zeros((nwalkers, mc.ndim))
    mid = 0.5 * (mc.bounds[:, 0] + mc.bounds[:, 1])
    rng = np.random.default_rng(7)
    for i in range(nwalkers):
        row = mid * (1.0 + 0.05 * rng.normal(size=mc.ndim))
        theta2d[i] = np.clip(row, mc.bounds[:, 0] + 0.1, mc.bounds[:, 1] - 0.1)
    oob = [mc.bounds[:, 0] - 1.0, mc.bounds[:, 1] + 1.0]
    theta2d[nwalkers - 2] = oob[0]
    theta2d[nwalkers - 1] = oob[1]

    before_lp, before_blob = [], []
    for i in range(nwalkers):
        lp, ll = mc.log_priors_likelihood(theta2d[i])
        before_lp.append(float(np.asarray(lp) + np.asarray(ll)))
        before_blob.append(float(lp))
    before = np.array(before_lp)
    before_blob = np.array(before_blob)
    from pyorbit.samplers.pyorbit_emcee import _fbatch_2d_results

    res_pairs = _fbatch_2d_results(mc, theta2d)
    after = np.array([p[0] for p in res_pairs])
    # F1 semantics: non-finite rows carry the exact per-theta log_priors
    # (the production 2-D wrapper; no blind -inf masking any more)
    after_blob = np.array([p[1] for p in res_pairs])

    say("### PARITY (same theta2d, BEFORE per-theta vs AFTER lnprob_batch)")
    say("    BEFORE lnprob: " + " ".join("{0:.6g}".format(v) for v in before))
    say("    AFTER  lnprob: " + " ".join("{0:.6g}".format(v) for v in after))

    nb = np.isneginf(before)
    na = np.isneginf(after)
    inf_match = bool(np.array_equal(nb, na))
    finite = ~(nb | na)
    maxdiff = (
        float(np.max(np.abs(before[finite] - after[finite]))) if finite.any() else 0.0
    )

    bbad = ~np.isfinite(before_blob)
    aad = ~np.isfinite(after_blob)
    blob_match = bool((bbad == aad).all())
    bfinite = ~(bbad | aad)
    bmaxdiff = (
        float(np.max(np.abs(before_blob[bfinite] - after_blob[bfinite])))
        if bfinite.any()
        else 0.0
    )
    say(
        "    rows={0} -inf(B)={1} -inf(A)={2} -inf-set-match={3} lnprob maxdiff={4}".format(
            nwalkers, int(nb.sum()), int(na.sum()), inf_match, maxdiff
        )
    )
    say(
        "    blob(log_priors) -inf-set-match={0} maxdiff={1}".format(
            blob_match, bmaxdiff
        )
    )
    say(
        "    NOTE: rows {0} = OOB below-bounds, row {1} = OOB above-bounds".format(
            nwalkers - 2, nwalkers - 1
        )
    )

    ok = (
        (maxdiff < 1e-6)
        and (bmaxdiff < 1e-6)
        and inf_match
        and blob_match
        and (not np.isnan(after).any())
    )
    say("    PARITY GATE: " + ("PASS" if ok else "FAIL"))
    return {
        "maxdiff": maxdiff,
        "inf_set_match": inf_match,
        "blob_match": blob_match,
        "blob_maxdiff": bmaxdiff,
        "pass": bool(ok),
    }


def run_leg(tag, mc, use_batch, nwalkers):
    import numpy as np
    import emcee

    n0 = {"n_calls": 0}

    if use_batch:
        from pyorbit.samplers.pyorbit_emcee import _fbatch_2d_results

        def log_prob(theta):
            theta = np.asarray(theta)
            if theta.ndim == 2:
                n0["n_calls"] += 1
                # the production vectorized wrapper (pyorbit_emcee),
                # F1 semantics: non-finite rows carry the exact per-theta
                # log_priors (no blind -inf masking any more)
                return _fbatch_2d_results(mc, theta)
            lp, ll = mc.log_priors_likelihood(theta)
            return lp + ll, lp
    else:

        def log_prob(theta):
            theta = np.asarray(theta)
            if theta.ndim == 2:
                n0["n_calls"] += 1
                # the production per-theta wrapper, per walker
                out = []
                for i in range(theta.shape[0]):
                    lp, ll = mc.log_priors_likelihood(theta[i])
                    out.append((lp + ll, lp))
                return out
            lp, ll = mc.log_priors_likelihood(theta)
            return lp + ll, lp

    theta0 = np.array(
        [0.5 * (mc.bounds[i, 0] + mc.bounds[i, 1]) for i in range(mc.ndim)]
    )
    rng = np.random.default_rng(SEED)
    population = np.zeros((nwalkers, mc.ndim))
    for jj, val in enumerate(theta0):
        sigma = 1e-4 if val == 0.0 else abs(val) * 1e-4
        population[:, jj] = rng.normal(val, sigma, size=nwalkers)

    # emcee seeds its internal RandomState from the GLOBAL numpy state at
    # construction -> identical proposal stream for both legs (fair A/B)
    np.random.seed(SEED)
    sampler = emcee.EnsembleSampler(
        nwalkers,
        mc.ndim,
        log_prob,
        moves=[(emcee.moves.StretchMove(), 1.0)],
        vectorize=True,
    )

    # warm-up (untimed): absorbs first-use overhead of this leg's path
    population, prob, state, blobs = sampler.run_mcmc(
        population, WARMUP_STEPS, progress=False, skip_initial_state_check=True
    )

    walls = []
    for _ in range(TIMED_REPS):
        a = time.perf_counter()
        population, prob, state, blobs = sampler.run_mcmc(
            population, N_STEPS, progress=False, skip_initial_state_check=True
        )
        walls.append(time.perf_counter() - a)
    sps = [N_STEPS / w for w in walls]
    med_sp = float(statistics.median(sps))

    it = np.asarray(sampler.backend.iteration).ravel()
    acc = np.asarray(sampler.backend.accepted).ravel()
    acc_fr = float(np.sum(acc) / float(np.sum(it)))
    say(
        "### LEG {0} (vec=True, {1})".format(
            tag,
            "lnprob_batch [1 dispatch/step]"
            if use_batch
            else "per-theta loop [{0} dispatches/step]".format(nwalkers),
        )
    )
    say(
        "    median_steps_per_s={0:.3f}  walls_s=[{1}]  reps={2}".format(
            med_sp, ", ".join("{0:.2f}".format(w) for w in walls), TIMED_REPS
        )
    )
    say("    acceptance={0:.4f} (steps_run={1})".format(acc_fr, int(np.sum(it))))
    say(
        "    log_prob batches seen by emcee={0} (1 per step expected)".format(
            n0["n_calls"]
        )
    )
    say("    accept_per_walker=" + " ".join("{0:.3f}".format(v) for v in acc / it))
    return {
        "median_steps_per_s": med_sp,
        "walls_s": walls,
        "acceptance": acc_fr,
        "n_calls": n0["n_calls"],
        "final_population": population,
        "accept_per_walker": acc / it,
    }


def main():
    import numpy as np

    results = {"n_points": N_POINTS, "nsteps": N_STEPS, "seed": SEED}

    say(
        "### emcee_v2 A/B: N={0} pts, steps={1}x{2}+warmup{3}, seed={4}, "
        "both legs vectorize=True".format(
            N_POINTS, TIMED_REPS, N_STEPS, WARMUP_STEPS, SEED
        )
    )

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    mc = build_mc(N_POINTS)
    ndim = mc.ndim
    nwalkers = max(2 * (ndim + 1), 4 * ndim)
    results["ndim"] = ndim
    results["nwalkers"] = nwalkers
    say("    ndim={0} nwalkers={1}".format(ndim, nwalkers))

    try:
        ok, why = mc._lnprob_batch_supported()
    except Exception as exc:
        ok, why = False, "gate raised: {0}".format(type(exc).__name__)
    results["batch_supported"] = {"ok": bool(ok), "why": why}
    say(
        "    container gate _lnprob_batch_supported() -> ok={0} why={1}".format(ok, why)
    )

    import jax

    say("    jax backend: {0}".format(jax.default_backend()))

    parity = parity_check(mc, nwalkers)
    results["parity"] = parity
    if not parity["pass"]:
        say("### HARD PARITY GATE FAILED -- STOPPING (do not report a speedup)")
        results["status"] = "parity_gate_failed"
        _finish(results)
        sys.exit(1)

    leg_before = run_leg("BEFORE", mc, use_batch=False, nwalkers=nwalkers)
    leg_after = run_leg("AFTER", mc, use_batch=True, nwalkers=nwalkers)

    diff_pop = float(
        np.max(np.abs(leg_before["final_population"] - leg_after["final_population"]))
    )
    say(
        "### CHAIN PARITY (both legs vec=True, identical seed/stream): "
        "maxdiff(final populations) = {0} (expect 0.0)".format(diff_pop)
    )
    say(
        "    acc BEFORE={0:.6f}  acc AFTER={1:.6f}  (lnprob values identical -> "
        "decisions identical)".format(leg_before["acceptance"], leg_after["acceptance"])
    )

    speedup = leg_after["median_steps_per_s"] / leg_before["median_steps_per_s"]
    say("### SPEEDUP AFTER/BEFORE = {0:.4g}x steps/s".format(speedup))

    fallback_lines = []
    capture = sys.argv[1] if len(sys.argv) > 1 else os.path.join(OUT, "stdout.log")
    try:
        with open(capture) as fh:
            for line in fh:
                if "lnprob_batch: per-theta fallback" in line:
                    fallback_lines.append(line.strip())
    except FileNotFoundError:
        pass
    results["fallback_lines"] = fallback_lines
    say("### lnprob_batch internal fallback triggered {0}x".format(len(fallback_lines)))
    for line in fallback_lines[:5]:
        say("    " + line)

    results["leg_before"] = {
        k: v for k, v in leg_before.items() if k != "final_population"
    }
    results["leg_after"] = {
        k: v for k, v in leg_after.items() if k != "final_population"
    }
    results["chain_maxdiff"] = diff_pop
    results["speedup"] = speedup
    results["status"] = "pass" if diff_pop == 0.0 else "chains_differ"
    _finish(results)


def _finish(results):
    with open(os.path.join(OUT, "results.json"), "w") as fh:
        json.dump(results, fh, indent=2, default=str)
    say("### DONE -> {0}/results.json".format(OUT))


if __name__ == "__main__":
    main()
