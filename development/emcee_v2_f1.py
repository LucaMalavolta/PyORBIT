#!/usr/bin/env python
"""F-batch F1 test: the recorded per-theta blob must be IDENTICAL in both
(per-theta and vectorized) modes on EVERY row -- the FINITE log_priors for
the finite-prior/-inf-likelihood class of rows included (the F1 fix: the old
``blobs[bad] = -inf`` masking in the 2-D wrapper diverged from per-theta
mode on exactly that class).

Also unit-checks the two in-scope polishes:
  * the emcee major-version guard (FIX 2): major != 3 -> gates OFF + one log
    line, no crash, no dependency change;
  * the lnprob_batch per-theta fallback log throttle (Obs d): first
    occurrence logged, repeats suppressed (no stdout flood).

HARD GATES (this script exits 1 if any fail):
  A. SEAM UNIT: a faithful seam stub producing the F1 row class -> the
     production 2-D function ``_fbatch_2d_results`` must return, per row,
     exactly the per-theta (lnprob, log_priors):
        row A (finite):        ( -3.0,  -1.0 )
        row B (F1 class):      ( -inf,  -3.0 )  <- blob FINITE, not -inf
        row C (out-of-bounds): ( -inf,  -inf )
     and the result must be consumable by emcee 3 (State iteration
     contract: list of (log_prob, *blob) tuples).
     Before the F1 fix this SAME test FAILS on row B's blob (masked -inf).
  B. REAL CONTAINER: 6-param kepler RV (gpu: true); the theta batch
     (finite-e0.13 / below-bounds / above-bounds / circular-e0
     / parabolic-e1.0 / hyperbolic-e1.05) must be lnprob- AND blob-
     identical per-theta on every row, plus a zero-error (sigma -> 0)
     probe; prints the honest "### RESULT F1-class ... : Y/N" for
     finite-prior/-inf-likelihood rows (expected: N -- the Part A
     seam unit remains the authoritative F1 coverage).
  C. GUARD: _emcee_vectorize_gate('4.0') -> off + message; '3.1.6' -> on;
     malformed -> off. No exception, no import of another emcee.
  D. THROTTLE: 150 forced per-theta fallback calls -> exactly 2 log lines
     (occurrence 1 and 100), not 150.

Run:
  source ~/anaconda3/etc/profile.d/conda.sh && conda activate pyorbit-gpu
  python development/emcee_v2_f1.py
Output: /tmp/opencode/pyorbit_f1_fix/ (stdout.log, results.json)
"""

import io
import json
import os
import sys

OUT = "/tmp/opencode/pyorbit_f1_fix"
os.makedirs(OUT, exist_ok=True)
LOG = open(os.path.join(OUT, "stdout.log"), "a")


def say(*args):
    msg = " ".join(str(a) for a in args)
    print(msg)
    LOG.write(msg + "\n")
    LOG.flush()


# --------------------------------------------------------------------------
# A. SEAM UNIT -- the production function against a faithful seam stub
# --------------------------------------------------------------------------
def part_a_seam_unit():
    import numpy as np

    from pyorbit.samplers.pyorbit_emcee import _fbatch_2d_results

    class FakeBatchable(object):
        """Faithful seam stub: same API and same per-row semantics as the
        production container's lnprob_batch / _lnprob_batch_priors /
        log_priors_likelihood, with a deterministic row in EACH class:

          row A x=0.5  finite:                 (lp, ll) = (-1.0, -2.0)
          row B x=1.5  F1 CLASS:                (lp, ll) = (-3.0, -inf)
                       finite priors, -inf likelihood
          row C x=11.0 out-of-bounds:          (lp, ll) = (-inf, -inf)
        """

        include_priors = True

        def check_bounds(self, theta):
            return bool(theta[0] < 10.0)

        def log_priors_likelihood(self, theta):
            x = float(theta[0])
            if x < 10.0:
                if x < 1.0:
                    return -1.0, -2.0
                return -3.0, -np.inf  # F1 class: finite lp, -inf ll
            return -np.inf, -np.inf

        def lnprob_batch(self, theta_list):
            rows = [float(t[0]) for t in theta_list]
            out = np.empty(len(rows), dtype=np.double)
            for i, x in enumerate(rows):
                lp, ll = self.log_priors_likelihood(np.array([x]))
                out[i] = lp + ll
            return out

        def _lnprob_batch_priors(self, theta_list):
            rows = [float(t[0]) for t in theta_list]
            out = np.empty(len(rows), dtype=np.double)
            for i, x in enumerate(rows):
                lp, _ll = self.log_priors_likelihood(np.array([x]))
                out[i] = lp
            return out

        def _lnprob_batch_supported(self):
            return True, None

    mc = FakeBatchable()
    theta2d = np.array([[0.5], [1.5], [11.0]])

    # per-theta ground truth (the production per-theta wrapper semantics)
    per_theta = []
    for row in theta2d:
        lp, ll = mc.log_priors_likelihood(row)
        per_theta.append((float(lp) + float(ll), float(lp)))

    res = _fbatch_2d_results(mc, theta2d)
    say("### A. SEAM UNIT (production _fbatch_2d_results vs per-theta)")
    for i, (got, want) in enumerate(zip(res, per_theta)):
        say(
            "    row {0}: wrapper=({1!r}, {2!r})  per-theta=({3!r}, {4!r}){5}".format(
                "ABC"[i],
                got[0],
                got[1],
                want[0],
                want[1],
                "" if got == want else "   << MISMATCH",
            )
        )

    lnprobs = np.array([r[0] for r in res])
    blobs = np.array([r[1] for r in res])
    want_lp = np.array([r[0] for r in per_theta])
    want_blob = np.array([r[1] for r in per_theta])

    f1_row = 1  # the finite-prior / -inf-likelihood row
    batched_priors = np.asarray(mc._lnprob_batch_priors(theta2d))

    # (a) lnprob identical in both modes (both -inf on row B)
    lnprob_match = bool(
        np.array_equal(lnprobs, want_lp)
        and np.isneginf(lnprobs[f1_row])
        and np.isneginf(want_lp[f1_row])
    )
    # (b) blob identical in both modes; row B's blob is the FINITE log_priors
    blob_match = bool(np.array_equal(blobs, want_blob) and np.isfinite(blobs[f1_row]))
    # the stub genuinely exercises the divergence: the batched prior for row
    # B IS finite -- the OLD masking code set the blob to -inf here, the
    # F1 fix keeps the per-theta finite value
    divergence_exercised = bool(np.isfinite(batched_priors[f1_row]))

    say("    (a) lnprob identical (both -inf on F1 row): " + str(lnprob_match))
    say(
        "    (b) blob identical; F1-row blob finite log_priors={0!r}".format(
            blobs[f1_row]
        )
    )
    say(
        "    F1 divergence exercised (batched prior finite, old code masked -inf): "
        + str(divergence_exercised)
    )

    # emcee-3 consumption contract: emcee's vectorize path, run exactly as it
    # does at EVERY sampled step, must accept the list of (log_prob, *blob)
    # tuples and preserve the per-theta blobs -- including a FINITE blob on a
    # row whose log_prob is -inf (exactly the F1 class). We drive the
    # deterministic ``compute_log_prob`` path (no StretchMove position drift),
    # i.e.   results = log_prob_fn(p2d)
    #        blob    = [l[1:] for l in results if len(l) > 1]  (len-1 dims
    #        squeezed, -inf log_prob rows kept)
    emcee_ok = False
    emcee_note = "emcee import failed"
    try:
        import emcee

        def log_prob(theta):
            theta = np.asarray(theta)
            if theta.ndim == 2:
                return _fbatch_2d_results(mc, theta)
            lp, ll = mc.log_priors_likelihood(theta)
            return lp + ll, lp

        sampler = emcee.EnsembleSampler(
            len(theta2d),
            1,
            log_prob,
            moves=[(emcee.moves.StretchMove(), 1.0)],
            vectorize=True,
        )
        got_lp, got_blob = sampler.compute_log_prob(theta2d)
        got_lp = np.asarray(got_lp, dtype=np.double).reshape(-1)
        got_blob = np.asarray(got_blob, dtype=np.double).reshape(-1)
        emcee_ok = bool(
            got_lp.shape == want_lp.shape
            and bool(np.array_equal(got_lp, want_lp))
            and got_blob.shape == want_blob.shape
            and bool(np.array_equal(got_blob, want_blob))
            and bool(np.isfinite(got_blob[f1_row]))
        )
    except Exception as exc:  # noqa: BLE001
        emcee_note = "{0}: {1}".format(type(exc).__name__, exc)

    pass_a = bool(lnprob_match and blob_match and divergence_exercised and emcee_ok)
    if not emcee_ok:
        say("    emcee-3 consumption: {0}".format(emcee_note))
    else:
        say(
            "    emcee-3 consumption: blobs recorded == per-theta blobs incl. finite F1-row blob: OK"
        )
    say("    A VERDICT: " + ("PASS" if pass_a else "FAIL"))
    return {
        "lnprob_match": bool(lnprob_match),
        "blob_match": bool(blob_match),
        "divergence_exercised": bool(divergence_exercised),
        "emcee_consumed": emcee_ok,
        "emcee_note": emcee_note if not emcee_ok else "ok",
        "pass": pass_a,
    }


def main():
    import numpy as np

    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if repo not in sys.path:
        sys.path.insert(0, repo)

    results = {}
    say("emcee_v2 F1  (wrapper blob parity + guard + throttle; OUT={0})".format(OUT))

    a = part_a_seam_unit()
    results["a_seam_unit"] = a
    if not a["pass"]:
        _finish(results, "a_seam_unit_fail")

    # ----------------------------------------------------------------------
    # B. REAL CONTAINER -- 6-param kepler RV, wrapper vs per-theta on every row
    # ----------------------------------------------------------------------
    import yaml

    import pyorbit.common.gpu_backend as gpu_backend
    import pyorbit.subroutines.input_parser as ip
    import pyorbit.subroutines.results_analysis as ra
    from pyorbit.classes.model_container_emcee import ModelContainerEmcee
    from pyorbit.samplers.pyorbit_emcee import _fbatch_2d_results

    tref = 2457000.0

    def build_rv_mc(subdir, n_pts, n_zero_sigma=0, init_gpu=False):
        """Tiny synthetic kepler RV container for the probes below.

        ``n_zero_sigma`` data rows at ~30% of the epochs carry sigma = 0.0
        (the sigma -> 0 candidate for a finite-prior / -inf-likelihood row);
        ``init_gpu`` runs the one-way GPU backend init (called exactly once,
        before the first container is built).
        """
        workdir = os.path.join(OUT, subdir)
        os.makedirs(workdir, exist_ok=True)
        rng = np.random.default_rng(1234)
        epoch = tref - 150.0 + np.sort(rng.uniform(0.0, 300.0, n_pts))
        P_true, K_true = 21.0, 18.0
        phase = 2.0 * np.pi * ((epoch - tref) / P_true)
        value = K_true * np.cos(phase) + 0.5 * rng.normal(size=n_pts)
        error = np.full(n_pts, 1.2)
        if n_zero_sigma:
            i0 = int(round(0.3 * n_pts))
            error[i0 : i0 + n_zero_sigma] = 0.0
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
                            "e": [0.0, 1.1],  # wide: e in (1.0, 1.1) stays in-bounds
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
        if init_gpu:
            gpu_backend.initialize(config_in["parameters"])
        mc = ModelContainerEmcee()
        ip.pars_input(config_in, mc, input_datasets=None)
        mc.model_setup()
        mc.boundaries_setup()
        mc.initialize_logchi2()
        return mc

    def probe_batch(mc, tag, theta2d, labels):
        """F1 contract on a REAL container: the wrapper's per-row
        (lnprob, blob) must equal the per-theta (lnprob, log_priors) on
        every row. Returns (parity_ok, rows); rows are json-serializable."""
        per_theta = []
        for row in theta2d:
            lp, ll = mc.log_priors_likelihood(row)
            per_theta.append((float(lp) + float(ll), float(lp)))
        res = _fbatch_2d_results(mc, theta2d)
        rows = []
        all_same = True
        for lab, row, got, want in zip(labels, theta2d, res, per_theta):
            same = bool(got[0] == want[0] and got[1] == want[1])
            all_same = bool(all_same and same)
            f1 = bool(not np.isfinite(want[0]) and np.isfinite(want[1]))
            rows.append(
                {
                    "label": lab,
                    "in_bounds": bool(mc.check_bounds(row)),
                    "wrapper_lnprob": str(got[0]),
                    "wrapper_blob": str(got[1]),
                    "per_theta_lnprob": str(want[0]),
                    "per_theta_log_priors": str(want[1]),
                    "same": same,
                    "f1_class": f1,
                }
            )
            say(
                "    {0} [{1}] in-bounds={2}: wrapper=({3!r}, {4!r})"
                "  per-theta=({5!r}, {6!r}) same={7}{8}".format(
                    tag,
                    lab,
                    rows[-1]["in_bounds"],
                    got[0],
                    got[1],
                    want[0],
                    want[1],
                    same,
                    "" if same else "   << MISMATCH",
                )
            )
        n_f1 = sum(1 for r in rows if r["f1_class"])
        say(
            "    {0} parity on all {1} rows: {2}; "
            "F1-class rows (finite log_priors, -inf lnprob): {3}".format(
                tag, len(rows), all_same, n_f1
            )
        )
        return bool(all_same), rows

    mc = build_rv_mc("mc_real", 600, init_gpu=True)
    td = ra.get_theta_dictionary(mc)
    say("### B. REAL CONTAINER ndim={0} params={1}".format(mc.ndim, list(td)))
    ok, why = mc._lnprob_batch_supported()
    say("    batch gate: ok={0} why={1}".format(ok, why))
    say(
        "    effective e bounds (conf requested [0.0, 1.1]): {0}".format(
            mc.common_models["p"].bounds.get("e")
        )
    )

    mid = 0.5 * (mc.bounds[:, 0] + mc.bounds[:, 1])
    theta_below = mid.copy()
    theta_below[:] = mc.bounds[:, 0] - 1.0
    theta_above = mid.copy()
    theta_above[:] = mc.bounds[:, 1] + 1.0
    theta_finite = mid.copy()
    theta_finite[td["p_sre_coso"]] = 0.3  # e = 0.3^2 + 0.2^2 = 0.13: a true ellipse
    theta_finite[td["p_sre_sino"]] = 0.2
    theta_circ = mid.copy()
    theta_circ[td["p_sre_coso"]] = 0.0  # e = 0: CIRCULAR -- the old row was
    theta_circ[td["p_sre_sino"]] = 0.0  # mislabelled "degenerate e=1"
    theta_par = mid.copy()
    theta_par[td["p_sre_coso"]] = np.sqrt(0.5)  # e = 1.0: parabolic edge
    theta_par[td["p_sre_sino"]] = np.sqrt(0.5)
    theta_hyp = mid.copy()
    theta_hyp[td["p_sre_coso"]] = np.sqrt(1.05 / 2.0)  # e = 1.05: hyperbolic
    theta_hyp[td["p_sre_sino"]] = np.sqrt(1.05 / 2.0)
    theta2d = np.array(
        [theta_finite, theta_below, theta_above, theta_circ, theta_par, theta_hyp],
        dtype=np.double,
    )
    labels = [
        "finite-e0.13",
        "below-bounds",
        "above-bounds",
        "circular-e0",
        "parabolic-e1.0",
        "hyperbolic-e1.05",
    ]
    pass_b1, rows_b = probe_batch(mc, "B", theta2d, labels)

    say("### B2. ZERO-ERROR ROW PROBE (sigma -> 0 candidate)")
    mc0 = build_rv_mc("mc_sigma0", 300, n_zero_sigma=2)
    mid2 = 0.5 * (mc0.bounds[:, 0] + mc0.bounds[:, 1])
    pass_b2, rows_0 = probe_batch(mc0, "B2", mid2.reshape(1, -1), ["in-bounds-mid"])

    f1_labels = [r["label"] for r in rows_b + rows_0 if r["f1_class"]]
    constructible = bool(f1_labels)
    say(
        "### RESULT F1-class (finite log_priors / -inf lnprob) constructible"
        " in-bounds on the real kepler container: {0}".format(
            "Y" if constructible else "N"
        )
    )
    if constructible:
        say(
            "    constructible on: {0} -- blob carried the exact per-row"
            " finite log_priors; parity PASS on all rows".format(f1_labels)
        )
    else:
        say(
            "    honest finding (N) -- the F1 class is in-bounds-unconstructible here because:"
        )
        say(
            "      - the per-theta container NaN guard maps any NaN prior/likelihood to (-inf, -inf);"
        )
        say(
            "      - e in {1.0, 1.05} (parabolic/hyperbolic) and the out-of-bounds rows land at (-inf, -inf) or finite -- never (finite log_priors, -inf lnprob);"
        )
        say(
            "      - sigma -> 0 only yields a finite contribution, env = 1/(sigma^2 + jitter^2) with jitter > 0 (B2: lnprob finite) -- not an F1 row;"
        )
        say(
            "      - Part A's seam unit (real emcee) remains the authoritative F1-fix coverage."
        )

    pass_b = bool(pass_b1 and pass_b2)
    say(
        "    B VERDICT (lnprob AND blob identical per-theta on every row): "
        + ("PASS" if pass_b else "FAIL")
    )
    results["b_real_container"] = {
        "pass": bool(pass_b),
        "rows": rows_b,
        "f1_class_rows": f1_labels,
    }
    results["b2_sigma0_probe"] = {"pass": bool(pass_b2), "rows": rows_0}
    if not pass_b:
        _finish(results, "b_real_container_fail")

    # ----------------------------------------------------------------------
    # C. EMCEE GUARD -- decision function with fake versions (no other emcee)
    # ----------------------------------------------------------------------
    from pyorbit.samplers.pyorbit_emcee import _emcee_vectorize_gate

    say("### C. EMCEE GUARD (unit, fake versions)")
    cases = {
        "3.1.6": True,
        "3": True,
        "4.0": False,
        "2.9": False,
        "": False,
        "garbage": False,
    }
    pass_c = True
    import emcee as _emcee  # the REAL installed emcee, for the runtime gate line

    for v, want_ok in cases.items():
        okc, msg = _emcee_vectorize_gate(v)
        good = bool(okc) == want_ok
        pass_c = pass_c and good
        say(
            "    gate('{0}') -> ok={1} msg='{2}' {3}".format(
                v, okc, msg, "" if good else "<< MISMATCH"
            )
        )
    say(
        "    installed emcee={0} -> gate ok={1}".format(
            _emcee.__version__, bool(_emcee_vectorize_gate(_emcee.__version__)[0])
        )
    )
    say(
        "    C VERDICT: "
        + ("PASS" if pass_c else "FAIL")
        + " (deps unchanged: pyproject.toml not edited)"
    )
    results["c_guard"] = {"pass": bool(pass_c)}
    if not pass_c:
        _finish(results, "c_guard_fail")

    # ----------------------------------------------------------------------
    # D. FALLBACK-LOG THROTTLE -- 150 forced per-theta fallback calls
    # ----------------------------------------------------------------------
    import pyorbit.classes.model_container_abstract as mca

    say("### D. FALLBACK-LOG THROTTLE (M==1 path x150)")
    mca._LNPROB_BATCH_FALLBACK_COUNT["n"] = 0  # deterministic start (test only)
    buf = io.StringIO()
    old_stdout = sys.stdout
    sys.stdout = buf
    # ONE 1-D row per call -> M == 1 -> the throttled per-theta fallback
    # (the counter increments once per lnprob_batch call into _per_theta).
    try:
        for _ in range(150):
            mc.lnprob_batch([mid])
    finally:
        sys.stdout = old_stdout
    lines = [
        line
        for line in buf.getvalue().splitlines()
        if "lnprob_batch: per-theta fallback" in line
    ]
    pass_d = bool(len(lines) == 2 and bool(lines[0].strip()))
    say(
        "    log lines out of 150 calls: {0} (expect 2: occurrence 1 and 100)".format(
            len(lines)
        )
    )
    for line in lines:
        say("      " + line)
    say("    D VERDICT: " + ("PASS" if pass_d else "FAIL"))
    results["d_throttle"] = {"pass": pass_d, "n_lines": len(lines)}
    if not pass_d:
        _finish(results, "d_throttle_fail")

    say("### OVERALL: PASS (all gates)")
    _finish(results, "pass")


def _finish(results, status):
    results["status"] = status
    with open(os.path.join(OUT, "results.json"), "w") as fh:
        json.dump(results, fh, indent=2, default=str)
    say("### DONE -> {0}/results.json (status={1})".format(OUT, status))
    LOG.close()
    if status != "pass":
        sys.exit(1)


if __name__ == "__main__":
    main()
