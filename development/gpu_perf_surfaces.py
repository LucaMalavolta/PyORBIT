#!/usr/bin/env python
"""PyORBIT-GPU T1: wall-time performance surfaces (CPU leg vs GPU leg).

Companion to the value-parity probes (``development/gpu_parity_*.py``):
numerical CPU<->GPU parity is already proven on 7 surfaces; this script only
measures **wall time**. Parity numbers are re-checked as cheap sanity
guards (a perf leg that silently diverged in value would invalidate it).

Conventions (identical across all surfaces)
-------------------------------------------
* Both legs run in the ``pyorbit-gpu`` conda env (jax 0.6.2 + CUDA), float64.
* "CPU leg" = the production ``gpu: false`` code path
  (numpy Kepler/batman/spleaf-C/np.polyval, as in ``gpu_parity_*``).
* "GPU leg" = the production ``gpu: true`` code path invoked **as production
  invokes it** (eager jax where production is eager, jitted where production
  jits -- no extra wrapping is added to either leg).
* Warm-up is mandatory: ``WARMUP`` (3) un-timed calls absorb JIT tracing /
  NVRTC / CUDA-kernel compile / first-call overhead, then ``REPS`` (7) timed
  reps; the **median** is reported with min/max. The first-ever call (cold)
  is timed separately for every leg (Kepler & transit per the task).
* ``nvidia-smi`` memory/utilization is snapshotted before/after each leg
  group (shared-GPU contention evidence; the card is shared with ollama and
  must never be touched).

Surface list
------------
  kepler   pyorbit.subroutines.kepler_exo.kepler_compute_rv_deltabjd (numpy)
           vs the RVkeplerian.compute fused gpu: true path (single-graph
           cached_jit chain, the production call form); N in {100, 1000, 10000}.
  transit  batman.TransitModel.light_curve (CPU) vs
           pyorbit.subroutines.transit_gpu.circular_transit (jitted);
           N in {1000, 10000} epochs per light curve. Per-eval the fitted
           t0 perturbed by 1e-6 d on both legs (production semantics: fitted
           parameters change every eval, so batman re-runs its geometry).
  esp      spleaf Cov.loglike (C) vs pyorbit.subroutines.esp_gpu.esp_loglike
           (single-block CUDA LDL kernel); N in {1000, 5000}, nharm=4.
  tinygp   same GP as B6; N in {1000, 5000}; JAX_PLATFORMS=cpu vs cuda in two
           clean subprocesses (B6 two-subprocess pattern).
  trend    np.polyval (CPU) vs production PolynomialTrend.compute GPU Horner
           (staged x0_dev fast path); N in {1000, 10000}.
  mc       end-to-end float(mc(theta)) on the parity-suite synthetic RV set
           (600 pts, 10x the parity-probe n=60): gpu:false vs a fresh gpu:true
           container, production ordering (gpu_backend.initialize() BEFORE
           container build), each leg a clean subprocess.
  emcee    T2: emcee A/B, both legs in ONE subprocess, 24 walkers x 30 steps,
           identical seed/config, emcee 3.1.6 vectorize=False (per-walker 1-D
           calls -- the container's check_bounds is not 2-D safe; emcee 2.x
           convention).

Output discipline (T0): all verbose output goes to files under
``/tmp/opencode/pyorbit_perf/``; stdout carries compact RESULT lines only.

Modes
-----
  (default)                  orchestrator: runs legs, prints the summary table
  --mode inproc              in-process surfaces (kepler/transit/esp/trend)
  --mode tinygp-leg --platform cpu|gpu
  --mode mc-leg --gpu cpu|gpu [--npoints N]
  --mode emcee               T2 (cpu leg then gpu leg, same process)

Clean-run timing controls (additive; defaults preserve prior behaviour):
  --reps N                    timed reps per leg (default 7; clean run passes 11)
  --warmup N                  un-timed warmup calls before timing (min 3)
  --logdir DIR                output/log directory (default
                              /tmp/opencode/pyorbit_perf); the clean-run
                              orchestrator points this at its own OUT dir
  --kepler-n "100,1000,10000" comma list of N for the kepler surface
  --tinygp-n "1000,5000"      comma list of N for the tinygp legs (the clean
                              run extends to 10000/20000 when VRAM allows; the
                              per-rep cost is O(N^3) so reps auto-scale down
                              for N >= 10000, keeping the budget realistic)

Run (orchestrator):
  source ~/anaconda3/etc/profile.d/conda.sh && conda activate pyorbit-gpu
  python development/gpu_perf_surfaces.py
"""

import argparse
import glob
import json
import math
import os
import statistics
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEV_DIR = os.path.join(REPO, "development")
SELF = os.path.abspath(__file__)
LOGDIR = "/tmp/opencode/pyorbit_perf"

WARMUP = 3
REPS = 7


# --- logging plumbing (T0: keep stdout small) ------------------------------
def _ensure_logdir():
    try:
        os.makedirs(LOGDIR, exist_ok=True)
        return LOGDIR
    except Exception:
        alt = os.path.join("/tmp", "pyorbit_perf")
        os.makedirs(alt, exist_ok=True)
        return alt


def log_append(name, text):
    if not text.endswith("\n"):
        text += "\n"
    with open(os.path.join(LOGDIR, name), "a", encoding="utf-8") as fh:
        fh.write(text)


def nvidia_snapshot(tag):
    """One nvidia-smi snapshot line -> contention.log (shared-card evidence)."""
    try:
        out = (
            subprocess.run(
                [
                    "nvidia-smi",
                    "--query-gpu=memory.used,utilization.gpu",
                    "--format=csv,noheader",
                ],
                capture_output=True,
                text=True,
                timeout=15,
            )
            .stdout.strip()
            .replace("\n", " | ")
        )
    except Exception as exc:  # noqa: BLE001
        out = "nvidia-smi unavailable: {0!r}".format(exc)
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    log_append("contention.log", "{0}  {1}: {2}".format(stamp, tag, out))
    return out


def _sync(obj):
    """Force host completion of GPU/device work before stopping the clock.

    JAX calls (jitted OR eager) return after *launching* work; the kernel may
    still be running when ``perf_counter`` stops. cupy's esp_loglike already
    returns a synced Python float (``result.get()[0]``). Production pays this
    same sync (numpy/``float`` materialise the device array), so measuring it
    is faithful.
    """
    mod = type(obj).__module__ or ""
    if mod.startswith("jax") or mod.startswith("jaxlib"):
        try:
            obj.block_until_ready()
            return
        except Exception:  # noqa: BLE001
            pass
    if mod.startswith("cupy"):
        try:
            obj.get()
            return
        except Exception:  # noqa: BLE001
            pass
    # host value (numpy / float): nothing to sync


def bench(label, fn, warmup=None, reps=None, sanity_note=""):
    """Cold first call + ``warmup`` un-timed + ``reps`` timed; median/min/max.

    ``warmup=None`` / ``reps=None`` resolve to the module globals at call
    time so the ``--reps``/``--warmup`` CLI flags (applied in ``main()``)
    take effect, including in child-leg subprocesses that re-parse them.
    """
    if warmup is None:
        warmup = WARMUP
    if reps is None:
        reps = REPS
    t0 = time.perf_counter()
    out0 = fn()
    _sync(out0)
    cold = time.perf_counter() - t0
    for _ in range(warmup):
        _sync(fn())
    times = []
    for _ in range(reps):
        a = time.perf_counter()
        _sync(fn())
        b = time.perf_counter()
        times.append(b - a)
    rec = {
        "label": label,
        "cold_s": cold,
        "med_s": float(statistics.median(times)),
        "min_s": min(times),
        "max_s": max(times),
        # raw per-rep wall times (s) for downstream stats (IQR/std) in the
        # clean-run aggregator; med/min/max above remain for compat.
        "times_s": list(times),
        "reps": reps,
        "warmup": warmup,
    }
    log_append(
        "bench.log",
        json.dumps(rec, sort_keys=True)
        + (("    # " + sanity_note) if sanity_note else ""),
    )
    return rec, out0


def fmt(x, unit="ms"):
    if x is None:
        return "-"
    v = x * 1000.0 if unit == "ms" else x
    if v == 0:
        return "0.0"
    if v < 0.01:
        return "{0:.4f}".format(v)
    if v < 1000:
        return "{0:.2f}".format(v)
    return "{0:.0f}".format(v)


def bench_pair(res, tag, cpu_fn, gpu_fn, sanity_note=""):
    nvidia_snapshot(tag + " pre")
    rc, _ = bench(tag + "  [CPU]", cpu_fn, sanity_note=sanity_note)
    rg, _ = bench(tag + "  [GPU]", gpu_fn, sanity_note=sanity_note)
    nvidia_snapshot(tag + " post")
    if rg["med_s"] > 0:
        sp = rc["med_s"] / rg["med_s"]
    elif rc["med_s"] > 0:
        sp = float("inf")
    else:
        sp = float("nan")
    rc["speedup_gpu_vs_cpu"] = sp
    rg["speedup_gpu_vs_cpu"] = sp
    res.append(rc)
    res.append(rg)


def print_res_line(rec, surface, n, extra=""):
    print(
        "RESULT inproc surface={0} N={1} leg={2} cold_ms={3} med_ms={4} "
        "min_ms={5} max_ms={6} reps={7} {8}".format(
            surface,
            n,
            rec["label"].split("[")[-1].strip(" ]"),
            fmt(rec["cold_s"]),
            fmt(rec["med_s"]),
            fmt(rec["min_s"]),
            fmt(rec["max_s"]),
            rec["reps"],
            extra,
        )
    )


# --- in-process surfaces -----------------------------------------------------
def run_inproc(res, kepler_ns=None):
    if kepler_ns is None:
        kepler_ns = (100, 1000, 10000)
    import numpy as np

    np.seterr(all="ignore")

    import jax

    try:
        jax.config.update("jax_enable_x64", True)
    except Exception:
        pass
    log_append(
        "bench.log",
        "### inproc header: python={0} jax={1} backend={2} x64={3} seed-policy=fixed".format(
            sys.version.split()[0],
            jax.__version__,
            jax.default_backend(),
            bool(jax.config.jax_enable_x64),
        ),
    )

    # ------------------------------------------------------------- Kepler RV
    from pyorbit.subroutines import kepler_exo, kepler_exo_jax
    from types import SimpleNamespace

    import pyorbit.common.gpu_backend as gpu_backend
    from pyorbit.models import keplerian_modelling as kmod

    K, P, ECC, OMEGA, MLONG = 18.0, 21.0, 0.25, 120.0, 45.0

    # Fused production leg (F1/F2 fusion): an RVkeplerian built without
    # AbstractModel.__init__ (B5/trend pattern), invoked exactly as production
    # invokes it -- RVkeplerian.compute(obj, pv, dataset) with the GPU enabled
    # -- which routes through the single-graph cached_jit chain. The variant
    # below (K-only, no Tc, no inclination, not apodized) is numerically
    # identical to the numpy leg, so the A/B stays meaningful (same math:
    # omega=OMEGA, Omega=0.0 on both legs).
    kobj = object.__new__(kmod.RVkeplerian)
    kobj.use_time_inferior_conjunction = False
    kobj.use_mass = kobj.use_scaled_mass = kobj.use_stellar_scaled_mass = False
    kobj.compute_inclination = False
    kobj.compute_scaled_semimajor_axis = True
    # CPU leg passes omega_deg=OMEGA with Omega_deg=0 (default); mirror that
    # exactly so the fused chain and the numpy leg compute the same orbit.
    kobj.default_Omega = 0.0

    for N in kepler_ns:
        bjd = np.linspace(-600.0, 600.0, N)
        pv = dict(P=P, K=K, e=ECC, omega=OMEGA, mean_long=MLONG)
        dset = SimpleNamespace(x=bjd, x0=bjd, Tref=0.0)

        def cpu_fn(b=bjd):
            gpu_backend._enabled = False
            return kepler_exo.kepler_compute_rv_deltabjd(b, K, P, MLONG, ECC, OMEGA)

        # Fused GPU production path: gpu_enabled() -> the single-graph
        # cached_jit chain (what production calls under gpu: true).
        def gpu_fn(pv=pv, dset=dset):
            gpu_backend._enabled = True
            try:
                return kmod.RVkeplerian.compute(kobj, dict(pv), dset)
            finally:
                gpu_backend._enabled = False

        cv = np.asarray(cpu_fn()).ravel()
        gv = np.asarray(gpu_fn()).ravel()
        if cv.shape != gv.shape or not (
            np.all(np.isfinite(cv)) and np.all(np.isfinite(gv))
        ):
            log_append(
                "bench.log", "### KEPLER N={0} SANITY FAIL shape/finite".format(N)
            )
            print("RESULT inproc surface=kepler N={0} STATUS=sanity-fail".format(N))
            continue
        mdiff = float(np.max(np.abs(cv - gv)))
        tag = "kepler N={0}".format(N)
        bench_pair(
            res,
            tag,
            cpu_fn,
            gpu_fn,
            sanity_note="parity max|diff|={0:.3e} (numpy vs fused chain)".format(mdiff),
        )
        for r in res[-2:]:
            print_res_line(r, "kepler", N, "maxdiff={0:.2e} fused".format(mdiff))
        # reference leg: *bare jitted* RV stage (the fused production chain runs
        # this same stage inside its single cached_jit graph)
        kjit = jax.jit(kepler_exo_jax.kepler_compute_rv_deltabjd_jax)

        def gpu_fn_jit(b=bjd):
            return kjit(b, K, P, MLONG, ECC, OMEGA)

        rj, _ = bench(
            tag + "  [GPU-jit-reference]",
            gpu_fn_jit,
            sanity_note="bare-jitted stage reference (production path is the fused single-graph chain above)",
        )
        print_res_line(rj, "kepler-jitref", N, "jitted reference")
        res.append(rj)

    # ---------------------------------------------------------------- transit
    import batman
    import jax.numpy as jnp

    from pyorbit.subroutines.transit_gpu import circular_transit

    TP, TCT, RR, AARR, INC, U, EXP_D, SFAC = (
        3.5,
        0.15,
        0.10,
        10.0,
        87.0,
        (0.1, 0.3),
        (180.0 / 86400.0),
        2,
    )
    for N in (1000, 10000):
        times = np.linspace(0.0, 2.0 * TP, N)

        def _mk_params(i):
            p = batman.TransitParams()
            p.t0 = (
                TCT + i * 1e-6
            )  # per-eval fitted-parameter update (production semantics)
            p.per = TP
            p.rp = RR
            p.a = AARR
            p.inc = INC
            p.ecc = 0.0
            p.w = 90.0
            p.limb_dark = "quadratic"
            p.u = np.array(U)
            return p

        base = _mk_params(0)
        model = batman.TransitModel(
            base, times, supersample_factor=SFAC, exp_time=EXP_D, nthreads=1
        )
        k = [0]

        def cpu_fn(t=times, m=model, kc=k):
            p = _mk_params(kc[0])
            kc[0] += 1
            return m.light_curve(p)

        def gpu_fn(t=times, kc=k):
            out = circular_transit(
                jnp.asarray(t),
                TP,
                TCT + kc[0] * 1e-6,
                RR,
                AARR,
                INC * np.pi / 180.0,
                jnp.asarray(np.array(U)),
                EXP_D,
                samples=SFAC,
            )
            kc[0] += 1
            return out

        cv = np.asarray(cpu_fn()).ravel()
        gv = np.asarray(gpu_fn()).ravel()
        if cv.shape != gv.shape or not (
            np.all(np.isfinite(cv)) and np.all(np.isfinite(gv))
        ):
            log_append(
                "bench.log", "### TRANSIT N={0} SANITY FAIL shape/finite".format(N)
            )
            print("RESULT inproc surface=transit N={0} STATUS=sanity-fail".format(N))
            continue
        # B1 proved the quadrature floor is <= 1e-3; a perf leg must stay far
        # below that (different geometry on purpose: t0 perturbed per call).
        mdiff = float(np.max(np.abs(cv - gv)))
        tag = "transit N={0}".format(N)
        bench_pair(
            res,
            tag,
            cpu_fn,
            gpu_fn,
            sanity_note="max|diff| (perturbed geometry) = {0:.3e}".format(mdiff),
        )
        for r in res[-2:]:
            print_res_line(
                r, "transit", N, "sf={0} exp={1:.0f}s".format(SFAC, EXP_D * 86400.0)
            )

    # ------------------------------------------------------------------- ESP
    from spleaf import cov as sp_cov, term as sp_term

    from pyorbit.subroutines.esp_gpu import esp_loglike
    import pyorbit.common.gpu_backend as gpu_backend

    PROT, PDEC, OAMP, ALPHA, NHARM, SIGMA = 30.0, 100.0, 0.35, 1.0, 4, 5.0
    for N in (1000, 5000):
        t = np.linspace(0.0, 6.0 * PROT, N)
        y = np.sin(2.0 * np.pi * t / PROT) + np.random.default_rng(N).standard_normal(N)
        e2 = np.full(N, SIGMA * SIGMA)

        D = sp_cov.Cov(
            t,
            err=sp_term.Error(np.full(N, SIGMA)),
            GP=sp_term.ESPKernel(ALPHA, PROT, PDEC, OAMP, nharm=NHARM),
        )
        pvec = np.array([ALPHA, PROT, PDEC, OAMP])

        def cpu_fn(D=D, p=pvec, y=y):
            D.set_param(p, D.param)
            return D.loglike(y)

        def gpu_fn(t=t, y=y, e2=e2):
            return esp_loglike(t, y, e2, PROT, PDEC, OAMP, ALPHA, 0.0, NHARM)

        v_cpu = float(cpu_fn())
        # gpu_backend._enabled is required so esp_loglike's get_device() resolves
        gpu_backend._enabled = True
        try:
            v_gpu = float(gpu_fn())
        finally:
            gpu_backend._enabled = False
        if not (math.isfinite(v_cpu) and math.isfinite(v_gpu)):
            log_append(
                "bench.log", "### ESP N={0} SANITY FAIL non-finite loglike".format(N)
            )
            print("RESULT inproc surface=esp N={0} STATUS=sanity-fail".format(N))
            continue
        dlp = abs(v_cpu - v_gpu)
        # value parity needs gpu_backend._enabled so get_device() resolves a GPU
        gpu_backend._enabled = True
        res_before = len(res)
        tag = "esp N={0}".format(N)
        try:
            bench_pair(
                res,
                tag,
                cpu_fn,
                gpu_fn,
                sanity_note="|dloglike|={0:.3e} (B3 parity PASS<=1e-6)".format(dlp),
            )
        finally:
            gpu_backend._enabled = False
        for r in res[res_before:]:
            print_res_line(r, "esp", N, "dharm={0} dloglike={1:.2e}".format(NHARM, dlp))

    # ----------------------------------------------------------------- trend
    from types import SimpleNamespace

    from pyorbit.models.polynomial_trend import PolynomialTrend

    obj = object.__new__(PolynomialTrend)  # skip AbstractModel.__init__ (B5 pattern)
    obj.order = 5
    obj.starting_order = 1
    obj.time_interval = 1.0
    TREF = 2458000.0
    CS = [0.30, -0.20, 0.10, -0.05, 0.02]  # c1..c5 (moderate, no 1e4 blowup)
    pv = {"x_zero": TREF}
    for i in range(1, 6):
        pv["poly_c" + str(i)] = CS[i - 1]

    def _cpu_trend(ds):
        gpu_backend._enabled = False
        return PolynomialTrend.compute(obj, dict(pv), ds)

    def _gpu_trend(ds):
        gpu_backend._enabled = True
        try:
            return PolynomialTrend.compute(obj, dict(pv), ds)
        finally:
            gpu_backend._enabled = False

    for N in (1000, 10000):
        x = TREF - 5.0 + np.linspace(0.0, 10.0, N)
        ds_cpu = SimpleNamespace(x=x, Tref=TREF, x0_dev=None)
        ds_gpu = SimpleNamespace(x=x, Tref=TREF, x0_dev=jnp.asarray(x - TREF))

        cv = np.asarray(_cpu_trend(ds_cpu)).ravel()
        gv = np.asarray(_gpu_trend(ds_gpu)).ravel()
        if cv.shape != gv.shape or not (
            np.all(np.isfinite(cv)) and np.all(np.isfinite(gv))
        ):
            log_append(
                "bench.log", "### TREND N={0} SANITY FAIL shape/finite".format(N)
            )
            print("RESULT inproc surface=trend N={0} STATUS=sanity-fail".format(N))
            continue
        mdiff = float(np.max(np.abs(cv - gv)))
        tag = "trend N={0}".format(N)
        cpu_fn = lambda: _cpu_trend(ds_cpu)
        gpu_fn = lambda: _gpu_trend(ds_gpu)
        res_before = len(res)
        try:
            bench_pair(
                res,
                tag,
                cpu_fn,
                gpu_fn,
                sanity_note="parity max|diff|={0:.3e} (B5 PASS<=1e-10)".format(mdiff),
            )
        finally:
            gpu_backend._enabled = False
        for r in res[res_before:]:
            print_res_line(
                r, "trend", N, "order=5 x0_dev fast path, maxdiff={0:.2e}".format(mdiff)
            )

    return res


# --- subprocess legs ---------------------------------------------------------
def tinygp_leg_body(platform, ns=None):
    """--mode tinygp-leg: B6 setup (identical to the parity probe), timed here."""
    if ns is None:
        ns = (1000, 5000)
    import numpy as np

    np.seterr(all="ignore")
    import jax

    try:
        jax.config.update("jax_enable_x64", True)
    except Exception:
        pass
    log_append(
        "bench.log",
        "### tinygp-leg platform={0}: python={1} jax={2} backend={3}".format(
            platform, sys.version.split()[0], jax.__version__, jax.default_backend()
        ),
    )
    import tinygp
    from tinygp import kernels

    for N in ns:
        # tinygp log_probability is O(N^3) (Cholesky). For the big-N
        # extension points scale the timed reps down (keep >=3) so the
        # clean-run budget stays realistic; per-rep stats survive in times_s.
        reps_n = REPS if N < 10000 else max(3, int(REPS * 5000.0 / N))
        t = np.sort(np.random.default_rng(0).uniform(0.0, 10.0, N))
        y = np.random.default_rng(1).normal(size=N)
        gp = tinygp.GaussianProcess(
            kernels.ExpSquared(scale=1.5)
            * kernels.ExpSineSquared(scale=2.5, gamma=0.5),
            t,
            diag=1e-3,
            mean=0.0,
        )

        def call(gp=gp, y=y):
            return gp.log_probability(y)

        rec, val = bench("tinygp {0} N={1}".format(platform, N), call, reps=reps_n)
        v = float(np.asarray(val).ravel()[0]) if np.asarray(val).size else float(val)
        print(
            "RESULT tinygp platform={0} N={1} value={2} cold_ms={3} med_ms={4} "
            "min_ms={5} max_ms={6} reps={7}".format(
                platform,
                N,
                v,
                fmt(rec["cold_s"]),
                fmt(rec["med_s"]),
                fmt(rec["min_s"]),
                fmt(rec["max_s"]),
                rec["reps"],
            )
        )
    if platform == "cpu":
        # the paired GPU leg runs in the other subprocess; compare there
        pass
    return 0


def _build_mc_leg(n_points, gpu_flag):
    """Port of ``development/gpu_cpu_parity.py::_build_mc`` (n scalable,
    production ordering: ``gpu_backend.initialize`` BEFORE the container is
    built). Data/config generator is identical to the parity-suite synthetic."""
    import numpy as np
    import yaml

    import pyorbit.subroutines.input_parser as ip
    import pyorbit.common.gpu_backend as gpu_backend
    from pyorbit.classes.model_container_emcee import ModelContainerEmcee

    workdir = os.path.join(LOGDIR, "mc_leg" + ("_gpu" if gpu_flag else "_cpu"))
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
        "parameters": {"Tref": tref, "gpu": bool(gpu_flag)},
        "solver": {"emcee": {"npop_mult": 4, "nsteps": 10, "nburn": 5, "thin": 1}},
    }
    conf_path = os.path.join(workdir, "config.yaml")
    with open(conf_path, "w") as fh:
        yaml.safe_dump(conf, fh)

    config_in = ip.yaml_parser(conf_path)

    # ---- production ordering (pyorbit_run): initialize() BEFORE build -------
    gpu_backend.initialize(config_in["parameters"])

    mc = ModelContainerEmcee()
    ip.pars_input(config_in, mc, input_datasets=None)
    mc.model_setup()
    mc.boundaries_setup()
    mc.initialize_logchi2()
    return mc


def mc_leg_body(leg_flag, n_points):
    mc = _build_mc_leg(n_points, leg_flag)
    import numpy as np

    theta = np.array(
        [0.5 * (mc.bounds[i, 0] + mc.bounds[i, 1]) for i in range(mc.ndim)]
    )

    def call():
        return float(mc(theta))

    rec, val = bench(
        "mc {0} n={1}".format("gpu" if leg_flag else "cpu", n_points), call
    )
    v = float(val)
    if not math.isfinite(v):
        print(
            "RESULT mc leg={0} n={1} STATUS=non-finite value={2}".format(
                "gpu" if leg_flag else "cpu", n_points, v
            )
        )
        return 1
    import jax

    print(
        "RESULT mc leg={0} n={1} value={2} cold_ms={3} med_ms={4} min_ms={5} "
        "max_ms={6} reps={7} ndim={8} backend={9}".format(
            "gpu" if leg_flag else "cpu",
            n_points,
            v,
            fmt(rec["cold_s"]),
            fmt(rec["med_s"]),
            fmt(rec["min_s"]),
            fmt(rec["max_s"]),
            rec["reps"],
            mc.ndim,
            jax.default_backend(),
        )
    )
    return 0


def emcee_body():
    """T2: emcee A/B in ONE process (cpu leg first, then gpu leg)."""
    import numpy as np

    N_POINTS = 600
    NWALK_MULT, N_STEPS, SEED = 4, 30, 42

    def _run_leg(tag, gpu_flag):
        import emcee

        mc = _build_mc_leg(N_POINTS, gpu_flag)
        theta0 = np.array(
            [0.5 * (mc.bounds[i, 0] + mc.bounds[i, 1]) for i in range(mc.ndim)]
        )
        nwalkers = max(2 * (mc.ndim + 1), NWALK_MULT * mc.ndim)
        rng = np.random.default_rng(SEED)
        population = np.zeros((nwalkers, mc.ndim))
        for jj, val in enumerate(theta0):
            sigma = 1e-4 if val == 0.0 else abs(val) * 1e-4
            population[:, jj] = rng.normal(val, sigma, size=nwalkers)

        def log_prob(th):
            lp, ll = mc.log_priors_likelihood(np.atleast_1d(th))
            return float(np.asarray(lp) + np.asarray(ll))

        # NOTE: vectorize=False -- the container's check_bounds() is not safe
        # for emcee 3's default 2-D (nwalkers, ndim) batched call (1-D
        # per-walker convention, as emcee 2.x used in production).
        # emcee 3.1.6: no random_state kwarg on EnsembleSampler — it seeds
        # its internal RandomState from the *global* numpy state at
        # construction. Seed globally identically before each leg so both
        # legs draw the same proposal stream (fair A/B).
        np.random.seed(SEED)
        sampler = emcee.EnsembleSampler(
            nwalkers,
            mc.ndim,
            log_prob,
            vectorize=False,
        )
        # warm-up (un-timed): first eval of this process/leg
        sampler.run_mcmc(population, 5, progress=False, skip_initial_state_check=True)
        walls = []
        for _ in range(2):  # 2 timed runs of N_STEPS steps each
            a = time.perf_counter()
            sampler.run_mcmc(
                population, N_STEPS, progress=False, skip_initial_state_check=True
            )
            walls.append(time.perf_counter() - a)
        sps = [N_STEPS / w for w in walls]
        med_sp = float(statistics.median(sps))
        print(
            "RESULT emcee leg={0} n={1} ndim={2} nwalkers={3} steps_run={4} "
            "median_steps_per_s={5} walls_s={6} value_start={7}".format(
                tag,
                N_POINTS,
                mc.ndim,
                nwalkers,
                2 * N_STEPS,
                "{0:.3f}".format(med_sp),
                " ".join("{0:.2f}".format(w) for w in walls),
                theta0[0],
            )
        )
        return med_sp

    log_append(
        "bench.log",
        "### emcee T2: n={0} walkers={1}*ndim steps={2}/run x2 seed={3} vectorize=False".format(
            N_POINTS, NWALK_MULT, N_STEPS, SEED
        ),
    )
    sps_cpu = _run_leg("cpu", gpu_flag=False)
    sps_gpu = _run_leg("gpu", gpu_flag=True)
    print(
        "RESULT emcee speedup_gpu_vs_cpu_steps_per_s={0:.4g}".format(sps_gpu / sps_cpu)
    )
    return 0


# --- subprocess plumbing (B6 child-env pattern) -----------------------------
def child_env(platform, extra=None):
    env = dict(os.environ)
    env["JAX_ENABLE_X64"] = "true"
    if platform == "gpu":
        # jax 0.6.2's "gpu" alias resolves to rocm and aborts if rocm is
        # absent; name the concrete cuda backend (B6 finding).
        env["JAX_PLATFORMS"] = "cuda"
        # F2 (Fusion): shared allocator defaults (preallocate true/0.30 by
        # default; false/0.10 auto-fallback when < ~10 GiB is free) instead of
        # the hard-coded false/0.1 that OOM'd the N=10000 leg. PYORBIT_XLA_*
        # and existing XLA_PYTHON_CLIENT_* env still win.
        from pyorbit.common import gpu_backend as _gb

        env.update(_gb.xla_allocator_envs())
    else:
        env["JAX_PLATFORMS"] = "cpu"
    # NVRTC/CUDA fallback (B6 pattern)
    prefix = os.environ.get("CONDA_PREFIX", "")
    if prefix:
        extra_paths = [os.path.join(prefix, "lib")]
        for pydir in glob.glob(os.path.join(prefix, "lib", "python*")):
            extra_paths += glob.glob(
                os.path.join(pydir, "site-packages", "nvidia", "*", "lib")
            )
        extra_paths += glob.glob(os.path.join(prefix, "nvidia", "*", "lib"))
        present = [p for p in extra_paths if os.path.isdir(p)]
        if present:
            env["LD_LIBRARY_PATH"] = os.pathsep.join(
                present + [env.get("LD_LIBRARY_PATH", "")]
            )
    if extra:
        env.update(extra)
    return env


def run_subleg(mode_args, tag, log_name, timeout=1800, platform=None, extra_env=None):
    env = child_env(platform, extra=extra_env) if platform else dict(os.environ)
    if extra_env:
        env.update(extra_env)
    r = subprocess.run(
        [sys.executable, SELF] + mode_args,
        env=env,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    with open(os.path.join(LOGDIR, log_name), "a") as fh:
        fh.write(
            "## {0} (exit={1})\n## --- stdout ---\n{2}\n## --- stderr ---\n{3}\n".format(
                tag, r.returncode, r.stdout, r.stderr
            )
        )
    return r


def extract_results(r, prefix="RESULT "):
    out = []
    for line in (r.stdout or "").splitlines():
        if line.startswith(prefix):
            out.append(line)
    return out


# --- orchestrator ------------------------------------------------------------
def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--mode",
        default="all",
        choices=["all", "inproc", "tinygp-leg", "mc-leg", "emcee"],
    )
    parser.add_argument("--platform", default=None, choices=[None, "cpu", "gpu"])
    parser.add_argument("--gpu", default="cpu", choices=["cpu", "gpu"])
    parser.add_argument("--npoints", type=int, default=600)
    # --- clean-run timing controls (additive; defaults keep old behaviour) ---
    parser.add_argument(
        "--reps", type=int, default=7, help="timed reps per leg (default 7)"
    )
    parser.add_argument(
        "--warmup",
        type=int,
        default=3,
        help="un-timed warm-up calls before timing (min 3)",
    )
    parser.add_argument(
        "--logdir",
        default=None,
        help="log directory (default /tmp/opencode/pyorbit_perf)",
    )
    parser.add_argument(
        "--kepler-n",
        default="100,1000,10000",
        help="comma-separated N list for the kepler surface",
    )
    parser.add_argument(
        "--tinygp-n",
        default="1000,5000",
        help="comma-separated N list for the tinygp legs",
    )
    args = parser.parse_args(argv)

    global WARMUP, REPS, LOGDIR
    WARMUP = max(3, args.warmup)
    REPS = max(1, args.reps)
    if args.logdir:
        LOGDIR = args.logdir
    kepler_ns = tuple(int(x) for x in args.kepler_n.split(",") if x.strip())
    tinygp_ns = tuple(int(x) for x in args.tinygp_n.split(",") if x.strip())
    # flags forwarded to child-leg subprocesses so their globals match too
    common_flags = ["--reps", str(REPS), "--warmup", str(WARMUP), "--logdir", LOGDIR]

    logdir = _ensure_logdir()
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    log_append(
        "bench.log",
        "### === session {0} mode={1} python={2} ===".format(
            stamp, args.mode, sys.version.split()[0]
        ),
    )

    if args.mode == "tinygp-leg":
        # Pin the JAX platform for standalone runs too: JAX reads
        # JAX_PLATFORMS before its first device init, and tinygp_leg_body
        # imports jax lazily, so setting it here binds correctly. Without
        # this a lone `--mode tinygp-leg` run falls to the default backend
        # (cuda) even when --platform cpu is passed (observed cuSolver fail
        # at N=5000: mislabelled leg).
        plat = args.platform or "cpu"
        if plat == "gpu":
            os.environ["JAX_PLATFORMS"] = "cuda"
            # F2 (Fusion): shared allocator defaults (true/0.30 by default;
            # false/0.10 auto-fallback when < ~10 GiB is free), replacing the
            # hard-coded false/0.1 that OOM'd the N=10000 leg. PYORBIT_XLA_*
            # and any already-set XLA_PYTHON_CLIENT_* env still win.
            from pyorbit.common import gpu_backend as _gb

            os.environ.update(_gb.xla_allocator_envs())
        else:
            os.environ["JAX_PLATFORMS"] = "cpu"
        os.environ["JAX_ENABLE_X64"] = "true"
        return tinygp_leg_body(plat, tinygp_ns)
    if args.mode == "mc-leg":
        return mc_leg_body(args.gpu == "gpu", args.npoints)
    if args.mode == "emcee":
        return emcee_body()
    if args.mode == "inproc":
        run_inproc([], kepler_ns)
        return 0

    # ---- mode all: orchestrator -------------------------------------------
    res = []
    print("GPU perf surfaces (orchestrator); logs in {0}".format(logdir))
    nvidia_snapshot("orchestrator start")

    # 1) in-process surfaces (kepler / transit / esp / trend),
    #    subprocess-isolated: jaxlib preallocates most of the card's *free*
    #    VRAM per process (~5.2 GiB on this shared card). A long-lived
    #    parent would keep that arena and starve every child that launches
    #    afterwards (observed: all four later legs then failed CUDA context
    #    init while nvidia-smi read 31999 MiB / ~97 %). Isolation releases
    #    the arena before leg 2 starts.
    bench_path = os.path.join(logdir, "bench.log")

    def _bench_count():
        try:
            with open(bench_path, encoding="utf-8") as fh:
                return sum(1 for _ in fh)
        except OSError:
            return 0

    n0 = _bench_count()
    nvidia_snapshot("inproc pre")
    r_in = run_subleg(
        ["--mode", "inproc", "--kepler-n", args.kepler_n] + common_flags,
        "inproc",
        "inproc_full.log",
    )
    nvidia_snapshot("inproc post")
    if r_in.returncode != 0:
        print(
            "### inproc FAILED (exit={0}); log {1}".format(
                r_in.returncode, os.path.join(logdir, "inproc_full.log")
            )
        )
    else:
        for line in (r_in.stdout or "").splitlines():
            if line.startswith("RESULT "):
                print("  " + line)
        # Rebuild the res records from the JSON appended to bench.log
        # (bench() writes one JSON record per leg).
        try:
            with open(bench_path, encoding="utf-8") as fh:
                new_lines = fh.readlines()[n0:]
        except OSError:
            new_lines = []
        for line in new_lines:
            line = line.strip()
            if not line.startswith("{"):
                continue
            try:
                rec = json.loads(line.split("    #")[0])
            except ValueError:
                continue
            if rec.get("label", "").startswith(
                ("kepler ", "transit ", "esp ", "trend ")
            ):
                res.append(rec)
        # post-hoc speedup per cpu/gpu pair (same rule as bench_pair)
        pairs = {}
        for rec in res:
            lab = rec["label"]
            key = (
                lab.replace("  [CPU]", "")
                .replace("  [GPU]", "")
                .replace("  [GPU-jit-reference]", "")
                .strip()
            )
            if "  [CPU]" in lab:
                pairs.setdefault(key, {})["cpu"] = rec
            elif "  [GPU]" in lab:
                pairs.setdefault(key, {})["gpu"] = rec
        for parts in pairs.values():
            rc = parts.get("cpu")
            rg = parts.get("gpu")
            if rc is None or rg is None:
                continue
            if rg["med_s"] > 0:
                sp = rc["med_s"] / rg["med_s"]
            elif rc["med_s"] > 0:
                sp = float("inf")
            else:
                sp = float("nan")
            rc["speedup_gpu_vs_cpu"] = sp
            rg["speedup_gpu_vs_cpu"] = sp

    # 2) tinygp legs (two clean subprocesses, B6 pattern)
    for plat in ("cpu", "gpu"):
        nvidia_snapshot("tinygp-{0} pre".format(plat))
        r = run_subleg(
            ["--mode", "tinygp-leg", "--platform", plat, "--tinygp-n", args.tinygp_n]
            + common_flags,
            "tinygp-leg-{0}".format(plat),
            "tinygp_{0}_leg.log".format(plat),
            platform=plat,
        )
        nvidia_snapshot("tinygp-{0} post".format(plat))
        if r.returncode != 0:
            print(
                "### tinygp-leg {0} FAILED (exit={1}); log {2}".format(
                    plat,
                    r.returncode,
                    os.path.join(logdir, "tinygp_{0}_leg.log".format(plat)),
                )
            )
            continue
        for line in extract_results(r):
            print("  " + line)

    # 3) end-to-end mc legs (fresh containers, production ordering)
    for leg in ("cpu", "gpu"):
        nvidia_snapshot("mc-leg-{0} pre".format(leg))
        r = run_subleg(
            ["--mode", "mc-leg", "--gpu", leg, "--npoints", str(args.npoints)]
            + common_flags,
            "mc-leg-{0}".format(leg),
            "mc_{0}_leg.log".format(leg),
        )
        nvidia_snapshot("mc-leg-{0} post".format(leg))
        if r.returncode != 0:
            print(
                "### mc-leg {0} FAILED (exit={1}); log {2}".format(
                    leg,
                    r.returncode,
                    os.path.join(logdir, "mc_{0}_leg.log".format(leg)),
                )
            )
            continue
        for line in extract_results(r):
            print("  " + line)

    # 4) T2 emcee A/B (single process, both legs) -- optional, VRAM-gated
    nvidia_snapshot("emcee T2 pre")
    r = run_subleg(
        ["--mode", "emcee"] + common_flags, "emcee-T2", "emcee_T2.log", platform=None
    )
    nvidia_snapshot("emcee T2 post")
    if r.returncode != 0:
        print(
            "### emcee T2 FAILED/OOM (exit={0}); log {1}".format(
                r.returncode, os.path.join(logdir, "emcee_T2.log")
            )
        )
    else:
        for line in extract_results(r):
            print("  " + line)

    # 5) final compact table (stdout + file) --------------------------------
    lines = [
        "surface/N                    leg        cold_ms  med_ms  min_ms  max_ms  reps  speedup(G/C)"
    ]
    rows_sorted = sorted(
        res, key=lambda r: (r["label"], 0 if "[CPU]" in r["label"] else 1)
    )
    for r in rows_sorted:
        sp = r.get("speedup_gpu_vs_cpu")
        sp_s = (
            "-"
            if sp is None or math.isnan(sp)
            else ("{0:.2f}x".format(sp) if math.isfinite(sp) else ">1e3x")
        )
        lines.append(
            "{0:<28s} {1:<8s} {2:>9.3f} {3:>8.3f} {4:>8.3f} {5:>8.3f} {6:>5d}  {7}".format(
                r["label"].split("  [")[0].strip()[:28],
                r["label"].split("[")[-1].strip(" ]"),
                r["cold_s"] * 1000.0,
                r["med_s"] * 1000.0,
                r["min_s"] * 1000.0,
                r["max_s"] * 1000.0,
                r["reps"],
                sp_s,
            )
        )
    print()
    print("\n".join(lines))
    log_append("table.log", "\n".join(lines) + "\n")
    with open(os.path.join(logdir, "results.json"), "w") as fh:
        json.dump(res, fh, indent=2, default=str)
    print("log dir: {0}".format(logdir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
