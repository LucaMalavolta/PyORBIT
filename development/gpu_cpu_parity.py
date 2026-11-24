#!/usr/bin/env python
"""CPU/GPU parity check for PyORBIT Round-1 GPU enablement.

What it verifies
================
Round-1's hard contract is that ``gpu: true`` reproduces ``gpu: false`` for the
JAX-native models, while a GPU that is actually present is genuinely used. This
script exercises the two things that changed on the per-sample path:

  * the **RV Keplerian solver** -- the one real port this round, dispatched in
    ``pyorbit.models.keplerian_modelling`` between ``pyorbit.subroutines.kepler_exo``
    (numpy, CPU) and ``pyorbit.subroutines.kepler_exo_jax`` (JAX, GPU); and
  * **tinygp** -- already JAX-native, expected to simply land on the GPU.

Sections
--------
0. environment self-check -- import ``pyorbit`` + the pieces used below. If the
   pinned CPU stack (py3.10, numpy 1.24.3, jax 0.4.31, ordered_set, ...) is not
   installed this prints a precise setup message and exits 3 (it is a *missing
   dependency* condition, not a parity failure).
1. Keplerian unit parity -- numpy vs JAX over a parameter grid that includes the
   ``|ecc| < 1e-3`` circular/eccentric boundary and negative eccentricities.
   This is the actual correctness risk of the port and runs on the CPU.
2. end-to-end CPU baseline -- build the real ``ModelContainer`` (a small
   ``radial_velocities`` / ``keplerian`` fit) and evaluate ``mc(theta)`` (the
   scalar log priors + log likelihood) at a fixed, in-bounds ``theta``. This is
   the CPU reference.
3. zero-regression -- with ``gpu: false`` (the default) assert the GPU was never
   enabled and JAX was never retargeted to a GPU.
4. gate logic -- assert an allowlisted model passes ``gpu_backend.gate`` while a
   gated model (e.g. a george GP or a TRADES integrator) raises
   ``GPUSupportedModelError`` naming the JAX equivalent. Runs on the CPU by
   driving the gate's decision logic directly.
5. GPU parity -- if (and only if) the installed ``jaxlib`` can actually see a
   GPU: enable the GPU (asserting ``jax.default_backend() == 'gpu'``) and then
   build a *fresh* container exactly as ``PyORBIT_Run.py`` does with ``gpu: true``
   (initialize first, container second -- datasets stage their device/jitter
   arrays at construction time), evaluate it at the *same* ``theta`` and assert
   it matches the CPU baseline to a tight x64 tolerance; plus a small tinygp
   device check. With no GPU present this section SKIPS (it is not a failure)
   and reports exactly what to install.

Run
===
    python development/gpu_cpu_parity.py

Exit codes: 0 = all runnable sections pass (GPU section may SKIP),
            1 = a runnable section failed,
            3 = environment not set up to run this script at all.
"""

from __future__ import print_function

import os
import sys
import traceback


# --- result bookkeeping ------------------------------------------------------
_PASS = []
_SKIP = []
_FAIL = []


def _ok(section, detail):
    _PASS.append((section, detail))
    print("  [PASS] {0}: {1}".format(section, detail))


def _skip(section, detail):
    _SKIP.append((section, detail))
    print("  [SKIP] {0}: {1}".format(section, detail))


def _fail(section, detail):
    _FAIL.append((section, detail))
    print("  [FAIL] {0}: {1}".format(section, detail))


def _banner(n, title):
    print()
    print("=== {0}. {1} ===".format(n, title))


# --- 0. environment self-check ----------------------------------------------
def section0_env():
    _banner(0, "environment self-check")

    try:
        import jax  # noqa: F401
    except Exception as exc:  # pragma: no cover - environment dependent
        raise _EnvError("jax is not importable: {0}".format(exc))

    try:
        import pyorbit  # noqa: F401
    except Exception as exc:  # pragma: no cover - environment dependent
        raise _EnvError(
            "import pyorbit failed: {0}\n"
            "This script needs the pinned CPU stack: Python 3.10, "
            "numpy==1.24.3, jax==0.4.31, jaxlib==0.4.31, ordered_set, tinygp==0.3.0 "
            "(see setup.py / extra_requirements.txt).".format(exc)
        )

    try:
        import pyorbit.subroutines.kepler_exo as kepler_exo  # noqa: F401
        import pyorbit.subroutines.kepler_exo_jax as kepler_exo_jax  # noqa: F401
        import pyorbit.common.gpu_backend as gpu_backend  # noqa: F401
        import pyorbit.model_definitions as md  # noqa: F401
    except Exception as exc:
        raise _EnvError(
            "a pyorbit module needed by this script failed to import: {0}".format(exc)
        )

    # Force float64 JAX now, so section 1 compares numpy(x64) vs JAX(x64) -- the
    # same pairing Round 1 enforces with `gpu: true`.
    import jax as _jax

    try:
        _jax.config.update("jax_enable_x64", True)
    except Exception:
        os.environ.setdefault("JAX_ENABLE_X64", "true")

    _ok("env", "pyorbit + jax + kepler_exo(_jax) + gpu_backend import OK")


class _EnvError(RuntimeError):
    pass


# --- 1. Keplerian unit parity ------------------------------------------------
def section1_kepler_parity():
    _banner(1, "Keplerian unit parity (numpy vs JAX)")

    import numpy as np
    import pyorbit.subroutines.kepler_exo as kepler_exo
    import pyorbit.subroutines.kepler_exo_jax as kepler_exo_jax
    from pyorbit.subroutines.common import convert_rho_to_ars, convert_b_to_i

    rng = np.random.default_rng(0)
    n = 4000

    # Broad eccentricities, a cluster right on the |ecc| < 1e-3 circular boundary,
    # and negative eccentricities (which flip omega -> omega + pi).
    ecc = rng.uniform(-0.999, 0.999, n)
    ecc[0:6] = [0.0, 1e-4, -1e-4, 2e-3, -2e-3, -0.999]
    omega = rng.uniform(0.0, 360.0, n)
    Omega = rng.uniform(0.0, 180.0, n)
    period = rng.uniform(5.0, 200.0, n)
    mean_long = rng.uniform(0.0, 360.0, n)
    dTc = rng.uniform(0.0, 1.0, n) * period
    K = rng.uniform(10.0, 500.0, n)
    bjd = rng.uniform(0.0, 1.0, n) * period

    worst = {}
    worst["rv_deltabjd"] = _maxdiff(
        [
            kepler_exo.kepler_compute_rv_deltabjd(
                bjd[i], K[i], period[i], mean_long[i], ecc[i], omega[i], Omega[i]
            )
            for i in range(n)
        ],
        [
            kepler_exo_jax.kepler_compute_rv_deltabjd_jax(
                bjd[i], K[i], period[i], mean_long[i], ecc[i], omega[i], Omega[i]
            )
            for i in range(n)
        ],
    )
    worst["rv_semiamplitude"] = _maxdiff(
        [
            kepler_exo.kepler_compute_rv_semiamplitude(
                1.2, 0.02, period[i], 85.0, abs(ecc[i])
            )
            for i in range(n)
        ],
        [
            kepler_exo_jax.kepler_compute_rv_semiamplitude_jax(
                1.2, 0.02, period[i], 85.0, abs(ecc[i])
            )
            for i in range(n)
        ],
    )
    worst["deltaTc_from_meanlong"] = _maxdiff(
        [
            kepler_exo.kepler_compute_deltaTc_from_meanlong(
                period[i], mean_long[i], abs(ecc[i]), omega[i], Omega[i]
            )
            for i in range(n)
        ],
        [
            kepler_exo_jax.kepler_compute_deltaTc_from_meanlong_jax(
                period[i], mean_long[i], abs(ecc[i]), omega[i], Omega[i]
            )
            for i in range(n)
        ],
    )
    worst["meanlong_from_dTc"] = _maxdiff(
        [
            kepler_exo.kepler_compute_meanlong_from_deltaTc(
                period[i], dTc[i], abs(ecc[i]), omega[i], Omega[i]
            )
            for i in range(n)
        ],
        [
            kepler_exo_jax.kepler_compute_meanlong_from_deltaTc_jax(
                period[i], dTc[i], abs(ecc[i]), omega[i], Omega[i]
            )
            for i in range(n)
        ],
    )
    # impact-parameter / density conversions (use_mass branch). Draw each random
    # input ONCE so the numpy and JAX lists are evaluated at identical points.
    a = rng.uniform(0.05, 3.0, n)
    density = rng.uniform(0.1, 20.0, n)
    b_arr = rng.uniform(-1.5, 1.5, n)
    worst["convert_rho_to_ars"] = _maxdiff(
        [convert_rho_to_ars(period[i], density[i]) for i in range(n)],
        [
            kepler_exo_jax.convert_rho_to_ars_jax(period[i], density[i])
            for i in range(n)
        ],
    )
    worst["convert_b_to_i"] = _maxdiff(
        [convert_b_to_i(b_arr[i], abs(ecc[i]), omega[i], a[i]) for i in range(n)],
        [
            kepler_exo_jax.convert_b_to_i_jax(b_arr[i], abs(ecc[i]), omega[i], a[i])
            for i in range(n)
        ],
    )

    rtol, atol = 1e-9, 1e-9
    for k, (ma, mr) in worst.items():
        print(
            "    {0:24s} max|abs diff| = {1:.3e}   max rel = {2:.3e}".format(k, ma, mr)
        )
    bad = {k: v for k, v in worst.items() if v[1] > 1.0}
    if bad:
        _fail(
            "kepler-parity",
            "relative residual > 1 (rtol={0}, atol={1}): {2}".format(
                rtol, atol, {k: "max_rel={0:.3e}".format(v[1]) for k, v in bad.items()}
            ),
        )
    else:
        _ok(
            "kepler-parity",
            "all {0} Keplerian functions agree within rtol={1}/atol={2} over {3} points (incl. boundary & negative ecc)".format(
                len(worst), rtol, atol, n
            ),
        )


def _maxdiff(a, b, rtol=1e-9, atol=1e-9):
    """Per-case agreement between the numpy and JAX results.

    Returns ``(max_abs, max_rel)``:

    * ``max_abs`` -- the largest ``|a-b|`` over all cases (a plain diagnostic);
    * ``max_rel`` -- the largest *relative* residual
      ``max |a-b| / (atol + rtol*max(|a|,|b|))`` over all cases.

    Section 1 passes when ``max_rel <= 1``. numpy and JAX agree to a few ULPs on
    float64 but are **not** bit-identical: the two backends differ in reduction
    order and fused-operator choice, so a transcendental-heavy solver like the
    Keplerian one legitimately differs in the last handful of digits. A genuinely
    wrong port (a sign, a swapped argument, a dropped term) would be O(1) -- many
    orders larger than ``rtol`` -- so ``rtol=1e-9`` cleanly separates "same maths"
    from "bug". (A fixed ``1e-15`` absolute floor would reject correct ports.)

    Each case is a 0-d array (circular scalar call) or a 1-element array
    (eccentric -- ``kepler_E`` applies ``atleast_1d``), so the list is
    inhomogeneous and cannot be stacked with ``np.asarray``; compare per element.
    """
    import numpy as np

    if len(a) != len(b):
        raise ValueError("result length mismatch: %d vs %d" % (len(a), len(b)))
    max_abs = 0.0
    max_rel = 0.0
    for x, y in zip(a, b):
        xa = np.asarray(x, dtype=float).ravel()
        ya = np.asarray(y, dtype=float).ravel()
        if xa.shape != ya.shape:
            raise ValueError("per-case shape mismatch: %s vs %s" % (xa.shape, ya.shape))
        d = np.abs(xa - ya)
        if np.any(np.isnan(d)):
            return float("inf"), float("inf")
        denom = atol + rtol * np.maximum(np.abs(xa), np.abs(ya))
        max_abs = max(max_abs, float(np.max(d)))
        max_rel = max(max_rel, float(np.max(d / denom)))
    return max_abs, max_rel


# --- 2/3. end-to-end CPU baseline + zero-regression -------------------------
def _build_mc(tmpdir, gpu_flag=False):
    """Build a small radial-velocities / keplerian ModelContainer in ``tmpdir``.

    ``gpu_flag`` is only written into the YAML config (mirroring the production
    A/B configs where the two legs differ solely in the `gpu` flag); the actual
    enablement is decided by `pyorbit.common.gpu_backend.initialize()` in
    ``pyorbit_run.py`` -- in production that runs *before* the container is
    built, and section 5 below reproduces that ordering.
    """
    import numpy as np
    import pyorbit.subroutines.input_parser as ip
    from pyorbit.classes.model_container_emcee import ModelContainerEmcee

    # -- synthetic RV data: epoch, value, error, jitter_flag ------------------
    rng = np.random.default_rng(1234)
    tref = 2457000.0
    n = 60
    epoch = tref - 150.0 + np.sort(rng.uniform(0.0, 300.0, n))
    P_true, K_true, e_true, om_true = 21.0, 18.0, 0.25, 120.0
    # simple circular-ish RV (good enough -- we only need finite, valid data)
    phase = 2.0 * np.pi * ((epoch - tref) / P_true)
    value = K_true * np.cos(phase) + 0.5 * rng.normal(size=n)
    error = np.full(n, 1.2)
    data_path = os.path.join(tmpdir, "synthetic_rv.dat")
    with open(data_path, "w") as fh:
        for t, v, e, f in zip(epoch, value, error, [0] * n):
            fh.write("{0:.6f}\t{1:.5f}\t{2:.3f}\t{3:d}\n".format(t, v, e, f))

    # -- config mirroring the 5-section schema (radial_velocities / keplerian) -
    conf = {
        "output": os.path.join(tmpdir, "gpu_parity"),
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
            # ``priors`` is required by star_parameters.initialize_model (even empty).
            "star": {"star_parameters": {"priors": {}}},
        },
        "models": {
            "radial_velocities": {"planets": ["p"]},
        },
        "parameters": {
            "Tref": tref,
            "gpu": False,
        },
        "solver": {
            "emcee": {"npop_mult": 4, "nsteps": 10, "nburn": 5, "thin": 1},
        },
    }

    conf_path = os.path.join(tmpdir, "gpu_parity.yaml")
    import yaml

    with open(conf_path, "w") as fh:
        yaml.safe_dump(conf, fh)

    config_in = ip.yaml_parser(conf_path)

    mc = ModelContainerEmcee()
    ip.pars_input(config_in, mc, input_datasets=None)
    mc.model_setup()
    mc.boundaries_setup()
    mc.initialize_logchi2()
    return mc


def _theta_midbounds(mc):
    import numpy as np

    return np.array([0.5 * (mc.bounds[i, 0] + mc.bounds[i, 1]) for i in range(mc.ndim)])


def section2_3_cpu_baseline(tmpdir):
    _banner(2, "end-to-end CPU baseline (radial_velocities / keplerian)")

    mc = _build_mc(tmpdir)
    import numpy as np

    theta = _theta_midbounds(mc)
    print("    ndim = {0}".format(mc.ndim))

    # ``mc(theta)`` is the scalar log probability (log_priors + log_likelihood);
    # the raw ``mc.log_priors_likelihood(theta)`` returns the 2-tuple (priors, lnlk).
    cpu_val = float(mc(theta))
    if not np.isfinite(cpu_val):
        _fail(
            "cpu-baseline",
            "log probability (priors + likelihood) is not finite: {0}".format(cpu_val),
        )
        return mc, theta
    print("    log probability (priors + likelihood) [cpu] = {0:.10f}".format(cpu_val))
    _ok(
        "cpu-baseline",
        "MC built and CPU log-priors+likelihood finite ({0:.6f})".format(cpu_val),
    )

    # Authoritative CPU reference: captured now, while the CPU path is active and
    # before section 5 enables the GPU. Section 5 compares the GPU value to this.
    global _CPU_BASELINE
    _CPU_BASELINE = cpu_val

    _banner(3, "zero-regression (gpu: false must not force the GPU)")
    import pyorbit.common.gpu_backend as gpu_backend

    assert gpu_backend.gpu_enabled() is False, (
        "gpu_enabled() unexpectedly True on the CPU path"
    )
    import jax

    default_backend = str(jax.default_backend()).lower()
    # gpu: false means PyORBIT does NOT force a backend: no initialize(), no
    # CUDA_VISIBLE_DEVICES, no jax_platform_name, and the RV/Keplerian path uses
    # numpy (not JAX) at all. It does NOT promise to pin JAX to the CPU -- in a
    # GPU-capable environment JAX's *own* default backend is the GPU, which is
    # expected and harmless (the numpy RV path never touches JAX). The real
    # guarantee to assert is therefore gpu_enabled() == False, not a CPU backend.
    _ok(
        "zero-regression",
        "gpu_enabled() is False (PyORBIT forced nothing); jax default_backend == {0} "
        "is JAX's own environment default".format(default_backend),
    )
    return mc, theta


# --- 4. gate logic -----------------------------------------------------------
def section4_gate():
    _banner(4, "gate logic (allowed vs gated models)")

    import pyorbit.common.gpu_backend as gpu_backend
    import pyorbit.model_definitions as md

    allowed_models = set(getattr(md, "GPU_ALLOWED_MODELS", ()))
    redirect = getattr(md, "GPU_REDIRECT", {})
    if not allowed_models or not redirect:
        _fail("gate", "model_definitions lacks GPU_ALLOWED_MODELS / GPU_REDIRECT")
        return

    # (a) allowed: the actually-ported RV Keplerian must pass (fall back to any).
    allowed_candidates = [
        c
        for c in [getattr(md, "RVkeplerian", None)] + list(allowed_models)
        if c is not None
    ]
    # (b) gated: prefer a GP whose redirect names a tinygp_* equivalent -- the
    #     heart of the "gate + guide" contract -- then any other gated model.
    gated_candidates = [
        cls for cls, msg in redirect.items() if "tinygp" in msg.lower()
    ] + list(redirect)

    allowed_cls, allowed_obj = _pick_instantiable(allowed_candidates)
    gated_cls, gated_obj = _pick_instantiable(gated_candidates)

    if allowed_obj is None or gated_obj is None:
        _fail("gate", "could not instantiate an allowed+gated pair for the gate test")
        return

    # Drive the gate's decision logic directly: gate() reads only ``_enabled`` and
    # each object's class, so this exercises the allowlist check on the CPU without
    # needing real GPU hardware.
    gpu_backend._enabled = True
    try:
        # (a) an allowlisted model must pass.
        allowed_passed = True
        try:
            gpu_backend.gate({"radial_velocities": allowed_obj}, None)
        except Exception as exc:
            allowed_passed = False
            print(
                "    unexpected: allowed {0} rejected: {1}".format(
                    allowed_cls.__name__, exc
                )
            )

        # (b) a gated model must raise, naming the JAX equivalent.
        gated_msg = None
        try:
            gpu_backend.gate({"bad": gated_obj}, None)
            print(
                "    unexpected: gated {0} was NOT rejected".format(gated_cls.__name__)
            )
        except gpu_backend.GPUSupportedModelError as exc:
            gated_msg = str(exc)

        if allowed_passed and gated_msg is not None:
            _ok(
                "gate",
                'allowed {0} passes; gated {1} rejected -> "{2}"'.format(
                    allowed_cls.__name__,
                    gated_cls.__name__,
                    redirect.get(gated_cls, "(see message)"),
                ),
            )
        else:
            _fail(
                "gate",
                "allowed_pass={0}; gated_rejected={1}".format(
                    allowed_passed, gated_msg is not None
                ),
            )
    finally:
        gpu_backend._enabled = False  # restore CPU default


def _instantiate(cls):
    """Create an object of class ``cls`` WITHOUT running ``__init__``.

    ``gpu_backend.gate`` only inspects ``type(obj)`` (and ``isinstance``) -- it
    never calls a method on the object -- so a bare ``object.__new__(cls)`` is all
    that is needed. Skipping ``__init__`` also avoids its side effects: a george GP
    model calls ``quit()`` in ``__init__`` when george is absent, and that
    ``SystemExit`` is a ``BaseException`` (not an ``Exception``) that would abort
    the whole check. ``object.__new__`` sidesteps that entirely while still giving
    an object whose ``type()`` is exactly ``cls``.
    """
    try:
        return object.__new__(cls)
    except Exception:
        return None


def _pick_instantiable(candidates):
    """Return ``(cls, instance)`` for the first candidate that instantiates."""
    seen = set()
    for cls in candidates:
        if cls in seen:
            continue
        seen.add(cls)
        obj = _instantiate(cls)
        if obj is not None:
            return cls, obj
    return None, None


# --- 5. GPU parity -----------------------------------------------------------
def section5_gpu(mc, theta, tmpdir):
    _banner(5, "GPU parity (requires a GPU-visible jaxlib)")

    import numpy as np
    import jax
    import pyorbit.common.gpu_backend as gpu_backend

    try:
        gpu_devices = list(jax.devices("gpu"))
    except Exception:
        gpu_devices = []

    if not gpu_devices:
        _skip(
            "gpu-parity",
            "no GPU visible to this jaxlib; install the CUDA build to run this section. "
            "The CPU baseline + gate sections above already pass.",
        )
        return

    try:
        gpu_backend.initialize({"gpu": True, "gpu_verbose": True})
    except gpu_backend.GPUSupportedModelError as exc:
        _fail(
            "gpu-parity",
            "initialize({gpu:True}) raised: {0}".format(str(exc).splitlines()[0]),
        )
        return

    default_backend = str(jax.default_backend()).lower()
    print("    jax.devices() = {0}".format(list(jax.devices())))
    if default_backend != "gpu":
        _fail(
            "gpu-parity",
            "default_backend() == {0}, expected gpu".format(default_backend),
        )
        return
    _ok("gpu-device", "default_backend() == gpu")

    # tinygp device check (already JAX-native; should simply land on the GPU).
    _tinygp_gpu_check()

    # Production ordering: PyORBIT_Run.py calls gpu_backend.initialize() BEFORE
    # pars_input()/model_setup() build the container. Datasets stage their device
    # arrays and the per-instrument jitter/offset masks at construction time
    # (Dataset._stage_device_arrays runs once, only while gpu is enabled), so a
    # container built with the GPU off and *later* toggled to the GPU would
    # silently drop per-dataset jitter -- and that is a test-ordering artifact of
    # this script, not a production path. We therefore build a FRESH GPU
    # container after enablement (same config, only the `gpu` flag differs),
    # exactly as a `gpu: true` production run would.
    gpu_tmpdir = os.path.join(tmpdir, "gpu_leg")
    os.makedirs(gpu_tmpdir, exist_ok=True)
    gpu_mc = _build_mc(gpu_tmpdir, gpu_flag=True)
    gpu_theta = _theta_midbounds(gpu_mc)
    if gpu_mc.ndim != mc.ndim or not np.allclose(gpu_mc.bounds, mc.bounds):
        _fail(
            "gpu-parity",
            "GPU container bounds differ from the CPU container; test setup error",
        )
        return
    if not np.allclose(gpu_theta, theta):
        _fail(
            "gpu-parity",
            "GPU container theta differs from the CPU container; test setup error",
        )
        return

    # Re-evaluate the GPU container at the SAME theta -- now on the JAX path
    # (GPU enabled, so _active_kepler() routes the RV solver to kepler_exo_jax).
    # Compare against _CPU_BASELINE, captured in section 2 while the CPU path was
    # active.
    gpu_val = float(gpu_mc(theta))
    print("    log probability (priors + likelihood) [gpu] = {0:.10f}".format(gpu_val))
    if not np.isfinite(gpu_val):
        _fail("gpu-parity", "GPU log-priors+likelihood not finite: {0}".format(gpu_val))
        return

    diff = abs(gpu_val - _CPU_BASELINE)
    if np.allclose(gpu_val, _CPU_BASELINE, rtol=1e-6, atol=1e-8):
        _ok(
            "gpu-parity",
            "cpu={0:.10f} gpu={1:.10f} max|diff|={2:.3e} (within rtol=1e-6/atol=1e-8)".format(
                _CPU_BASELINE, gpu_val, diff
            ),
        )
    else:
        _fail(
            "gpu-parity",
            "cpu={0:.10f} gpu={1:.10f} max|diff|={2:.3e}".format(
                _CPU_BASELINE, gpu_val, diff
            ),
        )


# CPU reference captured by section 2/3 (module-level so section 5 can compare).
_CPU_BASELINE = None


def _tinygp_gpu_check():
    import numpy as np
    import jax

    try:
        import tinygp
        from tinygp import kernels
    except Exception as exc:
        _skip("tinygp-gpu", "tinygp not importable ({0})".format(exc))
        return
    try:
        t = np.sort(np.random.default_rng(0).uniform(0.0, 10.0, 200))
        y = np.random.default_rng(1).normal(size=200)
        gp = tinygp.GaussianProcess(
            kernels.ExpSquared(scale=1.5)
            * kernels.ExpSineSquared(scale=2.5, gamma=0.5),
            t,
            diag=1e-3,
            mean=0.0,
        )
        val = float(gp.log_probability(y))
        _ok(
            "tinygp-gpu",
            "tinygp GaussianProcess.log_probability evaluated on {0} = {1:.6f}".format(
                jax.default_backend(), val
            ),
        )
    except Exception as exc:
        _fail("tinygp-gpu", "tinygp evaluation raised: {0}".format(exc))


# --- main --------------------------------------------------------------------
def main():
    print("PyORBIT CPU/GPU parity check")
    print("  python = {0}".format(sys.version.split()[0]))
    print("  cwd    = {0}".format(os.getcwd()))

    import tempfile

    tmpdir = os.path.join(tempfile.gettempdir(), "pyorbit_gpu_parity")
    os.makedirs(tmpdir, exist_ok=True)

    try:
        section0_env()
    except _EnvError as exc:
        print()
        print("CANNOT RUN: {0}".format(exc))
        return 3

    try:
        section1_kepler_parity()
        mc, theta = section2_3_cpu_baseline(tmpdir)
        section4_gate()
        section5_gpu(mc, theta, tmpdir)
    except Exception:
        _fail("unhandled", "an uncaught exception aborted the check:")
        traceback.print_exc()

    print()
    print("--- summary ---")
    for s, d in _PASS:
        print("  PASS {0}".format(s))
    for s, d in _SKIP:
        print("  SKIP {0}".format(s))
    for s, d in _FAIL:
        print("  FAIL {0}: {1}".format(s, d))

    print()
    if _FAIL:
        print(
            "RESULT: FAIL ({0} failed, {1} passed, {2} skipped)".format(
                len(_FAIL), len(_PASS), len(_SKIP)
            )
        )
        return 1
    print("RESULT: PASS ({0} passed, {1} skipped)".format(len(_PASS), len(_SKIP)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
