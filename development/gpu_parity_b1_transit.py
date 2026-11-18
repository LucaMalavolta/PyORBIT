#!/usr/bin/env python
"""B1 transit-convention probe: GPU ``circular_transit`` vs batman CPU reference.

Motivation
==========
``batman_transit.Batman_Transit`` is a **unitary** model
(``abstract_transit`` sets ``unitary_model = True``) composed by
``common/dataset.py`` as::

    model = additive + (1 + unitary) * normalization

So the unitary term must be ~0 out of transit and dip (or rise) in transit.
The CPU branch of ``compute`` deliberately returns
``batman_model.light_curve(params) - 1``  (batman flux ~1 out of transit) to
produce that convention.  The GPU branch returns
``transit_gpu.circular_transit(...)`` directly, and
``transit_gpu.circular_transit`` returns the raw
``jaxoplanet.core.limb_dark.light_curve`` (also ~1 out of transit) **without**
subtracting 1.  If that is really true, the whole GPU photometric model is
shifted by +1.0 in relative flux versus the CPU (red flag R1).

This script settles it by execution: for a parameter grid it evaluates, in one
process,

  * GPU leg : ``pyorbit.subroutines.transit_gpu.circular_transit(...)``
  * CPU leg : ``batman.TransitModel(...).light_curve(params) - 1``

and reports the largest |GPU - (CPU-1)| per configuration.

Classification (fixed tolerances, never relaxed)
  * max|diff| <= 1e-3          -> PASS        (legitimate quadrature floor)
  * 1e-3 <  max|diff| < 0.5    -> SUSPICIOUS  (stop, report, do not fix)
  * max|diff| >= 0.5 (≈O(1))   -> R1 CONFIRMED

Runs in the ``pyorbit-gpu`` conda env (only env with both
``jaxoplanet`` and ``batman``).

Exit codes: 0 = PASS, 2 = SUSPICIOUS, 1 = R1 CONFIRMED (or any other failure).
"""

import sys

from itertools import product

# --- fixed grid (task-specified) ---------------------------------------------
PERIODS = (5.0, 10.0, 15.0)  # d
TC_MULTIPLIERS = (-2, -1, -0.5, -0.1, 0, 0.1, 0.5, 1, 2)  # x P/100 in days
R_RS = (0.02, 0.1, 0.3)
A_RS = (5.0, 15.0, 30.0)
INCL = (82.0, 85.0, 89.0)  # deg
U = (0.1, 0.3)  # quadratic limb darkening
SAMPLE_FACTORS = (1, 5, 10)
EXP_TIMES_S = (0, 60, 120)  # s (-> days like code_options)
N_TIMES = 251
PASS_MAX = 1e-3
SUSPICIOUS_MAX = 0.5


def _grid_desc():
    return (
        "grid: P={P}; (Tc-Tref) = {M} x P/100; R_Rs={RR}; a_Rs={AR}; i={I} deg; "
        "u=({u0},{u1}); sample_factor={SF}; exp_time={ET} s; "
        "times = {N} pts spanning one full period; "
        "{K} configurations"
    ).format(
        P=PERIODS,
        M=TC_MULTIPLIERS,
        RR=R_RS,
        AR=A_RS,
        I=INCL,
        u0=U[0],
        u1=U[1],
        SF=SAMPLE_FACTORS,
        ET=EXP_TIMES_S,
        N=N_TIMES,
        K=len(PERIODS)
        * len(TC_MULTIPLIERS)
        * len(R_RS)
        * len(A_RS)
        * len(INCL)
        * len(SAMPLE_FACTORS)
        * len(EXP_TIMES_S),
    )


def main():
    import numpy as np

    np.seterr(all="ignore")

    print("B1 transit-convention probe")
    print("  " + _grid_desc())

    # Enable float64 (the same way PyORBIT does with gpu: true).
    import jax

    try:
        jax.config.update("jax_enable_x64", True)
    except Exception:
        pass

    if "gpu" not in str(jax.default_backend()).lower():
        print(
            "RUNNING ON BACKEND {0}: a GPU is expected (pyorbit-gpu env). "
            "The convention question is backend-independent, so the numbers "
            "are still valid.".format(jax.default_backend())
        )

    from pyorbit.subroutines.transit_gpu import circular_transit
    import batman

    results = []  # (max_abs, tag) pairs

    model_cache = {}
    zero_exp_approx = 0  # configs using the batman exp_time floor (see below)
    for P in PERIODS:
        tc_offsets = tuple(float(m) * P / 100.0 for m in TC_MULTIPLIERS)
        # One shared time grid (days, Tref = 0), one full period: covers both
        # in-transit and out-of-transit times for every tc offset on this P.
        times = np.linspace(0.0, P, N_TIMES)
        base_params = batman.TransitParams()
        base_params.t0 = 0.0
        base_params.per = P
        base_params.rp = 0.1
        base_params.a = 10.0
        base_params.inc = 87.0
        base_params.ecc = 0.0
        base_params.w = 90.0
        base_params.limb_dark = "quadratic"
        base_params.u = np.array(U)

        for R, A, inc, sf, et_s in product(
            R_RS, A_RS, INCL, SAMPLE_FACTORS, EXP_TIMES_S
        ):
            key = (P, A, inc, sf, et_s)
            if key not in model_cache:
                # batman requires exp_time > 0 whenever supersampling is on.
                # With exp_time == 0 both legs reduce to point sampling
                # (GPU offsets are all zero -> mean of identical points),
                # so build the batman model with supersample_factor=1 there.
                sf_batman = sf if et_s > 0 else 1
                p = batman.TransitParams()
                p.per = P
                p.a = A
                p.inc = inc
                p.ecc = 0.0
                p.w = 90.0
                p.limb_dark = "quadratic"
                p.u = np.array(U)
                p.t0 = 0.0
                p.rp = 0.1
                model_cache[key] = batman.TransitModel(
                    p,
                    times,
                    supersample_factor=sf_batman,
                    exp_time=et_s / 86400.0,
                    nthreads=1,
                )
                if sf > 1 and et_s == 0:
                    zero_exp_approx += 1
            model = model_cache[key]

            for tc in tc_offsets:
                params = batman.TransitParams()
                params.t0 = tc
                params.per = P
                params.rp = R
                params.a = A
                params.inc = inc
                params.ecc = 0.0
                params.w = 90.0
                params.limb_dark = "quadratic"
                params.u = np.array(U)

                cpu = model.light_curve(params) - 1.0  # CPU convention
                gpu = circular_transit(
                    jax.numpy.asarray(times),
                    P,
                    tc,
                    R,
                    A,
                    inc * np.pi / 180.0,
                    jax.numpy.asarray(np.array(U)),
                    et_s / 86400.0,
                    samples=sf,
                )
                gpu = np.asarray(gpu, dtype=float)

                if cpu.shape != gpu.shape:
                    print("  shape mismatch {0} vs {1}".format(cpu.shape, gpu.shape))
                    return 1
                diff = np.abs(gpu - cpu)
                if np.any(~np.isfinite(diff)):
                    print(
                        "  NON-FINITE for {0}".format(
                            dict(P=P, tc=tc, R=R, A=A, inc=inc, sf=sf, et_s=et_s)
                        )
                    )
                    return 1
                m = float(np.max(diff))
                results.append(
                    (
                        m,
                        "P={0:5.1f} tc={1:+.4f} R_Rs={2:4.2f} "
                        "a_Rs={3:5.1f} i={4:4.1f} sf={5:2d} "
                        "et={6:3d}s max|diff|={7:.6e}".format(
                            P, tc, R, A, inc, sf, et_s, m
                        ),
                    )
                )

    results.sort(reverse=True, key=lambda r: r[0])

    print()
    if zero_exp_approx:
        print(
            "note: batman cannot supersample at exp_time=0, so {0} "
            "model-groups (sf>1, et=0 s) used point sampling on the CPU "
            "leg; the GPU leg degenerates to the same point there".format(
                zero_exp_approx
            )
        )
    print()
    print("--- worst 10 configurations (full table in log) ---")
    for m, tag in results[:10]:
        print("  " + tag)
    for m, tag in results[10:]:
        print("  " + tag)

    overall = results[0][0]
    print()
    if overall <= PASS_MAX:
        label = (
            "PASS (max|diff|={0:.6e} <= {1:.0e}: legitimate "
            "quadrature floor; no convention bug)".format(overall, PASS_MAX)
        )
        code = 0
    elif overall < SUSPICIOUS_MAX:
        label = (
            "SUSPICIOUS (max|diff|={0:.6e} in ({1:.0e}, {2:.0e}): "
            "stop and report, do not fix)".format(overall, PASS_MAX, SUSPICIOUS_MAX)
        )
        code = 2
    else:
        label = (
            "R1 CONFIRMED (max|diff|={0:.6e} >= {1:.0e}: GPU branch "
            "returns the raw light curve without the '-1' convention "
            "used by the CPU branch)".format(overall, SUSPICIOUS_MAX)
        )
        code = 1
    print("RESULT: {0}".format(label))
    return code


if __name__ == "__main__":
    sys.exit(main())
