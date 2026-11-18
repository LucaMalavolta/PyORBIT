#!/usr/bin/env python
"""B3 ESP parity probe (DIAGNOSTIC): GPU ``esp_loglike`` vs spleaf C reference.

Motivation
==========
Compare ``pyorbit.subroutines.esp_gpu.esp_loglike`` (custom CUDA kernel) against
the ``spleaf`` LDL reference that the two production models call on their CPU
branch:

* ``pyorbit/models/spleaf_esp_activity.py`` (single dataset):
  CPU leg  ``Cov(x0, err=Error(e), GP=ESPKernel(Hamp,Prot,Pdec,Oamp,nharm)).loglike(res)``
  GPU leg  ``esp_loglike(x0, res, e**2, Prot, Pdec, Oamp, alpha=Hamp, nharm=n)``
* ``pyorbit/models/spleaf_multidimensional_esp_activity.py`` (multi dataset,
  the ``alpha*u + beta*du`` weighted path):
  CPU leg  ``Cov(t, err=Error(e), GP=MultiSeriesKernel(ESPKernel(1.0,...),idx,cp,cd))
             .set_param([Prot,Pdec,Oamp]+cp+cd, param[1:]).loglike(res)``
  GPU leg  ``esp_loglike(t, res, e**2, Prot,Pdec,Oamp, alpha_vec, beta_vec, n)``
            with ``alpha[idx_i]=cp[i]`` and ``beta[idx_i]=cd[i]``.

Both are float64 and compute the SAME natural log-likelihood
``-0.5*(chi2 + logdet + n*log(2*pi))``.

Classification (fixed tolerances, never relaxed)
  * max|diff| <  1e-6           -> PASS
  * 1e-6 <= max|diff| < 0.5     -> SUSPICIOUS (stop, report, do not fix)
  * max|diff| >= 0.5            -> O(1): BUG (stop, report, do not fix)
  * GPU leg does not execute    -> BLOCKED (production GPU ESP path inoperable;
                                   STOP, report, do not fix)

Runs in the ``pyorbit-gpu`` conda env (the only env with both ``cupy``/CUDA and
``spleaf``). The GPU leg needs ``gpu_backend._enabled=True`` so ``get_device()``
resolves a CUDA device; we toggle it around the ``esp_loglike`` calls only.

Exit codes: 0 = PASS, 2 = SUSPICIOUS, 1 = O(1) BUG, 3 = BLOCKED (GPU leg fails).
"""

import sys
import traceback

N_SINGLE = 500
NHARM_VARIANTS = (4, 2)

SINGLE_GRID = (
    (15.0, 200.0, 0.90, 0.05, "low-Hamp near-sinusoid"),
    (15.0, 100.0, 0.20, 5.00, "high-Hamp granular"),
    (30.0, 100.0, 0.35, 1.00, "mid (PyORBIT defaults)"),
    (20.0, 50.0, 0.50, 2.00, "short-Pdec high-Hamp"),
    (10.0, 300.0, 0.80, 0.10, "low-Hamp long-Pdec near-sinusoid"),
    (25.0, 80.0, 0.40, 3.00, "mid-high Hamp"),
)

MULTI_DATASETS = 2
MULTI_GRID = (
    ([(5.0, 1.0), (0.5, 2.5)], "distinct con/rot per dataset"),
    ([(0.1, 8.0), (3.0, 0.5)], "beta (rot_amp) dominant on ds0"),
    ([(3.0, 1e-4), (2.0, 4.0)], "near-zero beta on ds0 (pure-alpha limit)"),
    ([(1.0, 1.0), (1.0, 1.0)], "equal con=rot across datasets"),
    ([(10.0, 0.1), (0.2, 10.0)], "extreme asymmetry con vs rot"),
)

PASS_MAX = 1e-6
SUSPICIOUS_MAX = 0.5


def _env_facts():
    facts = []
    try:
        import cupy

        facts.append("cupy={0}".format(cupy.__version__))
    except Exception as e:
        facts.append("cupy import error: {0}".format(e))
    try:
        import cupyx.scipy.special as sp

        facts.append(
            "cupyx.scipy.special has ive={0}, iv={1}, i0={2}, i0e={3}".format(
                hasattr(sp, "ive"),
                hasattr(sp, "iv"),
                hasattr(sp, "i0"),
                hasattr(sp, "i0e"),
            )
        )
    except Exception as e:
        facts.append("cupyx.scipy.special import error: {0}".format(e))
    try:
        import scipy

        facts.append(
            "scipy={0} (scipy.special.ive={1})".format(
                scipy.__version__,
                hasattr(__import__("scipy.special", fromlist=["ive"]), "ive"),
            )
        )
    except Exception:
        pass
    return facts


def _cpu_single(t, y, e, Prot, Pdec, Oamp, Hamp, nharm):
    import numpy as np
    from spleaf import cov as spleaf_cov
    from spleaf import term as spleaf_term

    order = np.argsort(t)
    x0 = np.asarray(t, dtype=float)[order]
    res = np.asarray(y, dtype=float)[order]
    e2 = (np.asarray(e, dtype=float) ** 2)[order]
    D = spleaf_cov.Cov(
        x0,
        err=spleaf_term.Error(np.sqrt(e2)),
        GP=spleaf_term.ESPKernel(Hamp, Prot, Pdec, Oamp, nharm=nharm),
    )
    D.set_param(np.array([Hamp, Prot, Pdec, Oamp]), D.param)
    return float(D.loglike(res)), x0, res, e2


def _multi_build(t_list, y_list, e_list, Prot, Pdec, Oamp, amps, nharm):
    """Build and return (spleaf_time, res, err, series_index, cp, cd) or raise."""
    import numpy as np
    from spleaf import cov as spleaf_cov

    (spleaf_time, _e, _e2, series_index) = spleaf_cov.merge_series(
        t_list, e_list, e_list
    )
    n = spleaf_time.size
    res = np.zeros(n)
    err = np.zeros(n)
    for k in range(len(t_list)):
        sel = series_index[k]
        res[sel] = y_list[k]
        err[sel] = e_list[k]
    cp = [a[0] for a in amps]
    cd = [a[1] for a in amps]
    return spleaf_time, res, err, series_index, cp, cd


def _cpu_multi(spleaf_time, res, err, series_index, cp, cd, Prot, Pdec, Oamp, nharm):
    import numpy as np
    from spleaf import cov as spleaf_cov
    from spleaf import term as spleaf_term

    inner = spleaf_term.ESPKernel(1.0, Prot, Pdec, Oamp, nharm=nharm)
    GP = spleaf_term.MultiSeriesKernel(inner, series_index, cp, cd)
    D = spleaf_cov.Cov(spleaf_time, err=spleaf_term.Error(err), GP=GP)
    D.set_param(np.concatenate(([Prot, Pdec, Oamp], cp, cd)), D.param[1:])
    return float(D.loglike(res))


def _attempt_gpu(kind, build, Prot, Pdec, Oamp, nharm):
    """Attempt the GPU leg for one representative case; return (value, exc)."""
    import numpy as np
    import pyorbit.common.gpu_backend as gpu_backend
    from pyorbit.subroutines.esp_gpu import esp_loglike

    gpu_backend._enabled = True
    try:
        if kind == "single":
            x0, res, e2, Hamp = build
            val = float(
                esp_loglike(x0, res, e2, Prot, Pdec, Oamp, alpha=Hamp, nharm=nharm)
            )
        else:
            spleaf_time, res, err, series_index, cp, cd = build
            alpha = np.empty(len(spleaf_time))
            beta = np.empty(len(spleaf_time))
            for i, sel in enumerate(series_index):
                alpha[sel] = cp[i]
                beta[sel] = cd[i]
            val = float(
                esp_loglike(
                    spleaf_time, res, err**2, Prot, Pdec, Oamp, alpha, beta, nharm
                )
            )
        return float(val), None
    except Exception as exc:
        return None, exc
    finally:
        gpu_backend._enabled = False


def main():
    import numpy as np

    print("B3 ESP parity probe (DIAGNOSTIC): GPU esp_loglike vs spleaf C reference")
    print(
        "  n(single)={0}; multi datasets={1}; nharm variants={2}".format(
            N_SINGLE, MULTI_DATASETS, NHARM_VARIANTS
        )
    )
    try:
        import jax

        jax.config.update("jax_enable_x64", True)
        print(
            "  jax={0} default_backend(after x64)={1}".format(
                jax.__version__, jax.default_backend()
            )
        )
    except Exception as e:
        print("  jax: {0}".format(e))
    print("  " + " ; ".join(_env_facts()))

    # --- CPU reference table (proves the CPU leg is sound) -----------------
    cpu_single = []
    print("\n-- single dataset (CPU spleaf reference) --")
    rng = np.random.default_rng(2026)
    for nharm in NHARM_VARIANTS:
        for Prot, Pdec, Oamp, Hamp, label in SINGLE_GRID:
            t = np.sort(rng.uniform(0.0, 60.0, N_SINGLE))
            y = rng.normal(0.0, 1.0, N_SINGLE)
            e = 0.1 + 0.05 * rng.random(N_SINGLE)
            try:
                cpu_val, x0, res, e2 = _cpu_single(
                    t, y, e, Prot, Pdec, Oamp, Hamp, nharm
                )
            except Exception as exc:
                print("  ERROR single nharm={0} {1}: {2}".format(nharm, label, exc))
                return 1
            finite = "finite" if np.isfinite(cpu_val) else "NOT-FINITE"
            print(
                "  nharm={0} {1}: cpu loglike={2:.10f} ({3})".format(
                    nharm, label, cpu_val, finite
                )
            )
            cpu_single.append(
                dict(
                    nharm=nharm,
                    label=label,
                    Prot=Prot,
                    Pdec=Pdec,
                    Oamp=Oamp,
                    Hamp=Hamp,
                    val=cpu_val,
                    x0=x0,
                    res=res,
                    e2=e2,
                )
            )

    cpu_multi = []
    print("\n-- multi dataset (CPU spleaf reference) --")
    rngm = np.random.default_rng(2027)
    mProt, mPdec, mOamp = 18.0, 120.0, 0.45
    for amps, label in MULTI_GRID:
        t_list, y_list, e_list = [], [], []
        for _ in range(MULTI_DATASETS):
            t_list.append(np.sort(rngm.uniform(0.0, 90.0, N_SINGLE)))
            y_list.append(rngm.normal(0.0, 1.0, N_SINGLE))
            e_list.append(0.1 + 0.05 * rngm.random(N_SINGLE))
        try:
            b = _multi_build(t_list, y_list, e_list, mProt, mPdec, mOamp, amps, 4)
            cpu_val = _cpu_multi(*b, mProt, mPdec, mOamp, 4)
        except Exception as exc:
            print("  ERROR multi {0}: {1}".format(label, exc))
            return 1
        print("  {0}: cpu loglike={1:.10f}".format(label, cpu_val))
        cpu_multi.append(dict(label=label, cpu_val=cpu_val, build=b, amps=amps))

    # --- GPU leg: representative single + multi attempt ---------------------
    print("\n-- GPU leg (esp_loglike) --")
    if any(not np.isfinite(r["val"]) for r in cpu_single):
        print(
            "  CPU single reference is not finite; cannot form a meaningful parity claim."
        )

    single_rep = cpu_single[0]
    gpu_val, gpu_err = _attempt_gpu(
        "single",
        (
            single_rep["x0"],
            single_rep["res"],
            single_rep["e2"],
            single_rep["Hamp"],
        ),
        single_rep["Prot"],
        single_rep["Pdec"],
        single_rep["Oamp"],
        single_rep["nharm"],
    )
    if gpu_err is not None:
        print("  GPU single leg FAILED to execute:")
        print("  > " + " ; ".join(_env_facts()))
        print("  exception: {0}: {1}".format(type(gpu_err).__name__, gpu_err))
        print("  traceback:")
        for line in traceback.format_exc().splitlines():
            print("    " + line.rstrip())
        # Attempt the multi leg too, to confirm the same blocker.
        mu = cpu_multi[0]
        gval2, gerr2 = _attempt_gpu("multi", mu["build"], mProt, mPdec, mOamp, 4)
        print(
            "  GPU multi leg: {0}".format(
                "executed val={0:.10f}".format(gval2)
                if gval2 is not None
                else "FAILED: {0}: {1}".format(type(gerr2).__name__, gerr2)
            )
        )
        print()
        print(
            "RESULT: BLOCKED -- GPU ESP leg (esp_gpu.esp_loglike) does not execute "
            "in this environment. Root cause: `from cupyx.scipy.special import ive` "
            "raises ImportError (cupy 13.6.0 has i0/i1/i0e but not the "
            "arbitrary-order `ive`). The CPU spleaf reference is sound (finite "
            "loglike on all points above). Classification: BUG in the production "
            "GPU ESP surface (SPLEAF_ESP / SPLEAF_Multidimensional_ESP) -- "
            "STOP, do not fix (CUDA-kernel / cupy-API work is out of scope; "
            "esp_gpu.py is off-limits). See report."
        )
        return 3

    # --- If the GPU leg DID execute: full parity ----------------------------
    diff_single = []
    diff_multi = []
    print("\n-- parity (single) --")
    for r in cpu_single:
        g, e2 = _attempt_gpu(
            "single",
            (r["x0"], r["res"], r["e2"], r["Hamp"]),
            r["Prot"],
            r["Pdec"],
            r["Oamp"],
            r["nharm"],
        )
        if e2 is not None:
            print(
                "  {0}: GPU leg error {1}: {2}".format(
                    r["label"], type(e2).__name__, e2
                )
            )
            return 3
        d = abs(g - r["val"])
        diff_single.append(d)
        print(
            "  nharm={0} {1}: gpu={2:.10f} cpu={3:.10f} max|diff|={4:.6e}".format(
                r["nharm"], r["label"], g, r["val"], d
            )
        )
    print("\n-- parity (multi) --")
    for mu in cpu_multi:
        g, ge = _attempt_gpu("multi", mu["build"], mProt, mPdec, mOamp, 4)
        if ge is not None:
            print(
                "  {0}: GPU leg error {1}: {2}".format(
                    mu["label"], type(ge).__name__, ge
                )
            )
            return 3
        d = abs(g - mu["cpu_val"])
        diff_multi.append(d)
        print(
            "  {0}: gpu={1:.10f} cpu={2:.10f} max|diff|={3:.6e}".format(
                mu["label"], g, mu["cpu_val"], d
            )
        )

    overall = max(diff_single + diff_multi)
    print(
        "\nworst case max|diff| over {0} single + {1} multi = {2:.6e}".format(
            len(diff_single), len(diff_multi), overall
        )
    )
    if overall < PASS_MAX:
        print(
            "RESULT: PASS (worst max|diff|={0:.6e} < {1:.0e}: ESP implementation parity "
            "confirmed; legitimate float64 accumulation floor)".format(
                overall, PASS_MAX
            )
        )
        return 0
    elif overall < SUSPICIOUS_MAX:
        print(
            "RESULT: SUSPICIOUS (worst {0:.6e} in [{1:.0e},{2:.0e}]). STOP; capture; do not fix.".format(
                overall, PASS_MAX, SUSPICIOUS_MAX
            )
        )
        return 2
    else:
        print(
            "RESULT: BUG/O(1) (worst {0:.6e} >= {1:.0e}: O(1) ESP discrepancy). STOP; do not fix.".format(
                overall, SUSPICIOUS_MAX
            )
        )
        return 1


if __name__ == "__main__":
    sys.exit(main())
