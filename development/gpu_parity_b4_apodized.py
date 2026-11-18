#!/usr/bin/env python
"""B4 apodized-RV parity probe.

Motivation
==========
``pyorbit/models/keplerian_modelling.py`` dispatches RV + apodization through
``_active_kepler()``:

* ``gpu: true``  -> ``jax.scipy.stats.norm.pdf`` for the apodizing term and the
  JAX ``kepler_exo_jax`` solver.
* ``gpu: false`` -> ``scipy.stats.norm.pdf`` and the numpy ``kepler_exo`` solver.

``ApodizedRVkeplerian.compute`` builds the output as::

    apo = norm.pdf(x, center, ts) / norm.pdf(center, center, ts)
    return apo * rv_deltabjd(x0, K, P, mean_long, e, omega, Omega)

This probe wires the production function **both ways** (``gpu_backend._enabled``
toggled, exactly like the parity harness) and also compares the ``norm_pdf``
term in isolation. Both legs are float64 under ``jax_enable_x64``.

The production CPU/RV parity (Section 1 of ``gpu_cpu_parity.py``) already bounds
the ``rv_deltabjd`` term to ~1e-12; B4 isolates the apodizing ``norm_pdf`` term
that the GPU path changes.

Classification (fixed tolerances, never relaxed)
  PASS        : max rel diff (full model) <= 1e-10 AND max|diff| (norm_pdf term)
                on the scale of the term is bounded by rel <= 1e-10
  SUSPICIOUS  : rel within (1e-10, 1e-6)   -> STOP, report, do not fix
  BUG         : rel >= 1e-6                -> STOP, report, do not fix

Runs in the ``pyorbit-gpu`` conda env (has jax CUDA + scipy).

Exit codes: 0 = PASS, 2 = SUSPICIOUS, 1 = BUG / any other failure.
"""

import sys


# --- fixed grid (task-specified) --------------------------------------------
REF = 78000.0  # Tref, day-like
TSPAN = 200.0  # day-like half-span for x
N = 400

P_GRID = (5.0, 20.0, 40.0)  # orbital period
E_GRID = (0.0, 0.2, 0.5)  # eccentricity
K_GRID = (1.0, 10.0)  # semi-amplitude (m/s)
CENTER_OFF = (-30.0, 0.0, 30.0)  # apo_center relative to Tref
TS_GRID = (5.0, 30.0, 200.0)  # apo_timescale
PASS_REL = 1e-10
SUS_REL = 1e-6


def _make_obj():
    """A minimal ApodizedRVkeplerian with just the attributes ``compute`` reads."""
    from pyorbit.models.keplerian_modelling import ApodizedRVkeplerian

    obj = object.__new__(ApodizedRVkeplerian)  # skip __init__
    obj.use_time_inferior_conjunction = False
    obj.use_mass = False
    obj.use_scaled_mass = False
    obj.use_stellar_scaled_mass = False
    obj.default_Omega = 0.0
    return obj, ApodizedRVkeplerian


def main():
    import numpy as np
    from types import SimpleNamespace

    np.random.seed(0)

    print("B4 apodized-RV parity probe")
    print(
        "  N={0}, P={1}, e={2}, K={3}, center_off={4}, ts={5}".format(
            N, P_GRID, E_GRID, K_GRID, CENTER_OFF, TS_GRID
        )
    )

    import jax

    try:
        jax.config.update("jax_enable_x64", True)
    except Exception:
        pass
    print(
        "  jax={0} backend(after x64)={1}".format(
            jax.__version__, jax.default_backend()
        )
    )

    from scipy.stats import norm as scipy_norm
    import pyorbit.common.gpu_backend as gpu_backend

    results = []  # (metric, diff, rel, tag)

    # =========================================================================
    # (1) norm_pdf term, in isolation
    # =========================================================================
    print("\n-- (1) norm_pdf term: jax.scipy.stats.norm.pdf vs scipy.stats.norm.pdf --")
    x = np.linspace(REF - TSPAN, REF + TSPAN, N)
    for coff in CENTER_OFF:
        center = REF + coff
        for ts in TS_GRID:
            gpu = np.asarray(
                jax.scipy.stats.norm.pdf(jax.numpy.asarray(x), center, ts), dtype=float
            )
            cpu = np.asarray(scipy_norm.pdf(x, center, ts), dtype=float)
            diff = np.abs(gpu - cpu)
            mag = np.abs(cpu)
            with np.errstate(invalid="ignore", divide="ignore"):
                rel = np.where(mag > 0, diff / mag, 0.0)
            mrel = float(np.max(rel))
            tag = "center={0:+.1f} {1:.1f}".format(ts, coff)
            results.append(("norm_pdf", float(np.max(diff)), mrel, tag))
            print(
                "  {0:22s} max|diff|={1:.6e} max rel={2:.6e}".format(
                    tag, np.max(diff), mrel
                )
            )

    # =========================================================================
    # (2) full ApodizedRVkeplerian.compute output, both ways
    # =========================================================================
    print(
        "\n-- (2) full model: ApodizedRVkeplerian.compute both ways (flag-toggled) --"
    )
    obj, ApodizedRVkeplerian = _make_obj()
    # phase-referenced times for the RV solve; x for the apodizing term
    x0 = x - REF
    dataset = SimpleNamespace(x=x, x0=x0, Tref=REF)

    rng = np.random.default_rng(2026)
    for P in P_GRID:
        for e in E_GRID:
            for K in K_GRID:
                omega = float(rng.uniform(0.0, 2.0 * np.pi))
                mean_long = float(rng.uniform(0.0, 2.0 * np.pi))
                for coff in CENTER_OFF:
                    center = REF + coff
                    for ts in TS_GRID:
                        pv = {
                            "P": P,
                            "K": K,
                            "e": e,
                            "omega": omega,
                            "mean_long": mean_long,
                            "apo_center": center,
                            "apo_timescale": ts,
                        }
                        gpu_backend._enabled = False
                        cpu_full = np.asarray(
                            ApodizedRVkeplerian.compute(obj, dict(pv), dataset),
                            dtype=float,
                        )
                        gpu_backend._enabled = True
                        try:
                            gpu_full = np.asarray(
                                ApodizedRVkeplerian.compute(obj, dict(pv), dataset),
                                dtype=float,
                            )
                        finally:
                            gpu_backend._enabled = False

                        if not (
                            np.all(np.isfinite(cpu_full))
                            and np.all(np.isfinite(gpu_full))
                        ):
                            print(
                                "  NON-FINITE full model P={0} e={1} K={2} c={3} ts={4}".format(
                                    P, e, K, center, ts
                                )
                            )
                            return 1
                        diff = np.abs(gpu_full - cpu_full)
                        mag = np.abs(cpu_full)
                        with np.errstate(invalid="ignore", divide="ignore"):
                            rel = np.where(mag > 0, diff / mag, 0.0)
                        mrel = float(np.max(rel))
                        tag = "P={0:4.1f} e={1:.2f} K={2:4.1f} c={3:+.1f} ts={4:3.1f}".format(
                            P, e, K, coff, ts
                        )
                        results.append(("full", float(np.max(diff)), mrel, tag))
                        print(
                            "  {0:44s} max|diff|={1:.6e} max rel={2:.6e}".format(
                                tag, np.max(diff), mrel
                            )
                        )

    # --- verdict -----------------------------------------------------------
    # (full model rel is the governing metric for B4's "full model output both ways";
    #  the norm_pdf-term rel is reported alongside.)
    full_rels = [r[2] for r in results if r[0] == "full"]
    worst_full = max(full_rels)
    worst_norm = max(r[2] for r in results if r[0] == "norm_pdf")
    # also capture absolute max diff for the report
    worst_full_abs = max(r[1] for r in results if r[0] == "full")
    worst_norm_abs = max(r[1] for r in results if r[0] == "norm_pdf")

    print()
    print(
        "worst full  : max rel={0:.6e}  (abs max|diff|={1:.6e})".format(
            worst_full, worst_full_abs
        )
    )
    print(
        "worst norm  : max rel={0:.6e}  (abs max|diff|={1:.6e})".format(
            worst_norm, worst_norm_abs
        )
    )

    if worst_full <= PASS_REL and worst_norm <= PASS_REL:
        print(
            "RESULT: PASS (full max rel={0:.6e} and norm max rel={1:.6e} <= {2:.0e}: "
            "jax vs scipy norm.pdf and the full apodized RV agree at/near float64).".format(
                worst_full, worst_norm, PASS_REL
            )
        )
        return 0
    if worst_full < SUS_REL and worst_norm < SUS_REL:
        print(
            "RESULT: SUSPICIOUS (worst full rel={0:.6e} in ({1:.0e},{2:.0e}). "
            "STOP; capture; do not fix.".format(worst_full, PASS_REL, SUS_REL)
        )
        return 2
    print(
        "RESULT: BUG (worst full rel={0:.6e} >= {1:.0e} OR norm rel out of band). "
        "STOP; capture the numbers above; do not fix.".format(worst_full, SUS_REL)
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
