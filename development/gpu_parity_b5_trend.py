#!/usr/bin/env python
"""B5 polynomial-trend parity probe.

Motivation
==========
``pyorbit/models/polynomial_trend.py :: PolynomialTrend.compute`` dispatches
on the GPU flag:

* ``gpu: true``  -- Horner evaluation on ``jax.numpy`` with the production
  coefficient array ``c[i] * x**i`` (increasing order), optionally using the
  staged ``dataset.x0_dev`` device fast path.
* ``gpu: false`` -- ``numpy.polynomial.polynomial.polyval`` (the production
  CPU reference; ``np.polyval`` requires the coefficient order reversed and
  evaluates the identical polynomial).

This probe wires the production function **both ways** (``gpu_backend._enabled``
toggled, exactly like the B4 parity harness), and covers **both** production
wiring variants of the GPU branch: the staged-``x0_dev`` fast path and the
plain host-``x_in`` path. Both legs are float64 under ``jax_enable_x64``.

Grid (task-specified)
  orders 1..5; coefficient sets spanning 1e-4..1e4 incl. alternating signs;
  x over [-3, 3] and [0, 1], 100 points each.

Classification (fixed tolerances, never relaxed)
  PASS        : max rel diff <= 1e-10
  SUSPICIOUS  : 1e-10 < rel < 1e-6   -> STOP, report, do not fix
  BUG         : rel >= 1e-6          -> STOP, report, do not fix

Runs in the ``pyorbit-gpu`` conda env (has jax CUDA).

Exit codes: 0 = PASS, 2 = SUSPICIOUS, 1 = BUG / any other failure.
"""

import sys

# --- fixed grid (task-specified) --------------------------------------------
TREF = 0.0  # dataset.Tref; x_zero = TREF so x_in == x on the probe grids
N = 100
GRID = ((-3.0, 3.0), (0.0, 1.0))
ORDERS = (1, 2, 3, 4, 5)
MAGS = (-4, -1, 2, 4)  # exponents spanning 1e-4 .. 1e4 (via the pos set)
PASS_REL = 1e-10
SUS_REL = 1e-6


def _coeff_sets(order):
    """Coefficient sets c1..c{order} spanning 1e-4..1e4 incl. alternating signs."""
    pos = [10.0 ** MAGS[(i - 1) % 4] for i in range(1, order + 1)]
    alt = [(-1.0) ** (i - 1) * 10.0 ** MAGS[(i - 1) % 4] for i in range(1, order + 1)]
    desc = [(-1.0) ** (i - 1) * 10.0 ** (5 - i) for i in range(1, order + 1)]
    return (("pos", pos), ("alt", alt), ("desc", desc))


def _make_obj(order):
    """A minimal PolynomialTrend with just the attributes ``compute`` reads."""
    from pyorbit.models.polynomial_trend import PolynomialTrend

    obj = object.__new__(PolynomialTrend)  # skip AbstractModel.__init__
    obj.order = order
    obj.starting_order = 1
    obj.time_interval = 1.0
    return obj


def main():
    import numpy as np
    from types import SimpleNamespace

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

    print("B5 polynomial-trend parity probe")
    print("  orders={0}, x-grids={1}, N={2}, Tref={3}".format(ORDERS, GRID, N, TREF))
    print("  GPU leg: PolynomialTrend.compute (gpu flag on, both wiring variants)")
    print("  CPU leg: PolynomialTrend.compute (gpu flag off -> polynomial.polyval)")

    import jax.numpy as jnp
    from pyorbit.models.polynomial_trend import PolynomialTrend
    import pyorbit.common.gpu_backend as gpu_backend

    results = []  # (abs_max_diff, max_rel, tag)
    ncases = 0

    for order in ORDERS:
        obj = _make_obj(order)
        for name, cs in _coeff_sets(order):
            pv = {"x_zero": TREF}
            for i, c in enumerate(cs, start=1):
                pv["poly_c" + repr(i)] = c
            for a, b in GRID:
                x = np.linspace(a, b, N)
                for variant in ("x0_dev", "x_in"):
                    ds = SimpleNamespace(x=x, Tref=TREF)
                    if variant == "x0_dev":
                        ds.x0_dev = jnp.asarray(x - TREF)
                    else:
                        ds.x0_dev = None

                    gpu_backend._enabled = False
                    cpu = np.asarray(
                        PolynomialTrend.compute(obj, dict(pv), ds), dtype=float
                    )
                    try:
                        gpu_backend._enabled = True
                        gpu = np.asarray(
                            PolynomialTrend.compute(obj, dict(pv), ds), dtype=float
                        )
                    finally:
                        gpu_backend._enabled = False

                    ncases += 1
                    if not (np.all(np.isfinite(cpu)) and np.all(np.isfinite(gpu))):
                        print(
                            "  NON-FINITE order={0} set={1} grid=[{2},{3}] {4}".format(
                                order, name, a, b, variant
                            )
                        )
                        return 1
                    diff = np.abs(gpu - cpu)
                    mag = np.abs(cpu)
                    with np.errstate(invalid="ignore", divide="ignore"):
                        rel = np.where(mag > 0, diff / mag, 0.0)
                    mrel = float(np.max(rel))
                    tag = "order={0} set={1:4s} x=[{2:5.1f},{3:5.1f}] {4}".format(
                        order, name, a, b, variant
                    )
                    results.append((float(np.max(diff)), mrel, tag))
                    print(
                        "  {0:38s} max|diff|={1:.6e} max rel={2:.6e}".format(
                            tag, np.max(diff), mrel
                        )
                    )

    worst = max(r[1] for r in results)
    worst_abs = max(r[0] for r in results)
    print()
    print(
        "ncases={0} worst: max rel={1:.6e}  (abs max|diff|={2:.6e})".format(
            ncases, worst, worst_abs
        )
    )

    if worst <= PASS_REL:
        print(
            "RESULT: PASS (worst max rel={0:.6e} <= {1:.0e}: GPU Horner polynomial "
            "agrees with numpy polyval at/near float64).".format(worst, PASS_REL)
        )
        return 0
    if worst < SUS_REL:
        print(
            "RESULT: SUSPICIOUS (worst max rel={0:.6e} in ({1:.0e},{2:.0e}). "
            "STOP; capture; do not fix.".format(worst, PASS_REL, SUS_REL)
        )
        return 2
    print(
        "RESULT: BUG (worst max rel={0:.6e} >= {1:.0e}. "
        "STOP; capture the numbers above; do not fix.".format(worst, SUS_REL)
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
