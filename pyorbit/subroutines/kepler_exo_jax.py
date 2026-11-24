"""JAX (GPU) port of the RV-Keplerian solver from ``pyorbit.subroutines.kepler_exo``.

This module is the *only* place in PyORBIT that re-implements the Keplerian
radial-velocity solver on JAX, and it is imported **lazily** -- only when the
``gpu: true`` flag is enabled (see ``pyorbit.common.gpu_backend``). The CPU path
(``gpu: false``, the default) never touches this file, so its behaviour is
byte-for-byte unchanged.

Every function here is a direct port of the corresponding numpy function in
``kepler_exo.py`` (and ``common.py``), with the two JAX-specific changes the
plan calls for:

* ``kepler_E`` -- the shipped numpy solver performs **at most one**
  Fulton/Halley 4th-order Newton step: its convergence bookkeeping
  ``nd = np.sum(convd is True)`` is a Python *identity* check on an ndarray,
  which is always False, so ``nd`` is always 0 and the ``while nd > 0`` loop
  exits after the single pass. To keep the ``gpu: true`` result **bit-for-bit
  identical** to ``gpu: false`` (Round 1's hard parity constraint) we replicate
  exactly that one-pass behaviour -- a single Halley step on the
  not-yet-converged elements, then stop -- rather than iterating to the
  ``1e-12`` tolerance. (Converging further would *disagree* with the CPU;
  see the note on ``kepler_E_jax`` and the latent-bug call-out there.)
* ``kepler_compute_rv_deltabjd`` -- the ``if abs(ecc) < 1e-3`` / ``if ecc < 0``
  branches become ``jax.lax.cond`` (the traced eccentricity cannot drive a
  Python ``if``).

All other functions are pure arithmetic and are a straight ``np`` -> ``jnp``
port. Under x64 (``gpu: true`` forces ``jax_enable_x64``) the results are
numerically identical to the numpy originals.
"""

import jax
import jax.numpy as jnp

import pyorbit.subroutines.constants as constants


def kepler_E_jax(M_in, ecc):
    """JAX port of ``kepler_exo.kepler_E`` (bit-for-bit, one-pass).

    Solves Kepler's equation ``E - e sin(E) = M`` for the eccentric anomaly
    ``E`` given the mean anomaly ``M`` and eccentricity ``e``.

    **Why one pass and not a convergence loop:** the shipped numpy
    ``kepler_E`` (``kepler_exo.py:84``) is *intended* to iterate to a ``1e-12``
    tolerance, but its loop-exit line ``nd = np.sum(convd is True)`` is a Python
    *identity* check (``ndarray is True`` is always False), so ``nd`` is always
    0 and the loop runs exactly **once**. The shipped CPU result therefore
    carries a residual of order ``1e-6`` at the highest eccentricities.

    Round 1's hard requirement is that ``gpu: true`` reproduces ``gpu: false``
    bit-for-bit, so we deliberately replicate that single-pass behaviour
    (one Halley step on the not-yet-converged elements, frozen where the
    initial guess already meets ``conv``) instead of converging further --
    converging would be *more* accurate but would **disagree** with the CPU.

    .. note:: The underlying one-iteration quirk is a latent bug in
       ``kepler_exo.py``. If you want the full ``1e-12`` accuracy it *must* be
       fixed in the **numpy** path and this JAX port **together**, in the same
       round, or ``gpu:false`` and ``gpu:true`` will diverge.
    """
    mx = jnp.atleast_1d(jnp.asarray(M_in)) % (2. * jnp.pi)

    # Using Fulton scheme
    k = 0.85
    e0 = mx + jnp.sign(jnp.sin(mx)) * k * ecc  # first guess at E

    conv = 1.0e-12  # convergence criterion
    # Elements whose initial guess already meets `conv` are left untouched,
    # exactly like the numpy's `convd` mask.
    not_conv = jnp.abs(e0 - ecc * jnp.sin(e0) - mx) > conv

    se0 = jnp.sin(e0)
    ce0 = jnp.cos(e0)
    f = e0 - ecc * se0 - mx
    f1 = 1. - ecc * ce0
    f2 = ecc * se0
    f3 = ecc * ce0
    u1 = -f / f1
    u2 = -f / (f1 + 0.5 * f2 * u1)
    u3 = -f / (f1 + 0.5 * f2 * u2 + (1. / 6.) * f3 * u2 ** 2.)
    u4 = -f / (f1 + 0.5 * f2 * u3 + (1. / 6.) * f3 * u3 ** 2 - (1. / 24.) * f2 * u3 ** 3)

    # One Halley step; freeze the already-converged elements (numpy one-pass).
    e0 = jnp.where(not_conv, e0 + u4, e0)
    ecan_tmp = e0 % (2. * jnp.pi)
    return ecan_tmp


def kepler_compute_rv_semiamplitude_jax(mass_primary, mass_secondary, period, inclination, eccentricity):
    """JAX port of ``kepler_exo.kepler_compute_rv_semiamplitude`` (pure arithmetic)."""
    return (2. * jnp.pi * constants.Gsi * constants.Msun / 86400.) ** (1. / 3.) \
           * (jnp.sin(inclination * constants.deg2rad) / jnp.sqrt(1. - eccentricity ** 2.)) * period ** (-1. / 3.) \
           * (mass_secondary * (mass_primary + mass_secondary) ** (-2. / 3.))


def kepler_compute_rv_deltabjd_jax(bjd_tref, rv_semiamplitude, period, mean_long, ecc, omega_deg, Omega_deg=0.0):
    """JAX port of ``kepler_exo.kepler_compute_rv_deltabjd``.

    The two eccentricity branches are expressed with ``jax.lax.cond`` so the
    function traces under ``jax.jit``.
    """
    omega = jnp.asarray(omega_deg) * constants.deg2rad
    Omega = jnp.asarray(Omega_deg) * constants.deg2rad
    MeAn = (2. * jnp.pi * (bjd_tref / period) + (mean_long * constants.deg2rad - omega - Omega)) % (2. * jnp.pi)
    # Normalise MeAn to >=1-D up front. The eccentric branch routes it through
    # kepler_E_jax -- which, like the numpy kepler_E, applies atleast_1d -- while
    # the circular branch returns it directly. jax.lax.cond type-checks *both*
    # branches at trace time and requires identical output shapes, so the raw
    # `()` scalar the circular branch would otherwise return mismatches the
    # eccentric branch's `(1,)`. Doing it here keeps both branches consistent for
    # scalar (unit-test) and array (real sampling) inputs alike, changing no value.
    MeAn = jnp.atleast_1d(MeAn)

    def circular(_):
        return (MeAn, jnp.asarray(0.), omega)

    def eccentric(_):
        def pos(_):
            return (ecc, omega)
        def neg(_):
            return (-ecc, omega + jnp.pi)

        ecc_eq, omega_e = jax.lax.cond(ecc < 0., neg, pos, None)
        # Eccentric Anomaly
        EccAn = kepler_E_jax(MeAn, ecc_eq)
        TrAn = 2. * jnp.arctan(jnp.sqrt((1.0 + ecc_eq) / (1.0 - ecc_eq)) * jnp.tan(EccAn / 2.0))
        return (TrAn, ecc_eq, omega_e)

    TrAn, ecc_eq, omega = jax.lax.cond(jnp.abs(ecc) < 1e-3, circular, eccentric, None)
    rv = rv_semiamplitude * (jnp.cos(TrAn + omega) + ecc_eq * jnp.cos(omega))
    return rv


def kepler_compute_deltaTc_from_meanlong_jax(period, mean_long, ecc=0., omega_deg=90., Omega_deg=0.0):
    """JAX port of ``kepler_exo.kepler_compute_deltaTc_from_meanlong`` (pure arithmetic)."""
    omega = omega_deg * constants.deg2rad
    Omega = Omega_deg * constants.deg2rad
    TrAn = jnp.pi / 2 - omega
    EccAn = 2. * jnp.arctan(jnp.sqrt((1.0 - ecc) / (1.0 + ecc)) * jnp.tan(TrAn / 2.0))
    MeAn = EccAn - ecc * jnp.sin(EccAn)
    MeAn_Tref = mean_long * constants.deg2rad - omega - Omega
    return (MeAn - MeAn_Tref) / (2. * jnp.pi) * period % period


def kepler_compute_meanlong_from_deltaTc_jax(period, delta_Tc, ecc, omega_deg, Omega_deg=0.0):
    """JAX port of ``kepler_exo.kepler_compute_meanlong_from_deltaTc`` (pure arithmetic)."""
    omega = omega_deg * constants.deg2rad
    Omega = Omega_deg * constants.deg2rad
    TrAn = jnp.pi / 2 - omega
    EccAn = 2. * jnp.arctan(jnp.sqrt((1.0 - ecc) / (1.0 + ecc)) * jnp.tan(TrAn / 2.0))
    MeAn = EccAn - ecc * jnp.sin(EccAn)
    return ((omega + Omega + MeAn - delta_Tc / period * 2. * jnp.pi) * constants.rad2deg) % 360.


def convert_rho_to_ars_jax(P, rho):
    """JAX port of ``common.convert_rho_to_ars`` (pure arithmetic)."""
    return jnp.power(constants.Gsi * (constants.d2s * constants.d2s) * (P ** 2)
                     * rho * constants.rho_Sun / (3. * jnp.pi), 1. / 3.)


def convert_b_to_i_jax(b, e, o, a):
    """JAX port of ``common.convert_b_to_i``.

    The Python list-comprehension clamp on ``arccos_argument`` becomes a single
    ``jnp.clip`` (works for both scalar and array inputs).
    """
    o_rad = o * constants.deg2rad
    rho_e = (1. - e ** 2) / (1. + e * jnp.sin(o_rad))
    arccos_argument = b / a / rho_e
    arccos_argument = jnp.clip(arccos_argument, -1., 1.)
    return jnp.arccos(arccos_argument) * 180. / jnp.pi
