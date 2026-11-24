from __future__ import print_function

from more_itertools import prepend
from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import AbstractModel

import pyorbit.subroutines.kepler_exo as kepler_exo

from scipy.stats import norm

from types import SimpleNamespace

import pyorbit.common.gpu_backend as gpu_backend


def _active_kepler():
    """Dispatch the Keplerian / impact-parameter / apodization calls to the
    active backend (Round 1 GPU enablement).

    Returns a small namespace exposing the functions this module calls:

        rv_deltabjd, rv_semiamplitude, meanlong_from_dTc,
        deltaTc_from_meanlong, convert_rho_to_ars, convert_b_to_i, norm_pdf

    * ``gpu: true``  -> the JAX variants from
      ``pyorbit.subroutines.kepler_exo_jax`` and ``jax.scipy.stats.norm.pdf``.
      ``kepler_exo_jax`` is imported **lazily** here so the CPU path never loads
      a JAX-dependent module.
    * ``gpu: false`` (the default) -> the numpy ``kepler_exo`` functions and
      ``scipy.stats.norm.pdf`` exactly as before, so the CPU path is unchanged.

    The two sets are numerically identical under x64 (see
    ``development/gpu_cpu_parity.py``), so results do not depend on the backend.
    """
    if gpu_backend.gpu_enabled():
        import jax
        import pyorbit.subroutines.kepler_exo_jax as kepler_exo_jax

        return SimpleNamespace(
            rv_deltabjd=kepler_exo_jax.kepler_compute_rv_deltabjd_jax,
            rv_semiamplitude=kepler_exo_jax.kepler_compute_rv_semiamplitude_jax,
            meanlong_from_dTc=kepler_exo_jax.kepler_compute_meanlong_from_deltaTc_jax,
            deltaTc_from_meanlong=kepler_exo_jax.kepler_compute_deltaTc_from_meanlong_jax,
            convert_rho_to_ars=kepler_exo_jax.convert_rho_to_ars_jax,
            convert_b_to_i=kepler_exo_jax.convert_b_to_i_jax,
            norm_pdf=jax.scipy.stats.norm.pdf,
        )
    return SimpleNamespace(
        rv_deltabjd=kepler_exo.kepler_compute_rv_deltabjd,
        rv_semiamplitude=kepler_exo.kepler_compute_rv_semiamplitude,
        meanlong_from_dTc=kepler_exo.kepler_compute_meanlong_from_deltaTc,
        deltaTc_from_meanlong=kepler_exo.kepler_compute_deltaTc_from_meanlong,
        convert_rho_to_ars=convert_rho_to_ars,
        convert_b_to_i=convert_b_to_i,
        norm_pdf=norm.pdf,
    )


def _make_fused_rv_chain(variant):
    """F1 (Fusion): build -- or fetch from the ``gpu_backend.cached_jit`` graph
    cache -- a **single** ``jax.jit``-compiled graph for the Kepler-RV chain,
    so under ``gpu: true`` the whole chain

        meanlong_from_dTc -> (convert_rho_to_ars / convert_b_to_i)
        -> rv_semiamplitude -> rv_deltabjd   (+ optionally the apodizing factor)

    runs as one device-side program: intermediates never leave the GPU and no
    per-evaluation host<->device round trips remain between the stages.

    ``variant`` is a 5-tuple of *static* structural flags (cache key):

        (use_Tc, mass_family, compute_inclination, ars_from_rho, apodized)

    with ``mass_family`` in 'none' | 'mass' | 'scaled_mass' |
    'stellar_scaled_mass'. Every Python branch below is resolved over those
    static flags **at trace time** (``lax.cond``-safe: no data-dependent
    Python branching in the traced body). All per-evaluation values are
    *arguments*, so changing theta reuses the compiled graph without
    re-tracing.

    The one-pass Halley quirk of ``kepler_E_jax`` is preserved exactly (it is
    the documented parity anchor with the numpy solver) -- do NOT fix it here.
    The ``gpu: false`` path never enters this code.

    Argument layout (unpositioned values not read for a given variant are
    passed as harmless constants):

        x                 RV times (dataset.x0 or x0_input)
        x_aper            apodization times (dataset.x; = x when not apodized)
        P, c_par          period; Tc when use_Tc, else mean_long
        e, omega, Omega   eccentricity / argument of pericenter / node
        Tref              reference time (use_Tc branch only)
        K                 semi-amplitude (mass_family == 'none' only)
        mass, M_Me, Me_Ms, M_Ms          mass-family inputs
        b, i              impact parameter / inclination (incl inputs)
        density, a_Rs     scaled-a or a_Rs input (incl inputs)
        apo_center, apo_timescale        apodizing Gaussian (apodized only)
    """
    import jax

    import pyorbit.subroutines.kepler_exo_jax as kxj
    import pyorbit.subroutines.constants as kconst

    use_Tc, mass_family, incl, ars_from_rho, apodized = variant

    def chain(
        x,
        x_aper,
        P,
        c_par,
        e,
        omega,
        Omega,
        Tref,
        K,
        mass,
        M_Me,
        Me_Ms,
        M_Ms,
        b,
        i,
        density,
        a_Rs,
        apo_center,
        apo_timescale,
    ):
        if use_Tc:
            mean_long = kxj.kepler_compute_meanlong_from_deltaTc_jax(
                P, c_par - Tref, e, omega, Omega
            )
        else:
            mean_long = c_par

        if mass_family == "none":
            K_rv = K
        else:
            if mass_family == "mass":
                planet_mass = M_Me / kconst.Msear
            elif mass_family == "scaled_mass":
                planet_mass = Me_Ms * mass / kconst.Msear
            else:
                planet_mass = M_Ms * mass
            if incl:
                a_eff = kxj.convert_rho_to_ars_jax(P, density) if ars_from_rho else a_Rs
                i_eff = kxj.convert_b_to_i_jax(b, e, omega, a_eff)
            else:
                i_eff = i
            K_rv = kxj.kepler_compute_rv_semiamplitude_jax(
                mass, planet_mass, P, i_eff, e
            )

        rv = kxj.kepler_compute_rv_deltabjd_jax(x, K_rv, P, mean_long, e, omega, Omega)

        if apodized:
            pdf = jax.scipy.stats.norm.pdf
            rv = rv * (
                pdf(x_aper, apo_center, apo_timescale)
                / pdf(apo_center, apo_center, apo_timescale)
            )
        return rv

    key = "rv_fused|" + "|".join(str(v) for v in variant)
    return gpu_backend.cached_jit(chain, key)


class RVkeplerian(AbstractModel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = "radial_velocities"

        self.list_pams_common = OrderedSet(
            [
                "P",  # Period
                "K",  # RV semi-amplitude
                "e",  # eccentricity, uniform prior - to be fixed
                "omega",  # argument of pericenter
            ]
        )

        self.use_time_inferior_conjunction = False
        self.use_mass = False
        self.use_scaled_mass = False
        self.use_stellar_scaled_mass = False

        """ Default is to use the inclination value as provided.
        Otherwise, compute it using the semimajor axis derived from orbital period and stellar density"""
        self.compute_inclination = False
        self.compute_scaled_semimajor_axis = True
        self.compute_semimajor_axis = True

    # F1 (Fusion) instance caches. Both, once populated, hold the jit/vmap
    # callables assembled by :func:`_make_fused_rv_chain`; those objects are
    # local closures that ``pickle`` cannot capture, while the samplers'
    # checkpoints (``pyde_save_to_pickle`` / emcee reload) pickle the whole
    # model container. Neither cache holds any per-evaluation values -- both
    # are rebuilt lazily and deterministically from the structural
    # ``variant`` key alone on the next ``compute`` / ``compute_batch`` call,
    # so dropping them from a pickle round-trip is lossless. On the
    # ``gpu: false`` path neither attribute is ever created, so this pair is
    # a graceful no-op there.
    _fused_cache_attrs = ("_fused_chain_cache", "_fused_batch_cache")

    def __getstate__(self):
        """Drop the (unpicklable) fusion caches from the instance state; the
        lazy cache-miss path rebuilds them on first use after unpickle."""
        state = vars(self)
        return {
            key: value
            for key, value in state.items()
            if key not in self._fused_cache_attrs
        }

    def __setstate__(self, state):
        """Restore the instance state; leave any fusion cache absent so the
        existing cache-miss ``getattr(..., None)`` path rebuilds it."""
        vars(self).update(state)
        for attr in self._fused_cache_attrs:
            if attr in vars(self):
                del vars(self)[attr]

    def initialize_model(self, mc, **kwargs):

        # Copying the property to the class for faster access
        self.use_time_inferior_conjunction = getattr(
            mc.common_models[self.planet_ref],
            "use_time_inferior_conjunction",
            self.use_time_inferior_conjunction,
        )
        self.use_mass = getattr(
            mc.common_models[self.planet_ref], "use_mass", self.use_mass
        )
        self.use_scaled_mass = getattr(
            mc.common_models[self.planet_ref], "use_scaled_mass", self.use_scaled_mass
        )
        self.use_stellar_scaled_mass = getattr(
            mc.common_models[self.planet_ref],
            "use_stellar_scaled_mass",
            self.use_stellar_scaled_mass,
        )
        self.compute_inclination = getattr(
            mc.common_models[self.planet_ref],
            "compute_inclination",
            self.compute_inclination,
        )
        self.compute_scaled_semimajor_axis = getattr(
            mc.common_models[self.planet_ref],
            "compute_scaled_semimajor_axis",
            self.compute_scaled_semimajor_axis,
        )

        if (
            mc.common_models[self.planet_ref].parametrization[:8] == "Ford2006"
            and mc.common_models[self.planet_ref].orbit != "circular"
        ):
            self.list_pams_common.discard("e")
            self.list_pams_common.discard("omega")

            self.list_pams_common.update(["e_coso"])
            self.list_pams_common.update(["e_sino"])

        elif (
            mc.common_models[self.planet_ref].parametrization[:8] != "Standard"
            and mc.common_models[self.planet_ref].orbit != "circular"
        ):
            # 'Eastman2013' is the standard choice
            self.list_pams_common.discard("e")
            self.list_pams_common.discard("omega")

            self.list_pams_common.update(["sre_coso"])
            self.list_pams_common.update(["sre_sino"])

        if self.use_time_inferior_conjunction:
            self.list_pams_common.update(["Tc"])
        else:
            self.list_pams_common.update(["mean_long"])

        if self.use_mass or self.use_scaled_mass or self.use_stellar_scaled_mass:
            self.list_pams_common.discard("K")

            if self.use_mass:
                self.list_pams_common.update(["M_Me"])
            elif self.use_scaled_mass:
                self.list_pams_common.update(["Me_Ms"])
            elif self.use_stellar_scaled_mass:
                self.list_pams_common.update(["M_Ms"])
            self.list_pams_common.update(["mass"])

            if self.compute_inclination:
                """ b is the impact parameter """
                self.list_pams_common.update(["b"])
                if self.compute_scaled_semimajor_axis:
                    self.list_pams_common.update(["density"])
                else:
                    self.list_pams_common.update(["a_Rs"])
            else:
                """ i is the orbital inclination (in degrees) """
                self.list_pams_common.update(["i"])
        else:
            self.list_pams_common.update(["K"])

        if mc.common_models[self.planet_ref].use_longitude_of_nodes:
            self.list_pams_common.update(["Omega"])

        try:
            self.default_Omega = mc.common_models[self.planet_ref].default_Omega
        except AttributeError:
            self.default_Omega = 0.00

    def compute(self, parameter_values, dataset, x0_input=None):

        if gpu_backend.gpu_enabled():
            # F1 (Fusion): single-graph GPU chain; CPU path below unchanged.
            return self._fused_compute(parameter_values, dataset, x0_input)

        kx = _active_kepler()

        Omega = parameter_values.get("Omega", self.default_Omega)

        if self.use_time_inferior_conjunction:
            mean_long = kx.meanlong_from_dTc(
                parameter_values["P"],
                parameter_values["Tc"] - dataset.Tref,
                parameter_values["e"],
                parameter_values["omega"],
                Omega,
            )
        else:
            mean_long = parameter_values["mean_long"]

        if self.use_mass or self.use_scaled_mass or self.use_stellar_scaled_mass:
            if self.compute_inclination:
                if self.compute_scaled_semimajor_axis:
                    parameter_values["a_Rs"] = kx.convert_rho_to_ars(
                        parameter_values["P"], parameter_values["density"]
                    )

                parameter_values["i"] = kx.convert_b_to_i(
                    parameter_values["b"],
                    parameter_values["e"],
                    parameter_values["omega"],
                    parameter_values["a_Rs"],
                )

            if self.use_mass:
                planet_mass = parameter_values["M_Me"] / constants.Msear
            elif self.use_scaled_mass:
                planet_mass = (
                    parameter_values["Me_Ms"]
                    * parameter_values["mass"]
                    / constants.Msear
                )
            elif self.use_stellar_scaled_mass:
                planet_mass = parameter_values["M_Ms"] * parameter_values["mass"]

            rv_semiamplitude = kx.rv_semiamplitude(
                parameter_values["mass"],
                planet_mass,
                parameter_values["P"],
                parameter_values["i"],
                parameter_values["e"],
            )

        else:
            rv_semiamplitude = parameter_values["K"]

        if x0_input is None:
            return kx.rv_deltabjd(
                dataset.x0,
                rv_semiamplitude,
                parameter_values["P"],
                mean_long,
                parameter_values["e"],
                parameter_values["omega"],
                Omega,
            )
        else:
            return kx.rv_deltabjd(
                x0_input,
                rv_semiamplitude,
                parameter_values["P"],
                mean_long,
                parameter_values["e"],
                parameter_values["omega"],
                Omega,
            )

    def _fused_compute(self, parameter_values, dataset, x0_input=None, apodized=False):
        """F1 (Fusion): ``gpu: true`` evaluation of the Kepler-RV chain through
        the single jit-compiled graph from :func:`_make_fused_rv_chain` (see its
        docstring for the variant/parity rules). All per-evaluation values go
        in as device arguments; nothing here syncs back to the host -- the
        single host materialisation for the evaluation happens when the
        consumer extracts the final value. ``apodized=True`` replicates the
        ApodizedRVkeplerian eager semantics (inclination never computed, i
        pinned to 90 deg; the apodizing factor evaluated at ``dataset.x``).

        Called **only** in the ``gpu_enabled()`` branches of ``compute``; the
        ``gpu: false`` path is byte-for-byte unchanged."""
        import jax.numpy as jnp

        pv = parameter_values
        use_Tc = bool(self.use_time_inferior_conjunction)
        if apodized:
            # the eager apodized mass-branch pins i = 90.0 and never reads
            # b / a_Rs / density -- mirror that exactly.
            incl = False
            ars_from_rho = False
            i = 90.0
        else:
            incl = bool(self.compute_inclination)
            ars_from_rho = bool(self.compute_scaled_semimajor_axis)
            i = pv.get("i", 0.0)
        if self.use_mass:
            mass_family = "mass"
        elif self.use_scaled_mass:
            mass_family = "scaled_mass"
        elif self.use_stellar_scaled_mass:
            mass_family = "stellar_scaled_mass"
        else:
            mass_family = "none"
        variant = (use_Tc, mass_family, incl, bool(ars_from_rho), bool(apodized))

        chains = getattr(self, "_fused_chain_cache", None)
        if chains is None:
            chains = {}
            self._fused_chain_cache = chains
        chain = chains.get(variant)
        if chain is None:
            chain = _make_fused_rv_chain(variant)
            chains[variant] = chain

        if x0_input is None:
            x0_dev = getattr(dataset, "x0_dev", None)
            x = x0_dev if x0_dev is not None else jnp.asarray(dataset.x0)
            if apodized:
                x_dev = getattr(dataset, "x_dev", None)
                xa = x_dev if x_dev is not None else jnp.asarray(dataset.x)
            else:
                xa = x
        else:
            x = jnp.asarray(x0_input)
            xa = jnp.asarray(dataset.x) if apodized else x

        return chain(
            x,
            xa,
            pv["P"],
            pv["Tc"] if use_Tc else pv["mean_long"],
            pv["e"],
            pv["omega"],
            pv.get("Omega", self.default_Omega),
            dataset.Tref,
            pv.get("K", 0.0),
            pv.get("mass", 0.0),
            pv.get("M_Me", 0.0),
            pv.get("Me_Ms", 0.0),
            pv.get("M_Ms", 0.0),
            pv.get("b", 0.0),
            i,
            pv.get("density", 0.0),
            pv.get("a_Rs", 0.0),
            pv.get("apo_center", 0.0),
            pv.get("apo_timescale", 0.0),
        )

    def compute_batch(self, parameter_values_list, dataset):
        """F1-batch: M Kepler-RV curves in a **single** jit-compiled GPU call.

        ``parameter_values_list`` is a list of M ``parameter_values`` dicts
        (same keys/semantics as :meth:`compute`). The per-theta values are
        stacked into (M,) argument columns and the single-theta fused chain
        from :func:`_make_fused_rv_chain` is ``jax.vmap``-ed over them, so the
        whole (M, N) RV surface is produced by ONE jitted program: the
        fixed per-call host<->device dispatch cost is paid once, not M times
        (proved in /tmp/opencode/pyorbit_batch_proof/batch_proof.py: 73-205x
        for M=24, maxdiff 5.8e-11 vs M separate ``compute`` calls).

        * ``gpu: true``  -> returns a (M, N) JAX device array; nothing syncs
          to the host here -- the single materialisation happens when the
          consumer reads the result (``.block_until_ready()`` /
          ``np.asarray``), not per theta.
        * ``gpu: false`` -> the batching gain only exists on-device, so this
          delegates to the existing single-theta loop and returns a list of
          M numpy arrays (CPU path byte-for-byte unchanged).
        """
        pvs = list(parameter_values_list)
        if not pvs:
            raise ValueError("compute_batch: parameter_values_list must not be empty")

        if not gpu_backend.gpu_enabled():
            return [self.compute(pv, dataset) for pv in pvs]

        import jax
        import jax.numpy as jnp

        use_Tc = bool(self.use_time_inferior_conjunction)
        """ ApodizedRVkeplerian advertises the apodizing parameters: the GPU
        path must then mirror ``_fused_compute(apodized=True)`` (factor
        evaluated at dataset.x, i pinned to 90 deg, inclination never
        computed). """
        apodized = "apo_center" in getattr(self, "list_pams_common", ())
        if apodized:
            incl = False
            ars_from_rho = False
            i_vals = [90.0] * len(pvs)
        else:
            incl = bool(self.compute_inclination)
            ars_from_rho = bool(self.compute_scaled_semimajor_axis)
            i_vals = [pv.get("i", 0.0) for pv in pvs]
        if self.use_mass:
            mass_family = "mass"
        elif self.use_scaled_mass:
            mass_family = "scaled_mass"
        elif self.use_stellar_scaled_mass:
            mass_family = "stellar_scaled_mass"
        else:
            mass_family = "none"
        variant = (use_Tc, mass_family, incl, bool(ars_from_rho), bool(apodized))

        chains = getattr(self, "_fused_chain_cache", None)
        if chains is None:
            chains = {}
            self._fused_chain_cache = chains
        chain = chains.get(variant)
        if chain is None:
            chain = _make_fused_rv_chain(variant)
            chains[variant] = chain

        # (M,) per-theta argument columns -- same keys/defaults as
        # ``_fused_compute`` (``pv["..."]`` where it indexes, ``pv.get``
        # where it falls back) so a theta missing a key fails identically.
        # Assembled as one (16, M) block -> a SINGLE host->device transfer
        # per call (16 separate uploads would each pay a dispatch round trip).
        cols = [
            [pv["P"] for pv in pvs],
            [pv["Tc"] if use_Tc else pv["mean_long"] for pv in pvs],
            [pv["e"] for pv in pvs],
            [pv["omega"] for pv in pvs],
            [pv.get("Omega", self.default_Omega) for pv in pvs],
            [pv.get("K", 0.0) for pv in pvs],
            [pv.get("mass", 0.0) for pv in pvs],
            [pv.get("M_Me", 0.0) for pv in pvs],
            [pv.get("Me_Ms", 0.0) for pv in pvs],
            [pv.get("M_Ms", 0.0) for pv in pvs],
            [pv.get("b", 0.0) for pv in pvs],
            i_vals,
            [pv.get("density", 0.0) for pv in pvs],
            [pv.get("a_Rs", 0.0) for pv in pvs],
            [pv.get("apo_center", 0.0) for pv in pvs],
            [pv.get("apo_timescale", 0.0) for pv in pvs],
        ]
        x0_dev = getattr(dataset, "x0_dev", None)
        x = x0_dev if x0_dev is not None else jnp.asarray(dataset.x0)
        if apodized:
            x_dev = getattr(dataset, "x_dev", None)
            xa = x_dev if x_dev is not None else jnp.asarray(dataset.x)
        else:
            xa = x

        batches = getattr(self, "_fused_batch_cache", None)
        if batches is None:
            batches = {}
            self._fused_batch_cache = batches
        batched = batches.get(variant)
        if batched is None:
            """ The 16 per-theta columns enter as ONE (M, 16) block (a single
            host->device transfer and a single jit argument); the columns are
            split *inside* the traced body, where slicing is free, and
            ``chain`` is vmapped over theta. The shared per-dataset inputs
            (x, apodizing times, Tref) enter as plain unbatched arguments, so
            the graph is correct for any dataset. """

            def _single(
                x_a,
                xa_a,
                Tref,
                P,
                cpar,
                e,
                omega,
                Omega,
                K,
                mass,
                MMe,
                MeMs,
                MMs,
                b_,
                i_,
                density,
                aRs,
                apo_c,
                apo_t,
            ):
                return chain(
                    x_a,
                    xa_a,
                    P,
                    cpar,
                    e,
                    omega,
                    Omega,
                    Tref,
                    K,
                    mass,
                    MMe,
                    MeMs,
                    MMs,
                    b_,
                    i_,
                    density,
                    aRs,
                    apo_c,
                    apo_t,
                )

            def _batched_body(x_a, xa_a, Tref, cols):
                (
                    P,
                    cpar,
                    e,
                    omega,
                    Omega,
                    K,
                    mass,
                    MMe,
                    MeMs,
                    MMs,
                    b_,
                    i_,
                    density,
                    aRs,
                    apo_c,
                    apo_t,
                ) = (cols[:, k] for k in range(16))
                return jax.vmap(_single, in_axes=(None, None, None) + (0,) * 16)(
                    x_a,
                    xa_a,
                    Tref,
                    P,
                    cpar,
                    e,
                    omega,
                    Omega,
                    K,
                    mass,
                    MMe,
                    MeMs,
                    MMs,
                    b_,
                    i_,
                    density,
                    aRs,
                    apo_c,
                    apo_t,
                )

            batched = jax.jit(_batched_body)
            batches[variant] = batched

        return batched(x, xa, dataset.Tref, jnp.asarray(cols).T)


class ApodizedRVkeplerian(RVkeplerian):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.list_pams_common.update(["apo_center"])
        self.list_pams_common.update(["apo_timescale"])

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        try:
            min_apo_center = min(
                mc.common_models[self.planet_ref].default_bounds["apo_center"][0],
                np.min(dataset.x),
            )
            max_apo_center = max(
                mc.common_models[self.planet_ref].default_bounds["apo_center"][1],
                np.max(dataset.x),
            )
        except KeyError:
            min_apo_center = np.min(dataset.x)
            max_apo_center = np.max(dataset.x)

        mc.common_models[self.planet_ref].default_bounds["apo_center"] = [
            min_apo_center,
            max_apo_center,
        ]
        return

    def compute(self, parameter_values, dataset, x0_input=None):

        if gpu_backend.gpu_enabled():
            # F1 (Fusion): single-graph GPU chain (apodizing factor included);
            # CPU path below unchanged.
            return self._fused_compute(
                parameter_values, dataset, x0_input, apodized=True
            )

        kx = _active_kepler()

        Omega = parameter_values.get("Omega", self.default_Omega)

        if self.use_time_inferior_conjunction:
            mean_long = kx.meanlong_from_dTc(
                parameter_values["P"],
                parameter_values["Tc"] - dataset.Tref,
                parameter_values["e"],
                parameter_values["omega"],
                Omega,
            )
        else:
            mean_long = parameter_values["mean_long"]

        if self.use_mass or self.use_scaled_mass or self.use_stellar_scaled_mass:
            if self.use_mass:
                planet_mass = parameter_values["M_Me"] / constants.Msear
            elif self.use_scaled_mass:
                planet_mass = (
                    parameter_values["Me_Ms"]
                    * parameter_values["mass"]
                    / constants.Msear
                )
            elif self.use_stellar_scaled_mass:
                planet_mass = parameter_values["M_Ms"] * parameter_values["mass"]

            rv_semiamplitude = kx.rv_semiamplitude(
                parameter_values["mass"],
                planet_mass,
                parameter_values["P"],
                90.0,
                parameter_values["e"],
            )
        else:
            rv_semiamplitude = parameter_values["K"]

        if x0_input is None:
            # Apodization factor
            apo_function = kx.norm_pdf(
                dataset.x,
                parameter_values["apo_center"],
                parameter_values["apo_timescale"],
            ) / kx.norm_pdf(
                parameter_values["apo_center"],
                parameter_values["apo_center"],
                parameter_values["apo_timescale"],
            )

            return apo_function * kx.rv_deltabjd(
                dataset.x0,
                rv_semiamplitude,
                parameter_values["P"],
                mean_long,
                parameter_values["e"],
                parameter_values["omega"],
                Omega,
            )
        else:
            # Apodization factor
            apo_function = kx.norm_pdf(
                dataset.x,
                parameter_values["apo_center"],
                parameter_values["apo_timescale"],
            ) / kx.norm_pdf(
                parameter_values["apo_center"],
                parameter_values["apo_center"],
                parameter_values["apo_timescale"],
            )

            return apo_function * kx.rv_deltabjd(
                dataset.x0,
                rv_semiamplitude,
                parameter_values["P"],
                mean_long,
                parameter_values["e"],
                parameter_values["omega"],
                Omega,
            )


class TransitTimeKeplerian(AbstractModel):
    model_class = "transit_times"

    def __init__(self, *args, **kwargs):
        super(TransitTimeKeplerian, self).__init__(*args, **kwargs)
        self.use_time_inferior_conjunction = False

        self.list_pams_common = OrderedSet(["P"])  # Period

        self.list_pams_dataset = OrderedSet()

    def initialize_model(self, mc, **kwargs):

        if mc.common_models[self.planet_ref].use_time_inferior_conjunction:
            self.list_pams_common.update(["Tc"])
            self.use_time_inferior_conjunction = True
            # Copying the property to the class for faster access
        else:
            self.list_pams_common.update(["mean_long"])
            self.list_pams_common.update(["e"])
            self.list_pams_common.update(["omega"])
            # mean longitude = argument of pericenter + mean anomaly at Tref

        try:
            self.default_Omega = mc.common_models[self.planet_ref].default_Omega
        except AttributeError:
            self.default_Omega = 0.00

    def compute(self, parameter_values, dataset, x0_input=None):

        kx = _active_kepler()

        # TODO check this function, I believe it's a bit crazy.... URGENT
        Omega = parameter_values.get("Omega", self.default_Omega)

        if self.use_time_inferior_conjunction:
            delta_T = (
                parameter_values["Tc"]
                - np.floor(
                    (parameter_values["Tc"] - dataset.Tref) / parameter_values["P"]
                )
                * parameter_values["P"]
            )
        else:
            delta_T = dataset.Tref + kx.deltaTc_from_meanlong(
                parameter_values["P"],
                parameter_values["mean_long"],
                parameter_values["e"],
                parameter_values["omega"],
                Omega,
            )

        if x0_input is None:
            return (
                np.floor(dataset.x0 / parameter_values["P"]) * parameter_values["P"]
                + delta_T
            )
        else:
            return (
                np.floor(x0_input / parameter_values["P"]) * parameter_values["P"]
                + delta_T
            )
