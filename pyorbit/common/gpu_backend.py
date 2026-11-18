"""
GPU backend foundation for PyORBIT.

This module is the single place that decides whether PyORBIT runs its JAX-native
models (tinygp GPs, pytransit transits, and the JAX port of the RV Keplerian
solver) on a GPU, and it enforces the "gate + guide" contract:

  * ``gpu: false`` (the default) -> this module is a complete no-op. PyORBIT
    behaves exactly as before: pure CPU, no JAX device change, no model rejected,
    no new import side effects.
  * ``gpu: true`` -> JAX is forced to float64 and its *default* backend is
    retargeted to the GPU **before any sampling happens**. If a GPU is not
    actually visible to the installed ``jaxlib``, an actionable
    :class:`GPUSupportedModelError` is raised (never a silent CPU fallback).
    Models that are not JAX-native are rejected up-front by :func:`gate`, with a
    message naming the JAX equivalent to use.

The module deliberately imports ``jax`` only inside functions so that importing
``pyorbit`` in a CPU-only environment is unaffected.
"""

from __future__ import print_function

import os

__all__ = [
    "GPUSupportedModelError",
    "initialize",
    "gpu_enabled",
    "get_device",
    "to_device",
    "jitted",
    "cached_jit",
    "xla_allocator_envs",
    "gate",
]


class GPUSupportedModelError(Exception):
    """Raised when a model is requested on the GPU that has no JAX backend, or
    when ``gpu: true`` is requested but JAX cannot actually reach a GPU."""

    pass


# --- module-level runtime state (all CPU by default) -----------------------
_enabled = False
_device = None
_verbose = False
_initialized = False
# F1 (Fusion): per-structural-signature cache of compiled single-graph
# functions, so fused chains are traced exactly once, not per evaluation.
_GRAPH_CACHE = {}


def gpu_enabled():
    """True once :func:`initialize` has retargeted JAX to the GPU in this process."""
    return _enabled


def _import_jax():
    # Imported here (not at module import time) so that a CPU-only environment
    # can import this module -- and therefore `pyorbit` -- without JAX installed.
    import jax
    import jax.numpy as jnp

    return jax, jnp


def initialize(parameters=None):
    """Apply the GPU configuration described by ``parameters``.

    ``parameters`` is the ``parameters:`` section of the PyORBIT YAML config
    (a plain dict). Recognised keys (all optional, all default to CPU):

      * ``gpu``          -- master switch (default ``False``)
      * ``gpu_device``   -- CUDA device index (default ``0``)
      * ``gpu_verbose``  -- log ``jax.devices()`` etc. (default ``False``)

    With ``gpu: false`` this function is a no-op (JAX is never touched). With
    ``gpu: true`` it forces float64, retargets the JAX default backend to the
    GPU, and raises :class:`GPUSupportedModelError` if a GPU is not reachable.
    """
    global _enabled, _device, _verbose, _initialized

    if parameters is None:
        parameters = {}

    gpu_flag = bool(parameters.get("gpu", False))
    gpu_verbose = bool(parameters.get("gpu_verbose", False))
    gpu_device = parameters.get("gpu_device", 0)

    if not gpu_flag:
        # CPU default: do NOT touch JAX config in any way.
        _enabled = False
        _device = None
        _verbose = False
        _initialized = True
        return _enabled

    _verbose = gpu_verbose
    # F2 (Fusion): XLA allocator defaults for GPU legs, applied before jaxlib
    # allocates any GPU memory. Preallocating a deterministic arena (~0.30 of
    # the card) fixes the tinygp GPU OOM observed at N=10k under
    # PREALLOCATE=false + MEM_FRACTION=0.1 while the card had 31 GiB free;
    # when less than ~10 GiB is free (shared card, e.g. ollama holding VRAM)
    # we fall back to on-demand 0.10. PYORBIT_XLA_PREALLOCATE /
    # PYORBIT_XLA_MEM_FRACTION (or the XLA_* names) override both.
    for _k, _v in xla_allocator_envs().items():
        os.environ.setdefault(_k, _v)
    jax, jnp = _import_jax()

    # Best-effort device selection. CUDA_VISIBLE_DEVICES is only meaningful
    # before JAX's CUDA context is created, which (in PyORBIT) happens on the
    # first GPU op during sampling, i.e. strictly after this call.
    if gpu_device is not None and "CUDA_VISIBLE_DEVICES" not in os.environ:
        os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_device)

    # Can this jaxlib actually see a GPU?
    devices_error = None
    try:
        gpu_devices = list(jax.devices("gpu"))
    except Exception as exc:
        gpu_devices = []
        # Keep the actionable error below, but do not swallow the original
        # failure: embed its text so the real cause (e.g. a bad pre-import
        # JAX_PLATFORMS key) is visible instead of a misleading "no GPU".
        devices_error = exc

    if not gpu_devices:
        err_detail = ""
        if devices_error is not None:
            err_detail = "    JAX reported: {0}: {1}\n".format(
                type(devices_error).__name__, devices_error
            )
        raise GPUSupportedModelError(
            "gpu: true was requested, but JAX cannot see any GPU.\n"
            + err_detail
            + "Most likely cause: the CPU-only 'jaxlib' wheel is installed.\n"
            "Fix -- install the PyORBIT [gpu] extra (modern JAX 0.6.2 CUDA stack), plus a CUDA driver:\n"
            "    pip install pyorbit-package[gpu]\n"
            "If your shell exports JAX_PLATFORMS=gpu, unset it first -- PyORBIT's pre-import hooks can't override an already-set bad export (on jaxlib 0.6.x it resolves to the absent 'rocm' backend).\n"
            "If a GPU is present but JAX still reports none, launch with:\n"
            "    export JAX_PLATFORMS=cuda\n"
            "    python PyORBIT_Run.py <sampler> <config>"
        )

    # Always float64 on the GPU (user decision: no float32 path).
    try:
        jax.config.update("jax_enable_x64", True)
    except Exception:
        os.environ.setdefault("JAX_ENABLE_X64", "true")

    # Retarget JAX's *default* backend to the GPU for subsequently-created
    # arrays. This must happen before any sampling, and it does: initialize() is
    # called immediately after the config is parsed and before the sampler ever
    # evaluates mc(theta). tinygp/pytransit kernels, jitted on first call during
    # sampling, inherit this GPU default.
    try:
        jax.config.update("jax_platform_name", "gpu")
    except Exception:
        # Fallback env paths run pre-import, where the concrete 'cuda' key is
        # required (on jaxlib 0.6.x 'JAX_PLATFORMS=gpu' resolves to the absent
        # 'rocm' implementation and kills backend init).
        os.environ.setdefault("JAX_PLATFORM_NAME", "cuda")
        os.environ.setdefault("JAX_PLATFORMS", "cuda")

    try:
        default_backend = str(jax.default_backend()).lower()
    except Exception as exc:  # pragma: no cover - defensive
        raise GPUSupportedModelError(
            "gpu: true was requested but JAX could not resolve a default backend "
            "({0}).".format(exc)
        )

    if default_backend != "gpu":
        raise GPUSupportedModelError(
            "gpu: true was requested and JAX sees a GPU, but the default backend "
            "resolved to '{0}'. The GPU would NOT be used; refusing to silently "
            "fall back to the CPU.\n"
            "Fix: install the CUDA jaxlib wheel, or set JAX_PLATFORMS=cuda and "
            "CUDA_VISIBLE_DEVICES before launching, or run:\n"
            "    python PyORBIT_Run.py <sampler> <config>\n"
            "Visible devices: {1}".format(default_backend, jax.devices())
        )

    _device = gpu_devices[0]
    _enabled = True
    _initialized = True

    if _verbose:
        print("[gpu_backend] GPU enabled")
        print("[gpu_backend]   default_backend = {0}".format(jax.default_backend()))
        print(
            "[gpu_backend]   jax_enable_x64  = {0}".format(
                bool(jax.config.jax_enable_x64)
            )
        )
        print("[gpu_backend]   devices         = {0}".format(list(jax.devices())))
        print("[gpu_backend]   device in use   = {0}".format(_device))

    return _enabled


def get_device():
    """The GPU device in use, or ``None`` when the GPU is disabled/unavailable."""
    if not _enabled:
        return None
    jax, _ = _import_jax()
    try:
        return list(jax.devices("gpu"))[0]
    except Exception:
        return None


def to_device(arr):
    """Move an array to the GPU when enabled; otherwise return it unchanged."""
    if not _enabled:
        return arr
    jax, _ = _import_jax()
    return jax.device_put(arr, _device)


def jitted(fn, *args, **kwargs):
    """Thin ``jax.jit`` wrapper.

    Returns ``fn`` unchanged when the GPU is disabled, so callers may use it
    unconditionally without paying a JIT cost on the CPU path.
    """
    if not _enabled:
        return fn
    jax, _ = _import_jax()
    return jax.jit(fn, *args, **kwargs)


def cached_jit(fn, key):
    """F1 (Fusion): ``jax.jit`` whose compiled graph is **cached by a
    structural signature** ``key`` and returned unchanged on later calls.

    Callers build one fused function per structural variant (Python flags,
    never per-evaluation values -- values are always *arguments*) and reuse
    the same compiled graph every evaluation: no re-tracing, no per-eval
    host<->device round trips between the chain's stages.

    Returns ``fn`` unchanged when the GPU is disabled (CPU path untouched).
    """
    if not _enabled:
        return fn
    cached = _GRAPH_CACHE.get(key)
    if cached is None:
        jax, _ = _import_jax()
        cached = jax.jit(fn)
        _GRAPH_CACHE[key] = cached
    return cached


def _free_gpu_mem_mib():
    """Best-effort free VRAM (MiB) of the first visible GPU via nvidia-smi;
    ``None`` when nvidia-smi is unavailable or the output is unparseable."""
    try:
        import subprocess

        out = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=memory.free",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            timeout=15,
        ).stdout.splitlines()
        s = out[0].strip() if out else ""
        if s.isdigit():
            return int(s)
    except Exception:
        pass
    return None


def _allocator_defaults():
    """(PREALLOCATE, MEM_FRACTION) chosen by current free VRAM (F2)."""
    free = _free_gpu_mem_mib()
    if free is not None and free < 10 * 1024:
        # Card nearly full (shared card, e.g. ollama): a preallocated 0.30
        # arena would OOM at context init, so stay on-demand and small.
        return "false", "0.10"
    # Free (or unknown): deterministic preallocated arena so large GPs can
    # allocate without on-demand fragmentation (the N=10k OOM root cause).
    return "true", "0.30"


def xla_allocator_envs():
    """Resolve the XLA allocator env defaults for a GPU run (F2).

    Precedence per value: ``PYORBIT_XLA_PREALLOCATE`` /
    ``PYORBIT_XLA_MEM_FRACTION`` overrides, then any already-set
    ``XLA_PYTHON_CLIENT_*`` env, then the free-VRAM-based defaults
    (``true``/``0.30``, or ``false``/``0.10`` when < ~10 GiB is free).
    Imported only for ``os`` + ``subprocess`` -- no JAX, so it is safe to
    call from a CPU-only process or before ``jaxlib`` is imported.
    """
    pre = os.environ.get(
        "PYORBIT_XLA_PREALLOCATE", os.environ.get("XLA_PYTHON_CLIENT_PREALLOCATE")
    )
    frac = os.environ.get(
        "PYORBIT_XLA_MEM_FRACTION", os.environ.get("XLA_PYTHON_CLIENT_MEM_FRACTION")
    )
    if pre is None or frac is None:
        dpre, dfrac = _allocator_defaults()
        if pre is None:
            pre = dpre
        if frac is None:
            frac = dfrac
    envs = {}
    if pre is not None:
        envs["XLA_PYTHON_CLIENT_PREALLOCATE"] = pre
    if frac is not None:
        envs["XLA_PYTHON_CLIENT_MEM_FRACTION"] = frac
    return envs


def _class_allowed(instance, allowed):
    """True if ``instance``'s class is in ``allowed`` or subclasses an allowed base."""
    if not allowed:
        return True  # no allowlist defined for this namespace -> do not reject
    if type(instance) in allowed:
        return True
    try:
        return isinstance(instance, tuple(allowed))
    except TypeError:
        return False


def _redirect_message(cls, redirect):
    """Actionable "use this JAX model instead" message for a rejected class."""
    if cls in redirect:
        return redirect[cls]
    n = cls.__name__.lower()
    if "celerite" in n or "spleaf" in n:
        return "use the tinygp_* equivalent"
    if "batman" in n:
        return "use pytransit_transit"
    if "dynamical" in n:
        return "the TRADES integrator is not JAX-portable in Round 1 (CPU-only)"
    if "rossiter" in n or "mclaughlin" in n:
        return "no JAX backend in Round 1 (CPU-only)"
    if "spiderman" in n:
        return "not JAX-native in Round 1 (CPU-only)"
    if (
        "gpgaussianprocess" in n
        or "gpframework" in n
        or "pyaneti" in n
        or "quasiperiodic" in n
    ) and "tiny" not in n:
        return "use a tinygp_* Gaussian process"
    return "not JAX-native in Round 1 (CPU-only)"


def gate(models, common_models=None):
    """Reject any model in the fit that is not JAX-native, failing fast.

    ``models`` and ``common_models`` are the ``{name: model_object}`` dicts owned
    by the ModelContainer. Membership is checked by the object's **class**
    (``type(obj)``), which is robust to PyORBIT's internal model naming: the keys
    of ``mc.models`` are internal names (e.g. ``radial_velocities_planet_1``),
    not the YAML keyword, whereas the class is what determines the backend.

    This is a no-op while the GPU is disabled, so the CPU path is untouched.
    """
    if not _enabled:
        return

    # Lazy import: pyorbit.model_definitions imports every model class (and thus,
    # via tinygp, jax). Importing it here -- at gate time, never at module import
    # time -- avoids a circular import (model_definitions -> keplerian_modelling
    # -> gpu_backend -> model_definitions). At this point model_definitions is
    # already fully loaded in sys.modules, so this is a cheap attribute fetch.
    import pyorbit.model_definitions as md

    allowed_models = set(getattr(md, "GPU_ALLOWED_MODELS", ()))
    allowed_common = set(getattr(md, "GPU_ALLOWED_COMMON", ()))
    redirect = getattr(md, "GPU_REDIRECT", {})

    problems = []
    if models is not None:
        for name, model in models.items():
            cls = type(model)
            if not _class_allowed(model, allowed_models):
                problems.append((name, cls, _redirect_message(cls, redirect)))
    if common_models is not None:
        for name, common in common_models.items():
            cls = type(common)
            if not _class_allowed(common, allowed_common):
                problems.append((name, cls, _redirect_message(cls, redirect)))

    if problems:
        lines = [
            "gpu: true was requested, but the following objects in the fit are not JAX-native:"
        ]
        for name, cls, msg in problems:
            lines.append("  * {0} ({1}): {2}".format(name, cls.__name__, msg))
        lines.append("Remove them from the fit, or set `gpu: false` to run on the CPU.")
        raise GPUSupportedModelError("\n".join(lines))
