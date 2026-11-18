#!/usr/bin/env python
"""B6 tinygp value parity probe (CPU leg vs GPU leg, SAME jax 0.6.2 in both).

Motivation
==========
tinygp is JAX-native: the same tinygp build must produce the same
``GaussianProcess.log_probability(y)`` whether JAX is retargeted to the CPU or
to the GPU. The whole ``pyorbit-gpu`` parity campaign is float64, so the two
legs are compared at/near float64 with a strict value tolerance.

This probe runs the **exact same** GP evaluation twice, each as a clean
subprocess with a different JAX platform in its environment (this keeps the
pre-import ``JAX_PLATFORMS`` setup clean instead of mutating jax in-process):

* CPU leg: ``JAX_ENABLE_X64=true``  ``JAX_PLATFORMS=cpu``
* GPU leg: ``JAX_ENABLE_X64=true``  ``JAX_PLATFORMS=cuda``
            (+ ``XLA_PYTHON_CLIENT_PREALLOCATE=false``,
              ``XLA_PYTHON_CLIENT_MEM_FRACTION=0.1``)

Both legs build the same object (deterministic seed, so the data is identical)
and print ``log_probability(y)``; this script reads both back and compares.

Setup (identical both legs)
* kernel   = ``ExpSquared(scale=1.5) * ExpSineSquared(scale=2.5, gamma=0.5)``
* t        = ``np.sort(np.random.default_rng(0).uniform(0.0, 10.0, 200))``
* y        = ``np.random.default_rng(1).normal(size=200)``
* GP       = ``tinygp.GaussianProcess(kernel, t, diag=1e-3, mean=0.0)``

Classification (fixed tolerances, never relaxed)
  PASS        : |delta| log_probability <= 1e-8
  SUSPICIOUS  : 1e-8 < |delta| < 1e-6  -> STOP, report, do not fix
  BUG         : |delta| >= 1e-6        -> STOP, report, do not fix

Runs in the ``pyorbit-gpu`` conda env (has jux 0.6.2 + CUDA).

Exit codes: 0 = PASS, 2 = SUSPICIOUS, 1 = BUG / any other failure.
"""

import os
import sys
import glob
import subprocess
import tempfile

# --- fixed setup (task-specified) ------------------------------------------
N = 200
KERNEL_DESC = "ExpSquared(1.5) * ExpSineSquared(2.5, gamma=0.5)"
PASS_TOL = 1e-8
SUS_TOL = 1e-6

# The inline leg: build the exact same seeded GP and print log_probability(y).
# Runs under a plain ``python``; only jax/numpy/tinygp are needed. The value is
# emitted on the ``LOGPROB`` line so the parent can parse it back out cleanly.
LEG_CODE = r"""
import numpy as np
import jax
try:
    jax.config.update("jax_enable_x64", True)
except Exception:
    pass
t = np.sort(np.random.default_rng(0).uniform(0.0, 10.0, 200))
y = np.random.default_rng(1).normal(size=200)
import tinygp
from tinygp import kernels
gp = tinygp.GaussianProcess(
    kernels.ExpSquared(scale=1.5) * kernels.ExpSineSquared(scale=2.5, gamma=0.5),
    t,
    diag=1e-3,
    mean=0.0,
)
val = float(gp.log_probability(y))
print("LEGJAX {0} {1}".format(jax.__version__, jax.default_backend()))
print("LOGPROB " + repr(val))
"""


def _logdir():
    """Where leg subprocess output is captured (never dumped to the chat)."""
    cand = "/tmp/opencode/pyorbit_gpu_verify"
    try:
        os.makedirs(cand, exist_ok=True)
        return cand
    except Exception:
        alt = os.path.join(tempfile.gettempdir(), "pyorbit_gpu_verify")
        os.makedirs(alt, exist_ok=True)
        return alt


# jax platform name to use for each logical leg. jax 0.6.2's "gpu" *alias*
# resolves to rocm and aborts if rocm is absent, so name the concrete "cuda"
# backend for the GPU leg (the installed CUDA stack). "cpu" is already valid.
_JAX_PLATFORM = {"cpu": "cpu", "gpu": "cuda"}


def _child_env(platform):
    env = dict(os.environ)
    env["JAX_ENABLE_X64"] = "true"
    env["JAX_PLATFORMS"] = _JAX_PLATFORM[platform]
    if platform == "gpu":
        env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
        env["XLA_PYTHON_CLIENT_MEM_FRACTION"] = "0.1"

    # Extend LD_LIBRARY_PATH with the conda env lib dir and any nvidia/*/lib
    # (NVRTC/CUDA libs), the fallback that resolved in the previous round.
    prefix = os.environ.get("CONDA_PREFIX", "")
    if prefix:
        extra = [os.path.join(prefix, "lib")]
        for pydir in glob.glob(os.path.join(prefix, "lib", "python*")):
            extra += glob.glob(
                os.path.join(pydir, "site-packages", "nvidia", "*", "lib")
            )
        extra += glob.glob(os.path.join(prefix, "nvidia", "*", "lib"))
        present = [p for p in extra if os.path.isdir(p)]
        if present:
            env["LD_LIBRARY_PATH"] = os.pathsep.join(
                present + [env.get("LD_LIBRARY_PATH", "")]
            )
    return env


def _parse_logprob(stdout):
    for line in stdout.splitlines():
        if line.startswith("LOGPROB "):
            return float(line[len("LOGPROB ") :].strip())
    return None


def _run_leg(platform, logdir):
    env = _child_env(platform)
    r = subprocess.run(
        [sys.executable, "-c", LEG_CODE],
        env=env,
        capture_output=True,
        text=True,
        timeout=600,
    )
    with open(os.path.join(logdir, "b6_{0}.leg.log".format(platform)), "w") as fh:
        fh.write(
            "## exit={0}\n## --- stdout ---\n{1}\n## --- stderr ---\n{2}\n".format(
                r.returncode, r.stdout, r.stderr
            )
        )
    return r


def main():
    print("B6 tinygp value parity probe (same jax 0.6.2 both legs)")
    print("  kernel = {0}".format(KERNEL_DESC))
    print(
        "  n={0}, t=rng(0).uniform(0,10) sorted, y=rng(1).normal, diag=1e-3, mean=0.0".format(
            N
        )
    )
    print("  CPU leg: JAX_PLATFORMS=cpu     GPU leg: JAX_PLATFORMS=cuda")

    logdir = _logdir()

    cpu_r = _run_leg("cpu", logdir)
    gpu_r = _run_leg("gpu", logdir)

    cpu_lp = _parse_logprob(cpu_r.stdout)
    gpu_lp = _parse_logprob(gpu_r.stdout)

    if (
        cpu_lp is None
        or gpu_lp is None
        or cpu_r.returncode != 0
        or gpu_r.returncode != 0
    ):
        if cpu_lp is None:
            print(
                "  CPU leg FAILED to evaluate log_probability (see {0})".format(
                    os.path.join(logdir, "b6_cpu.leg.log")
                )
            )
        if gpu_lp is None:
            print(
                "  GPU leg FAILED to evaluate log_probability (see {0})".format(
                    os.path.join(logdir, "b6_gpu.leg.log")
                )
            )
        print(
            "RESULT: BUG (a leg could not be evaluated; exit codes CPU={0} GPU={1}; "
            "capture the leg logs; do not fix.)".format(
                cpu_r.returncode, gpu_r.returncode
            )
        )
        return 1

    diff = abs(cpu_lp - gpu_lp)
    cpu_ver = cpu_r.stdout.splitlines()[0] if cpu_r.stdout else ""
    gpu_ver = gpu_r.stdout.splitlines()[0] if gpu_r.stdout else ""
    print("  {0}".format(cpu_ver))
    print("  CPU leg log_probability = {0}".format(cpu_lp))
    print("  {0}".format(gpu_ver))
    print("  GPU leg log_probability = {0}".format(gpu_lp))
    print("  |delta| log_probability = {0:.6e}".format(diff))

    if diff <= PASS_TOL:
        print(
            "RESULT: PASS (|delta| log_probability={0:.6e} <= {1:.0e}: tinygp CPU==GPU "
            "at/near float64).".format(diff, PASS_TOL)
        )
        return 0
    if diff < SUS_TOL:
        print(
            "RESULT: SUSPICIOUS (|delta|={0:.6e} in ({1:.0e},{2:.0e}). "
            "STOP; capture the leg values above; do not fix.".format(
                diff, PASS_TOL, SUS_TOL
            )
        )
        return 2
    print(
        "RESULT: BUG (|delta|={0:.6e} >= {1:.0e}. "
        "STOP; capture the leg values above; do not fix.".format(diff, SUS_TOL)
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
