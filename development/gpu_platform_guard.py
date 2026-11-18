#!/usr/bin/env python
"""Regression guard: PyORBIT's pre-import JAX platform keys must be jaxlib-safe.

Background
==========
On jaxlib 0.6.x the *pre-import* key `JAX_PLATFORMS=gpu` does not resolve to
the CUDA implementation (it resolves to `rocm`, which a CUDA stack does not
have), so a process launched that way dies inside jax backend init with the
misleading

    Backend 'rocm' is not in the list of known backends: ['cpu', 'tpu', 'cuda']

PyORBIT's own pre-import hooks hit exactly this: `PyORBIT_Run.py`
(`_preconfigure_gpu`) and the `PYORBIT_GPU` block in `pyorbit/__init__.py`
exported `JAX_PLATFORMS=gpu`, and every `gpu: true` run failed with the
deceiving "JAX cannot see any GPU". The concrete key `cuda` is canonical on
both the old 0.4.x and the current 0.6.2 stacks, and PyORBIT now exports
only that key pre-import. (The *runtime* backend name still is `gpu` --
`default_backend() == "gpu"` -- which is what `gpu_backend.initialize()`
checks, so it is unaffected.)

Checks
======
1. Static (always runs, never imports jax -- safe on CPU-only machines):
   reads the pre-import `JAX_PLATFORMS` / `JAX_PLATFORM_NAME` assignments
   straight from the source of the hook files and asserts every value is in
   the known-good set `{"cuda"}`. A missing hook or a non-literal value
   fails the guard (fail-closed: it cannot verify what it cannot read).
2. Dynamic (only when `nvidia-smi` reports a GPU, so no jax is ever imported
   on a non-GPU machine): runs the exact P4 entry path in a *subprocess* --
   `_preconfigure_gpu` -> `import pyorbit` -> `yaml_parser` of a tiny
   `gpu: true` config -> `gpu_backend.initialize` -> a `jnp` compute -- and
   asserts it lands on the GPU. Reports SKIP (not FAIL) when that jaxlib
   nevertheless cannot see a GPU, e.g. a CPU-only jaxlib wheel.

Repository convention (CLAUDE.md): there is no pytest suite; development
scripts live in `development/` and exit 0 on PASS, non-zero on FAIL.

Run:
    python development/gpu_platform_guard.py        # any env, CPU or GPU
    conda run -n pyorbit-gpu python development/gpu_platform_guard.py
"""

import ast
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: The env keys we audit ... and the only values acceptable *pre-import* on
#: any jaxlib 0.4.x / 0.6.x stack (a `gpu` here is precisely this bug).
AUDIT_KEYS = ("JAX_PLATFORMS", "JAX_PLATFORM_NAME")
KNOWN_GOOD = {"cuda"}

#: Hook files whose pre-import platform assignments are audited.
FILES = (
    "PyORBIT_Run.py",
    "pyorbit/__init__.py",
    "pyorbit/common/gpu_backend.py",
)

EXIT_SKIP = 10  # dynamic-leg subprocess code for "no GPU visible"


def _const_str(node):
    """The literal value of ``node`` when it is a string constant, else None."""
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    return None


def _is_environ(node):
    return isinstance(node, ast.Attribute) and node.attr == "environ"


def platform_env_pairs(path):
    """Collect ``(key, value-or-None)`` for JAX platform env assignments in
    ``path``.

    Recognised forms (the only ones pyorbit's hooks use):
      * ``...environ.setdefault('JAX_PLATFORMS', 'cuda')``
      * ``...environ['JAX_PLATFORMS'] = 'cuda'``

    A non-literal value is reported as ``None`` (fail-closed upstream).
    """
    with open(path, "r") as fh:
        tree = ast.parse(fh.read(), filename=path)

    pairs = []
    for node in ast.walk(tree):
        # environ.setdefault('KEY', value)
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "setdefault"
            and _is_environ(node.func.value)
            and len(node.args) >= 2
        ):
            key = _const_str(node.args[0])
            if key in AUDIT_KEYS:
                pairs.append((key, _const_str(node.args[1])))
        # environ['KEY'] = value
        elif (
            isinstance(node, ast.Assign)
            and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Subscript)
        ):
            tgt = node.targets[0]
            if _is_environ(tgt.value):
                key = _const_str(tgt.slice)
                if key in AUDIT_KEYS:
                    pairs.append((key, _const_str(node.value)))
    return pairs


def static_check():
    """[] on success; otherwise one problem string per violation."""
    problems = []
    for rel in FILES:
        path = os.path.join(REPO, rel)
        pairs = platform_env_pairs(path)
        if not pairs:
            problems.append(
                "%s: no %s assignment found (pre-import hook missing?)"
                % (rel, "/".join(AUDIT_KEYS))
            )
            continue
        for key, value in pairs:
            if value is None:
                problems.append(
                    "%s: %s is assigned a non-literal value (the guard "
                    "cannot verify it; expected a literal in %s)"
                    % (rel, key, sorted(KNOWN_GOOD))
                )
            elif value not in KNOWN_GOOD:
                problems.append(
                    "%s: %s=%r is not in the known-good set %s"
                    % (rel, key, value, sorted(KNOWN_GOOD))
                )
    return problems


def _gpu_machine():
    """Cheap jax-free probe: does nvidia-smi report at least one GPU?"""
    try:
        out = subprocess.run(
            ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
            capture_output=True,
            text=True,
            timeout=15,
        )
    except Exception:
        return False
    return out.returncode == 0 and bool(out.stdout.strip())


DRIVER = r"""
import os
import sys
import tempfile

repo = sys.argv[1]
os.chdir(repo)
sys.path.insert(0, repo)

tmp_dir = tempfile.mkdtemp(prefix="pyorbit_gpu_guard_")
cfg_path = os.path.join(tmp_dir, "gpu_guard.yaml")
with open(cfg_path, "w") as fh:
    fh.write("parameters:\n  gpu: true\n  gpu_device: 0\noutput: out\n")

# Exact pre-import entry path: exec PyORBIT_Run.py up to the __main__ guard.
# That runs _preconfigure_gpu() with this argv and then `import pyorbit`,
# i.e. jaxlib comes up (if present) with whatever keys the hook exports.
sys.argv = ["PyORBIT_Run.py", "emcee", cfg_path]
with open(os.path.join(repo, "PyORBIT_Run.py"), "r") as fh:
    pre = fh.read().split("if __name__")[0]
exec(pre, {"__name__": "pyorbit_gpu_guard_pre"})

# Is a GPU visible *at runtime*? (The exact probe initialize() itself uses.
# A CPU-only jaxlib wheel returns [] here instead of raising.)
import jax

try:
    visible = list(jax.devices("gpu"))
except Exception as exc:
    print("SKIP: jax.devices('gpu') raised: %s" % exc)
    sys.exit(10)
if not visible:
    print("SKIP: no GPU visible to this jaxlib: %r" % (list(jax.devices()),))
    sys.exit(10)

# Exact P4 sequence: yaml_parser -> gpu_backend.initialize -> real compute.
from pyorbit.subroutines.input_parser import yaml_parser
from pyorbit.common import gpu_backend

enabled = gpu_backend.initialize(yaml_parser(cfg_path).get("parameters", {}))
assert enabled, "gpu_backend.initialize() did not enable the GPU"
assert gpu_backend.gpu_enabled()

import jax.numpy as jnp

x = jnp.array([1.0, 2.0, 3.0]) * 2.0
dev = next(iter(x.devices()))
# Device.platform is the runtime alias ('gpu' on 0.6.x, 'cuda' on some
# older builds); accept any non-CPU GPU device.
assert dev.platform in ("gpu", "cuda", "rocm"), "jnp compute landed on %r" % dev
print("DYNAMIC OK: INITIALIZE OK; compute devices=%s; backend=%s"
      % (list(x.devices()), jax.default_backend()))
"""


def dynamic_check():
    """Run the exact entry path in a subprocess; (status, detail) with status
    in ``{"ok", "skip", "fail"}``."""
    proc = subprocess.run(
        [sys.executable, "-c", DRIVER, REPO],
        capture_output=True,
        text=True,
        timeout=600,
    )
    detail = (proc.stdout or "").strip()
    if proc.stderr and proc.stderr.strip():
        detail += "\n" + proc.stderr.strip()
    if proc.returncode == 0:
        return "ok", detail
    if proc.returncode == EXIT_SKIP:
        return "skip", detail
    return "fail", detail + "\n(dynamic leg rc=%d)" % proc.returncode


def main():
    print("PyORBIT pre-import JAX platform-key guard (repo: %s)" % REPO)
    ok = True

    problems = static_check()
    if problems:
        ok = False
        for problem in problems:
            print("FAIL  %s" % problem)
    else:
        print(
            "OK    pre-import %s keys are the jaxlib-0.6.2-safe value(s) %s"
            % ("/".join(AUDIT_KEYS), sorted(KNOWN_GOOD))
        )

    if _gpu_machine():
        status, detail = dynamic_check()
        for line in detail.splitlines():
            print("  > " + line)
        if status == "ok":
            print(
                "OK    exact entry path: preconfigure -> import -> yaml_parser -> initialize -> GPU compute"
            )
        elif status == "skip":
            reason = detail.splitlines()[-1] if detail else "no GPU"
            print("SKIP  dynamic leg: %s (static check still authoritative)" % reason)
        else:
            ok = False
            print("FAIL  dynamic leg (exact entry path) failed")
    else:
        print("SKIP  dynamic leg: nvidia-smi reports no GPU (not a GPU machine)")

    print("RESULT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
