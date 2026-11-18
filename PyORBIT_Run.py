from __future__ import print_function
import os
import sys


def _preconfigure_gpu():
    """Path-B insurance for GPU enablement (Round 1).

    When launched as ``python PyORBIT_Run.py <sampler> <config>`` we can read the
    config and set JAX's platform environment variables *before* ``import
    pyorbit`` (which transitively imports jaxlib and locks in its default
    platform). We do that only when the config explicitly requests ``gpu: true``;
    for ``gpu: false`` (the default) this function returns without touching the
    environment, so the CPU path is byte-for-byte unchanged.

    This is a best-effort pre-import hint. The authoritative enablement, the
    GPU-reachability check and the actionable error all still happen in
    ``pyorbit.common.gpu_backend.initialize()`` after the config is parsed by
    ``pyorbit_run()``.
    """
    # argv = [script, sampler, config_file, ...]
    if len(sys.argv) < 3:
        return
    conf = sys.argv[2]
    if conf in ("-h", "--help"):
        return
    if not os.path.exists(conf):
        return
    try:
        import yaml

        with open(conf, "r") as stream:
            data = yaml.safe_load(stream)
    except Exception:
        # Can't pre-read the config for some reason; fall back to the plain CPU
        # path. yaml_parser() will report any genuine problem during the run.
        return
    if not isinstance(data, dict):
        return
    params = data.get("parameters", {})
    if not isinstance(params, dict) or not params.get("gpu", False):
        return

    # gpu: true -- advertise the GPU to jaxlib before it is imported.
    # NOTE: use the concrete 'cuda' key, not the 'gpu' alias. On jaxlib 0.6.x
    # a pre-import 'JAX_PLATFORMS=gpu' resolves to the (absent) 'rocm'
    # implementation and kills backend init; 'cuda' is correct on both the
    # 0.4.x and 0.6.x stacks. The *runtime* name remains 'gpu'
    # (default_backend() == 'gpu'), so gpu_backend.initialize() is unaffected.
    os.environ["PYORBIT_GPU"] = "1"
    os.environ.setdefault("JAX_PLATFORMS", "cuda")
    os.environ.setdefault("JAX_PLATFORM_NAME", "cuda")
    os.environ.setdefault("JAX_ENABLE_X64", "true")
    gpu_device = params.get("gpu_device", 0)
    if gpu_device is not None and "CUDA_VISIBLE_DEVICES" not in os.environ:
        os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_device)


_preconfigure_gpu()
import pyorbit  # noqa: E402  (must come after _preconfigure_gpu)
import argparse  # noqa: F401


if __name__ == "__main__":
    # print(sys.__stderr__, sys.stderr)
    pyorbit.pyorbit_run()
