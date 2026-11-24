"""Top-level, picklable entry points for dynesty's parallel worker pool.

The per-point math already lives in ``ModelContainerDynesty`` (CPU-validated).
To fill the GPU with K concurrent evaluations we evaluate dynesty's independent
live points across K worker *processes*. ``dynesty.pool`` uses
``mp.get_context('spawn')``, so each worker is a fresh interpreter that builds
its own CUDA context -- no fork-after-CUDA.

Both ``loglike`` and ``prior_transform`` are pickled into the workers (via
``initargs``), so they must be top-level functions, not closures over the
master's container. Each worker rebuilds an equivalent ``ModelContainerDynesty``
from the plain (picklable) ``config_in`` dict the first time it is called and
caches the result per process. The master reuses the already-built container
for the (in-master) prior transform so it does not open a second JAX context
and contend for VRAM.
"""
from pyorbit.classes.model_container_dynesty import ModelContainerDynesty
from pyorbit.subroutines.input_parser import pars_input

_MASTER_MC = None   # master-only: the already-built container, reused for the in-master prior transform
_MC_CACHE = None    # per-process worker: lazily-built container


def set_master_mc(mc):
    """Register the master's already-built container so the in-master prior
    transform reuses it instead of opening a second (GPU) JAX context."""
    global _MASTER_MC
    _MASTER_MC = mc


def _build_mc(config_in, input_datasets):
    mc = ModelContainerDynesty()
    pars_input(config_in, mc, input_datasets)
    mc.output_directory = './' + config_in['output'] + '/dynesty/'
    if mc.nested_sampling_parameters['shutdown_jitter']:
        for dataset in mc.dataset_dict.values():
            dataset.shutdown_jitter()
    mc.model_setup()
    mc.boundaries_setup()
    mc.initialize_logchi2()
    mc.starting_points_setup()
    return mc


def _get_mc(config_in, input_datasets):
    global _MC_CACHE
    if _MASTER_MC is not None:
        return _MASTER_MC
    if _MC_CACHE is None:
        _MC_CACHE = _build_mc(config_in, input_datasets)
    return _MC_CACHE


def pool_loglike(theta, config_in=None, input_datasets=None):
    """Log-likelihood of a single point. Top-level so it is picklable by
    reference into spawn workers."""
    return _get_mc(config_in, input_datasets).dynesty_call(theta)


def pool_priors(cube, config_in=None, input_datasets=None):
    """Unit-cube -> parameter transform. Top-level so it is picklable by
    reference into spawn workers."""
    return _get_mc(config_in, input_datasets).dynesty_priors(cube)
