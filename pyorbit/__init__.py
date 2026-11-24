#!/usr/bin/env python
# -*- coding: utf-8 -*-

# GPU enablement (Round 1) -- pre-import hook.
#
# If the launcher (``python PyORBIT_Run.py``, or a user's own script) has already
# decided this is a GPU run, it can set ``PYORBIT_GPU=1`` *before* importing
# pyorbit. We honor that here, before any of the modules below -- which
# transitively import jaxlib and lock in its default platform -- execute, so that
# jaxlib itself comes up on the GPU.
#
# This is purely optional insurance: the authoritative enablement, the
# GPU-reachability check and the actionable error all live in
# ``pyorbit.common.gpu_backend.initialize()``, which runs after the config is
# parsed. When ``PYORBIT_GPU`` is unset (the default, i.e. ``gpu: false``) this
# block does nothing and the CPU path is byte-for-byte unchanged.
import os as _pyorbit_os

if _pyorbit_os.environ.get("PYORBIT_GPU", "") not in (
    "",
    "0",
    "false",
    "False",
    "no",
    "off",
):
    # Concrete 'cuda' key, not the 'gpu' alias: on jaxlib 0.6.x a pre-import
    # 'JAX_PLATFORMS=gpu' resolves to the absent 'rocm' implementation and
    # kills backend init. 'cuda' works on both the 0.4.x and 0.6.x stacks.
    # (The runtime name stays 'gpu'; see gpu_backend.initialize().)
    _pyorbit_os.environ.setdefault("JAX_PLATFORMS", "cuda")
    _pyorbit_os.environ.setdefault("JAX_PLATFORM_NAME", "cuda")
    _pyorbit_os.environ.setdefault("JAX_ENABLE_X64", "true")

from .pyorbit_run import *
from .pyorbit_results import *
from .model_definitions import *
from .datatype_definitions import *
from .samplers.pyorbit_emcee import *
from .samplers.pyorbit_emcee_mpi import *
from .samplers.pyorbit_optimize import *
from .samplers.pyorbit_polychord import *
from .samplers.pyorbit_pyde import *
from .samplers.pyorbit_dynesty import *
from .samplers.pyorbit_dynesty_static import *
from .samplers.pyorbit_nestle import *
from .samplers.pyorbit_nautilus import *
from .samplers.pyorbit_ultranest import *
from .samplers.pyorbit_ultranest_stepsampler import *
from .samplers.pyorbit_multinest import *
from .samplers.pyorbit_getresults import *
from .subroutines.input_parser import yaml_parser


__version__ = "12.0.0b0"
