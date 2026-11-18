(gpu_usage)=

# Using the GPU

`PyORBIT` can run its per-sample compute (Keplerian solver, transit light curve, activity GPs, Gaussian processes, trends) on an NVIDIA GPU via a JAX/CUDA backend. This page explains how the GPU mode works, what it can and cannot run, the parameter knobs, the measured CPU↔GPU numerical agreement, the known behaviors and recent fixes, a complete worked example with the data it uses, and how to verify and troubleshoot everything.

:::{admonition} CPU-only users
:class: note
If you never set `gpu: true`, none of this applies to you: `gpu: false` is a **complete no-op** — the gate (`pyorbit/common/gpu_backend.py::gate`) returns immediately, no JAX platform is forced, and the code path is identical to a CPU-only install. GPU mode is fully opt-in.
:::

## How the GPU gate works

The whole feature is driven by one YAML flag (plus two optional knobs):

```yaml
parameters:
  gpu: true                # default: false
  gpu_device: 0            # optional, which CUDA device to target
  gpu_verbose: true        # optional, extra setup diagnostics
```

When `gpu: true`, `PyORBIT_Run.py` **before JAX is imported** sets `JAX_ENABLE_X64=true`, the JAX platform, and (if `gpu_device` is given and unset) `CUDA_VISIBLE_DEVICES`. Then the module `pyorbit/common/gpu_backend.py` (1) retargets JAX to the CUDA device you select, forces **float64** everywhere (`jax_enable_x64` — there is **no float32 path**), and (2) runs a **gate** over every model and common object in the fit at setup time:

- models whose class is in `GPU_ALLOWED_MODELS` (`pyorbit/model_definitions.py`) pass;
- everything else **fails fast** with a `GPUSupportedModelError` that lists each offending object and its guidance. The guidance is honest: it either names an **allowlisted** class (for example *"in this fork the GPU transit backend is Batman_Transit (already permitted on gpu: true)"*) or says plainly *"no GPU backend in this fork; run with `gpu: false`"*. (Earlier versions pointed at classes that are themselves gated; that is no longer the case.)

There is deliberately **no silent CPU fallback**: if part of your model chain cannot run on the GPU, the run stops at setup with an actionable message rather than quietly running slower on the CPU. This invariant is maintained by design.

One precise qualification to the "no-op" wording: `gpu: false` means **no PyORBIT GPU forcing** (the gate returns immediately and no JAX platform is forced); it does *not* pin JAX to the CPU. In a CUDA-jaxlib env (e.g. `pyorbit-gpu`), JAX-native models (`tinygp_*` GPs) default to the GPU backend even under bare JAX, so a same-env A/B "CPU" leg runs tinygp on the GPU while PyORBIT dispatch (Kepler solver, transit, trend) stays numpy.

:::{admonition} Why fail fast?
:class: tip
A "GPU" run that secretly evaluates a numba model on the CPU is the worst of both worlds (you pay the JAX overhead and lose the GPU). The gate makes an incompatible fit impossible to start, so you always know exactly where your computation is running.
:::

## Installing the GPU stack

The base package keeps wide JAX bounds so a CPU-only install always resolves. The `[gpu]` extra pins a CUDA-capable stack that actually drives a GPU:

```{code} bash
# CPU-only (default)
pip install -e .

# GPU-enabled (requires a recent NVIDIA CUDA driver)
pip install -e .[gpu]
```

`.[gpu]` pins (see `pyproject.toml`): `jax[cuda12]==0.6.2`, `ml-dtypes<0.6`, `cupy-cuda12x>=13,<14`, `jaxoplanet==0.1.0`, `spleaf==2.1.20`, `batman-package==2.5.3` (batman is the CPU reference kernel).

:::{admonition} Not the requirements file
:class: note
`extra_requirements.txt` pins the **older CPU stack**; it is *not* the GPU install path. The `[gpu]` extra is.
:::

:::{admonition} Blackwell (RTX PRO 4000/5000/6000)
:class: warning
The *CPU-only* `jaxlib` build (e.g. `jaxlib 0.4.31`) that a plain `pip install -e .` can resolve **cannot** drive a Blackwell GPU (compute capability 12.x). If you see `gpu: true was requested ... no GPU is reachable` or a CUDA driver/library error, you are almost certainly running the CPU-only environment. Use the environment installed with `.[gpu]` (or your equivalent `conda` env with the CUDA JAX stack) for any `gpu: true` run.
:::

Verify the environment before running a fit:

```{code} bash
nvidia-smi                                            # driver + visible VRAM
python -c "import jax; print(jax.devices())"          # must list a CUDA device
```

## Which models run on the GPU

The authoritative list is `GPU_ALLOWED_MODELS` in `pyorbit/model_definitions.py` (read it for the exact classes). In YAML keywords:

| YAML keyword(s) | Classes (allowlist) | GPU backend |
|---|---|---|
| RV Keplerian models — including apodized RV and transit-time variants | `RVkeplerian`, `ApodizedRVkeplerian`, `TransitTimeKeplerian` | JAX-native Keplerian solver (`pyorbit/subroutines/kepler_exo_jax`) |
| `batman_transit` | `Batman_Transit` | jaxoplanet circular transit (`pyorbit/subroutines/transit_gpu`). **Gate constraint: circular orbit + quadratic limb darkening only** — anything else (TTV, eccentric, non-quadratic LD) is rejected |
| `spleaf_esp`, `spleaf_multidimensional_esp` (and the `spleaf_exponentialsineperiodic` aliases), single- and multi-dataset | `SPLEAF_ESP`, `SPLEAF_Multidimensional_ESP` | splEaF ESP Gaussian process on a CUDA kernel (`pyorbit/subroutines/esp_gpu`) — functional in the official `[gpu]` environment since the fixes described in [Known behaviors & recent fixes](#known-behaviors--recent-fixes) |
| tinygp GPs — `tinygp_quasiperiodic`, `_derivative`, `_cosine`, `_squaredexponential`, `tinygp_sho`, `tinygp_matern32`, `tinygp_multidimensional_*`, `tinygp_multiquasiperiodic_trainedsquaredexponential` | the eleven `TinyGP...` classes | tinygp-native GPs (JAX-native, no dispatch needed) |
| `polynomial_trend` (per-dataset) | `PolynomialTrend` | JAX-native Horner evaluation |
| `offset`, `jitter`, dilution / normalization / subset variants | `Offset`, `Jitter`, `DilutionFactor`, `LocalDilutionFactor`, `NormalizationFactor`, `LocalNormalizationFactor`, `SubsetNormalizationFactor` | JAX-safe scalar arithmetic |

All common objects (planets, star, limb darkening, activity, ...) are parameter holders only (`GPU_ALLOWED_COMMON`) and are always allowed.

Models that are **not** GPU-native are rejected with explicit guidance (from the `GPU_REDIRECT` map):

| You wrote | Gate says |
|---|---|
| `quasiperiodic` (george) | use `TinyGaussianProcess_QuasiPeriodicActivity` (tinygp) instead |
| `gp_matern32` (celerite2) | use `TinyGaussianProcess_Matern32` (tinygp) instead |
| `spleaf_esp_slow` (+ `_exponentialsineperiodic_slow` aliases), `spleaf_multidimensional_esp_slow`/`_devel` | no GPU backend in this fork; run with `gpu: false` |
| `batman_transit` TTV / subset / eclipse variants | no GPU backend in this fork; run with `gpu: false` |
| RML fits (`RossiterMcLaughlin_*`), TRADES dynamical orbits, `Spiderman_Thermal` starry phase curves, `shared/local/subset_polynomial_trend`, sinusoids, correlation models, CHEOPS detrending, `SpectralRotation*` | no GPU backend in this fork; run with `gpu: false` |

## GPU parameters

All live in the YAML `parameters:` section:

```yaml
parameters:
  Tref: 2460000.00
  cpu_threads: 16

  # ---- GPU selection ----
  gpu: true                # enable the JAX/CUDA backend (default: false)
  gpu_device: 0            # which CUDA device to target
  gpu_verbose: true        # print JAX/retargeting diagnostics at setup

  # ---- GPU worker pool (dynesty only, see below) ----
  gpu_workers: 8           # K concurrent dynesty live-point workers (default: 1 = single process)
  gpu_vram_fraction: 0.1   # per-process VRAM cap, as a fraction of the card (default: see below)
```

- `gpu` — the master switch; `gpu: false` is a complete no-op.
- `gpu_device` — index of the CUDA device (JAX devices are enumerated at import time).
- `gpu_verbose` — extra setup diagnostics.
- `gpu_workers` — how many **parallel worker processes** dynesty uses to evaluate its independent live points. Default `1` is the single-process behaviour.
- `gpu_vram_fraction` — caps each process's JAX memory at `fraction × card_VRAM` (`XLA_PYTHON_CLIENT_MEM_FRACTION`). Default when unset: `min(0.30, max(0.05, 0.90 / (gpu_workers + 1)))`, so with `gpu_workers: 8` it defaults to `0.10` per process.

## Samplers & VRAM

**emcee.** Under `gpu: true`, emcee is **forced single-process** (`pyorbit/samplers/pyorbit_emcee.py`): forking a worker after the parent has used CUDA is unsafe, and JAX already parallelizes on the GPU. You will see `GPU enabled: single-process (use_threading_pool is ignored).`

**dynesty worker pool.** A transit+GP likelihood is small (tens to thousands of data points, a few dozen parameters). Evaluated one at a time, each `mc(theta)` call is dominated by host↔device round-trips and Python control flow, so the GPU sits idle between evaluations. The fix is to keep the GPU busy with **K concurrent evaluations** that share one card.

When `gpu: true` and `gpu_workers > 1`, `pyorbit_dynesty` hands dynesty's live points to a pool of `gpu_workers` **spawned** worker processes (not forked — a CUDA context cannot be inherited across `fork()`):

- each worker is a fresh interpreter that rebuilds an equivalent `ModelContainerDynesty` from the pickled `config_in` and its own CUDA context;
- the master process keeps one context for the prior transform (dynesty's `use_pool={'prior_transform': False}`), so the prior always runs on the master's copy;
- `XLA_PYTHON_CLIENT_PREALLOCATE=false` and `XLA_PYTHON_CLIENT_MEM_FRACTION=<gpu_vram_fraction>` are exported before the workers start (`pyorbit/samplers/pyorbit_dynesty.py`), so every process (master included) is VRAM-capped;
- checkpoint/restore is unchanged (dynesty strips the pool on pickling and re-attaches it on `restore`).

Precedence: **GPU spawn pool** (`gpu: true`, `gpu_workers > 1`) → CPU fork pool (`use_threading_pool`, CPU-only runs) → single process.

### VRAM budget & shared cards

With `gpu_workers: 8` and the default `0.10` per-process fraction, the 8-worker budget on this **32.6 GiB** card is `8` workers × (0.10 × 32.6 GiB ≈ 3.26 GiB) ≈ **26.1 GiB that must be free** before the run starts (the master holds its own 0.10 share in addition).

:::{admonition} Known limitation: high worker counts on shared cards
:class: warning
A documented **8-worker** dynesty leg needs ≈ **29 GiB of free VRAM**. On a *shared* ~32.6 GiB card with, say, ~5.4 GiB free (e.g. a local LLM server is resident), that leg is simply infeasible — the run cannot start. This is a known environmental limitation, not a bug: the authoritative CPU↔GPU equivalence evidence for this fork is the **unit-level parity suite and B1–B6 probes** (below), which run single-context. On a small or shared card, lower the budget instead (fewer workers and/or a smaller fraction), e.g. `gpu_workers: 4`, `gpu_vram_fraction: 0.15`, and check first:

```{code} bash
nvidia-smi   # look for other processes (e.g. an LLM server like ollama) holding VRAM
```
:::

## Worked example: TOI-1554 (Kepler-63 b)

A complete, self-contained CPU/GPU A/B fit lives in [`development/examples/kepler63b_gpu/`](https://github.com/tiziano1590/PyORBIT-gpu/tree/main/development/examples/kepler63b_gpu). It fits the radial velocities plus two activity indicators (Hα equivalent width, line core depth) with a multi-dataset quasi-periodic activity GP, and two independent light curves (Kepler + TESS) with a transit plus a low-order detrend — the classic "free transit depth" problem.

### The data

| File | Contents | Points | Epochs (BJD) |
|---|---|---|---|
| `RVS_PyORBIT.txt` | epoch, RV (m/s), error, `jitter_flag`, `offset_flag`, `subset_flag` | 62 | ~2460040–2460280 |
| `Halpha_PyORBIT.txt` | epoch, Hα equivalent width, error, flags | 62 | ~2460040–2460280 |
| `EW_PyORBIT.txt` | epoch, line-core depth, error, flags | 62 | ~2460040–2460280 |
| `KEPLER63_KEPLER.DAT` | epoch, relative flux, error, flags (normalized, activity-dominated) | 56,191 | ~2455002–2456424 |
| `KEPLER63_TESS.DAT` | epoch, relative flux, error, flags | 6,484 | ~2458683–2460533 |

### The model chain

- **Activity:** `tinygp_multidimensional_quasiperiodic` on RV + Hα + EW, with `hyperparameters_condition: True` and `rotation_decay_condition: True`, per-dataset `rot_amp`/`con_amp` (and a derivative term on RV/EW). This is the JAX-native GP — on the CPU leg the same config works because the model is JAX-native on CPU too (the CPU/GPU difference is the *backend*, not the math).
- **Transit:** `batman_transit` on each light curve, **circular orbit + quadratic limb darkening** (the two constraints the GPU transit kernel imposes), free `R_Rs` ∈ [0.01, 0.30]. On `gpu: true` this evaluates through the JAX `circular_transit` kernel; on `gpu: false` through numba `batman` — the same physical model.
- **Detrend:** one `polynomial_trend` (order 2) per light curve. A deliberate detrend choice: a QP activity GP over the full 62,675-point light curves would need a ~29 GiB kernel (plus two `jax.grad` passes) — more than the 32 GiB card — and an infeasible O(N³) solve on CPU; a quadratic trend removes the same slow activity drift at O(N) cost and isolates the transit depth.
- **Priors:** all free/internal parameters carry **wide** uniform priors (broadly non-informative), so the data — not tight priors — sets the result.

*Provenance: light-curve size (62,675 points), per-worker ≈3.6 GiB and ≈29 GiB totals, and ≈1.6 GiB Python overhead are inherited from the original upstream GPU documentation (the author's measurements on their machine); they were not re-measured in this fork's verification — the 8-worker A/B was infeasible at ~5.4 GiB free VRAM.*

### The parameter files

Two byte-identical configs, differing **only** in the `gpu` flag:

- `TOI1554_CPU_dynesty_phot_WIDE_Trend.yaml` — `gpu: false` (CPU reference leg)
- `TOI1554_GPU_dynesty_phot_WIDE_Trend.yaml` — `gpu: true`, `gpu_device: 0`, `gpu_workers: 8`, `gpu_vram_fraction: 0.1`

Solver settings (both legs): `nlive: 200`, `bound: multi`, `sample: rwalk`, `dlogz: 0.01`, `maxcall: 40000`, `use_threading_pool: False`, `recenter_bounds: True`.

### Running it

```{code} bash
cd development/examples/kepler63b_gpu

# GPU leg (run in the GPU environment; make sure ~29 GiB of VRAM is free first)
pyorbit_run dynesty TOI1554_GPU_dynesty_phot_WIDE_Trend.yaml
pyorbit_results dynesty TOI1554_GPU_dynesty_phot_WIDE_Trend.yaml -all

# CPU reference leg (same model, numba transit; single process)
pyorbit_run dynesty TOI1554_CPU_dynesty_phot_WIDE_Trend.yaml
pyorbit_results dynesty TOI1554_CPU_dynesty_phot_WIDE_Trend.yaml -all
```

The two legs are a true A/B: same data, same priors, same solver, same model classes — the only differences are the compute backend (numba vs JAX/CUDA) and the dynesty worker pool on the GPU leg.

### Monitoring the GPU while it runs

`check_gpu_usage.sh` (in the same folder) is a read-only live monitor:

```{code} bash
./check_gpu_usage.sh            # 1 s refresh, Ctrl-C to stop
./check_gpu_usage.sh --once     # single snapshot
./check_gpu_usage.sh -n 10      # 10 samples and exit
./check_gpu_usage.sh --log run.log
```

It prints per-process VRAM (so you can see the workers each staying under their cap) and flags any other process holding VRAM (e.g. a local LLM server) that is stealing from the run.

## Measured CPU↔GPU numerical agreement

The following numbers were **measured** in the parity campaign (logs in `/tmp/opencode/pyorbit_gpu_verify/`); cite them as-is. Where the legs use *different algorithms* (noted per row), "agreement" means the surfaces agree at the stated band for the same physics — it is **not** a blanket "float64 round-off" across every surface.

| Surface | Probe | Measured | PASS band |
|---|---|---|---|
| Full parity suite, end-to-end `mc(θ)` (Kepler RV only — `radial_velocities`/`keplerian`, n = 60), GPU env | `development/gpu_cpu_parity.py` (section 5) | `cpu = gpu = -274.4690725832`, max\|diff\| = **0.0** — bit-identical | rtol = 1e-6, atol = 1e-8 |
| B1: transit light curves, batman (adaptive quadrature) vs jaxoplanet (fixed-order quadrature) — 6561-config grid | `development/gpu_parity_b1_transit.py` | max\|diff\| = **7.794456e-09** | ≤ 1e-3 (legitimate quadrature floor) |
| B3: splEaF ESP activity GP, single- and multi-dataset incl. `con_amp`/`rot_amp` (CPU splEaF vs CUDA ESP kernel) | `development/gpu_parity_b3_esp.py` | worst max\|Δloglike\| = **1.742301e-07** | < 1e-6 → PASS |
| B4: apodized RV (jax.scipy vs scipy `norm.pdf` + full model) | `development/gpu_parity_b4_apodized.py` | worst rel = **7.897982e-13** | ≤ 1e-10 |
| B5: polynomial trend (GPU Horner vs `np.polyval`) | `development/gpu_parity_b5_trend.py` | **0.0**, 60/60 cases bit-identical | exact |
| B6: tinygp GP `log_probability` (same jax 0.6.2, `JAX_PLATFORMS=cpu` vs `cuda` legs) | `development/gpu_parity_b6_tinygp.py` | CPU −77984.59675386584 vs GPU −77984.59675386391, \|Δ\| = **1.935405e-09** | ≤ 1e-8 |

Full-suite status: in the CPU env (`pyorbit`) the suite exits 0 with **5 PASS / 1 SKIP** (section 5 SKIPs for "no GPU"); in the GPU env (`pyorbit-gpu`) it exits 0 with **8 PASS / 0 SKIP**.

**The 2026-11-18 clean run re-confirmed this table.** After the two certified blockers were fixed (the `JAX_PLATFORMS=gpu` alias bug and the fusion-chain pickle crash), the clean run re-ran the suite and the per-surface probes: all surfaces PASS with worst-Δ bit-identical to the values above (parity CPU env `pyorbit`: 5 passed, 1 skipped; parity GPU env `pyorbit-gpu`: 8 passed, 0 skipped). The full clean-run record — performance included — is in section 8 below.

## Known behaviors & recent fixes

Read this before re-hunting any of these — each was investigated and resolved during the parity campaign.

**One-pass Kepler solver (both legs, by design).** The CPU solver's loop-exit test `nd = np.sum(convd is True)` (`pyorbit/subroutines/kepler_exo.py` line 96) is a Python *identity* check on a NumPy array, so `nd` is always 0 and the `while` loop applies **exactly one** Halley update (Fulton scheme) — never iterating to the 1e-12 nominal tolerance. This quirk was **deliberately replicated** in `pyorbit/subroutines/kepler_exo_jax.py`, so parity is preserved **by construction** and the CPU and GPU legs agree to the bit. Consequence: absolute accuracy of the anomaly→RV chain is capped, with O(1e-6)-scale residuals at high eccentricity. Documented and **unfixed by design this round** — do not "fix" one leg without the other.

**ESP fix #1 — Bessel weights from cuPy (import crash).** `pyorbit/subroutines/esp_gpu.py` originally took the Bessel weights from `cupyx.scipy.special.ive`, but cuPy (official pin `cupy-cuda12x>=13,<14`) does **not** provide `ive`/`iv` — the GPU ESP path therefore crashed on import in the official `[gpu]` environment. Fixed: the weights are now computed with `scipy.special.ive` (identical to the math the CPU splEaF reference itself uses) and staged to the device as float64.

**ESP fix #2 — `-INFINITY` in the NVRTC kernel (CompileException).** The CUDA kernel used `out[0] = -INFINITY`; `nvcc` defines that C99 macro but **NVRTC** (the runtime compiler path) does not, so kernel compilation failed. Fixed: `-inf` is now passed **from the host as a kernel parameter** (`neg_inf`) — bit-identical IEEE double, semantics unchanged (log-likelihood sentinel for degenerate cases).

**Refuted: the "+1.0 transit offset" hypothesis.** Early A/B comparisons of the GPU transit suspected a convention offset of +1.0 between the JAX transit and batman (i.e. that jaxoplanet returns un-normalized flux while batman's reference is shifted by `- 1`). A dedicated 6561-configuration grid probe (`development/gpu_parity_b1_transit.py`) measured max\|diff\| = **7.794456e-09** — no offset of any magnitude exists. Root cause of the scare: jaxoplanet's public `light_curve` **already returns out-of-transit-zero**, matching batman's `- 1` convention exactly; the two implementations agree at the quadrature-floor level. Do not reintroduce a `- 1` term.

**Gate message wording (what changed).** The `GPUSupportedModelError` "no GPU reachable" message **no longer tells users to install `jax==0.4.31` wheels**; it now points to the `[gpu]` extra (`pip install pyorbit-package[gpu]`, JAX 0.6.2 CUDA stack). All `GPU_REDIRECT` messages now name only allowlisted classes or honestly say "no GPU backend in this fork; run with `gpu: false`" (previously some pointed at `pytransit_*`-style classes that are themselves gated).

## Verification recipe

Authoritative reproducibility procedure. Environment activation for interactive work: `conda activate pyorbit-gpu`.

**Full parity suite (both environments):**

```{code} bash
# CPU reference environment -- expect exit 0: 5 PASS / 1 SKIP (section 5 SKIPs: no GPU visible)
conda run -n pyorbit python development/gpu_cpu_parity.py

# GPU environment -- expect exit 0: 8 PASS / 0 SKIP
conda run -n pyorbit-gpu python development/gpu_cpu_parity.py
```

Sections that PASS in both: `env`, `kepler-parity`, `cpu-baseline`, `zero-regression`, `gate`. The GPU-only sections (`gpu-device`, `gpu-parity`, `tinygp-gpu`) additionally PASS in the GPU environment, with `gpu-parity` reporting `cpu = gpu = -274.4690725832`, max\|diff\| = 0.0 (band rtol = 1e-6 / atol = 1e-8).

**Per-surface B1–B6 probes** (run in the GPU environment; each compares CPU vs GPU internally):

```{code} bash
conda run -n pyorbit-gpu python development/gpu_parity_b1_transit.py   # 6561 configs; PASS at max|diff| <= 1e-3 (measured 7.794456e-09)
conda run -n pyorbit-gpu python development/gpu_parity_b3_esp.py        # PASS iff worst max|Δloglike| < 1e-6 (measured 1.742301e-07)
conda run -n pyorbit-gpu python development/gpu_parity_b4_apodized.py   # PASS iff worst rel <= 1e-10 (measured 7.897982e-13)
conda run -n pyorbit-gpu python development/gpu_parity_b5_trend.py      # 60/60 bit-identical (measured 0.0)
conda run -n pyorbit-gpu python development/gpu_parity_b6_tinygp.py     # PASS iff |Δ| <= 1e-8 (measured 1.935405e-09; runs JAX_PLATFORMS=cpu vs cuda legs itself)
conda run -n pyorbit-gpu python development/gpu_platform_guard.py       # permanent pre-import platform-key regression guard: expect RESULT: PASS (static leg always runs; dynamic leg runs when a GPU is visible, SKIPs on CPU-only jaxlib)
```

**Environment notes.**

- **`LD_LIBRARY_PATH`.** Ad-hoc shells (e.g. one-off `python` invocations that load NVRTC/CUDA libraries directly) may need the environment's `nvidia/*/lib` directories added to `LD_LIBRARY_PATH` (glob the `$CONDA_PREFIX/nvidia/*/lib` entries). Normal `pyorbit_run` invocations do not need this — JAX/cupy find the libraries via `site-packages`.
- **Low-VRAM environments.** Prefer `XLA_PYTHON_CLIENT_PREALLOCATE=false` plus a small `XLA_PYTHON_CLIENT_MEM_FRACTION` (the dynesty pool path already sets both before spawning workers).
- **JAX 0.6.2 platform-alias quirk (a real PyORBIT hazard, now neutralized).** On jaxlib 0.6.x the *pre-import* key `JAX_PLATFORMS=gpu` resolves to the `rocm` implementation, which a CUDA stack does not have — `import jax` in such a process dies with `Backend 'rocm' is not in the list of known backends: ['cpu', 'tpu', 'cuda']`. This is no longer "moot for pyorbit": every `gpu: true` entry path used to export exactly `JAX_PLATFORMS=gpu` from its own pre-import hooks (`PyORBIT_Run.py::_preconfigure_gpu` and the `PYORBIT_GPU` block in `pyorbit/__init__.py`) and failed with the misleading *"JAX cannot see any GPU"*. PyORBIT now uniformly exports the concrete `JAX_PLATFORMS=cuda` / `JAX_PLATFORM_NAME=cuda` pre-import — `cuda` is the canonical key on both the jaxlib 0.4.x and 0.6.2 stacks, and the *runtime* backend name stays `gpu` (`default_backend() == "gpu"`, which is what `gpu_backend.initialize()` checks). **Do not export `JAX_PLATFORMS=gpu` yourself on jaxlib 0.6.x** — use `cuda` (as probe B6 does with `JAX_PLATFORMS=cpu` vs `cuda`).

## Troubleshooting

**`GPUSupportedModelError` at setup.** Read the message: it lists each offending model/common object and, where a JAX equivalent is allowlisted, names it; otherwise it says *"no GPU backend in this fork; run with `gpu: false`"* (see the table above). Fix the YAML and re-run. If the message says JAX "cannot see any GPU", it now points at the `[gpu]` extra — install `.[gpu]` and re-check `python -c "import jax; print(jax.devices())"`.

**OOM / "Failed to allocate" / the run dies part-way.** Too many processes × too much per-process VRAM, or another CUDA consumer is present. Lower `gpu_workers` and/or `gpu_vram_fraction`, and close other CUDA applications (`nvidia-smi`). See also the shared-card warning under [Samplers & VRAM](#samplers--vram).

**`gpu: true` but no GPU is reachable / CUDA driver error.** You are in the CPU-only environment (CPU-only `jaxlib` cannot drive newer GPUs, Blackwell included) or the ad-hoc shell cannot resolve the CUDA libraries. Activate the environment installed with `.[gpu]`, and if still failing extend `LD_LIBRARY_PATH` as in [Verification recipe](#verification-recipe).

**`pyorbit_results ... -all` hangs.** This was a real bug, and it is not GPU-specific: `pyorbit_results` evaluates models in a worker pool, and `os.fork()` after JAX has initialized is a classic deadlock (JAX/XLA runs multithreaded C thread pools that a forked child cannot inherit, and a live CUDA context cannot be inherited at all). It bites `gpu: true` runs *and* CPU-only runs that use JAX-native models such as `tinygp_*`. It is fixed in this branch — when the process is multithreaded, has used the GPU, or has merely loaded JAX, the results workers now start with `spawn` (fresh interpreters) instead of `fork` (you will see a `NOTE: ... using 'spawn' instead of 'fork'` line), and the code still falls back to serial evaluation if the pool cannot be used. If you still see a hang, run `pyorbit_results` in a fresh process (it rebuilds everything from the checkpoint).

**Verifying GPU/CPU parity.** `development/gpu_cpu_parity.py` evaluates the same model at identical parameter values on both backends and reports the max absolute difference; `development/run_kepler63b_gpu_test.sh` drives the full A/B. The authoritative, per-surface evidence is the full suite plus the B1–B6 probes listed in [Verification recipe](#verification-recipe) — expect the PASS/band pairs quoted there.

## 8. Performance (measured wall time)

**The numbers in this section are the 2026-11-18 clean run** — the first clean post-fix measurement in this fork after two certified blockers were fixed (the `JAX_PLATFORMS=gpu` alias bug and the fusion-chain pickle crash). The earlier contested-round figures (shared-card performance medians measured under contention, the old emcee A/B, the old per-surface "GPU win" claims, and the break-even projections) have been **replaced** by these clean numbers and are no longer cited. Legs: the production `gpu: true` path in the `pyorbit-gpu` env (JAX 0.6.2 + cupy CUDA stack) vs the production `gpu: false` path in the `pyorbit` env — **float64 on both legs**. Surface numbers are **medians**; **speedup = CPU median / GPU median**, so a value **< 1× is a GPU loss** (launch-overhead-dominated) and a value > 1× is a GPU win. Host/device synchronisation is forced before the clock stops (JAX returns after *launching*), so production-level latency is what is measured.

### Precision (clean-run re-confirmation — all PASS)

All parity surfaces **PASS**, with worst-Δ **bit-identical to the earlier certified contested references** (quoted above in *Measured CPU↔GPU numerical agreement*; this confirms parity, not just "pass"):

| Surface | worst-Δ | band |
|---|---|---|
| b1 transit | 7.794e-09 | 1e-3 |
| b3 esp | 1.742e-07 | 1e-6 |
| b4 apodized | 7.898e-13 | 1e-10 |
| b5 trend | 0 | 1e-10 |
| b6 tinygp | 1.935e-09 | 1e-8 |

Full suite: parity CPU env (`pyorbit`) — **5 passed, 1 skipped**; parity GPU env (`pyorbit-gpu`) — **8 passed, 0 skipped**. Both previously identified blockers (the `JAX_PLATFORMS` alias, the fusion-chain pickle crash) remained fixed throughout the run.

### Model-eval surfaces (median ms)

| Surface (N) | CPU ms | GPU ms | speedup (CPU median / GPU median) |
|---|---|---|---|
| esp (1000) | 1.0 | 5.3 | 0.18× |
| esp (5000) | 4.6 | 25.7 | 0.18× |
| kepler (100) | — | — | 0.17× |
| kepler (1000) | — | — | 0.54× |
| kepler (10000) | — | — | 4.30× |
| kepler (100000) | 12.2 | 0.3 | 36.55× |
| transit (1000) | — | — | 0.27× |
| transit (10000) | — | — | 1.38× |
| trend (1000) | — | — | 0.01× |
| trend (10000) | — | — | 0.04× |

(`—` = the per-leg median is not separately recorded for that row; the speedup ratio is the record.)

**Reading:** the GPU wins only at large per-evaluation work; launch-overhead dominates at small N. The measured crossover is on the kepler solver: 0.17× (N=100) → 0.54× (N=1000) → 4.30× (N=10000) → **36.55× at N=100000** (12.2 ms CPU vs 0.3 ms GPU). ESP is flat at 0.18× across both sizes; trend stays deep in GPU loss (0.01×, 0.04×); transit is below parity at N=1000 (0.27×) and just above at N=10000 (1.38×).

### EMCEE A/B (P3) — the per-walker baseline, not a post-fix measurement

Per-walker path (`vectorize=False`), n=600: **CPU 686.84 steps/s vs GPU 48.59 steps/s → 0.071×.** Read this as the **per-walker baseline** for that path class — it is **not** a post-fix measurement of the production emcee path. The fused F-batch path (`vectorize=True`) reference — **~1.83–1.87×** — is **inherited from the certified `development/emcee_v2_ab.py` harness, not from this run**, and is structurally unavailable for GP-internal-likelihood workloads (see the P4 leg below, which is forced onto the per-walker path).

### Kepler-63b production A/B (P4) — the honest headline

| Leg | fit wall | PyDE time |
|---|---|---|
| CPU (`gpu: false`) | 25.9 s | 3.6 s |
| GPU (`gpu: true`) | 69.8 s | 9.1 s |

**≈2.7× in CPU's favor** — not a GPU speedup; this is the per-walker path class on a launch-overhead-dominated workload.

- **Both legs completed their fits (rc=0):** 500 steps, 68 walkers, 17 params, single process; configs identical except the `gpu` flag.
- **Both legs executed the per-walker emcee path**, because this workload's activity GP (GP/internal likelihood) disables the fused F-batch path:
  `emcee: F-batch vectorized lnprob not available (activity_gp: GP/internal likelihood) -- using the per-theta path`

  That is the same 0.071× path class as the P3 baseline above. The ≈2.7× CPU-favor is a **workload characteristic** (launch-overhead-dominated, small per-eval work) — **not** measured contention and **not** a broken path.
- **Posteriors mutually consistent:** P 9.452 vs 9.479 d (Δ 0.027 d); K 10.5 vs 10.6 m/s (Δ 0.1 m/s); Tc 2458690.175 vs 2458690.214 (Δ 0.039 d); all within overlapping 15–84% bands — tighter agreement than the earlier contested reference (there ΔK ≈ 1.0 m/s).
- **The GPU leg was genuinely on CUDA:** PyORBIT peaked at ~10.0 GB VRAM during the fit (dmon `sm_mean` 21.9%, `sm_max` 28%).

### Contention audit: nominal only, not measured contention

The driver's contention flag fired only nominally: ollama's HTTP endpoint was reachable (driver `CONTENDED=YES`) but it held **0 MiB VRAM** at the time — `nvidia-smi` showed **16 MiB used** on the card, with only the `gnome-shell` process (5 MiB). The dmon baseline was 16 MiB; **no second GPU tenant in any sampled window**. Plainly: **contested-flag nominal-only (ollama reachable, 0 VRAM); not measured contention.** Every performance number above is therefore uncontended.

**Where the GPU genuinely helps (this run):** large per-evaluation work — the kepler solver at N=10000 (4.30×) and N=100000 (36.55×), and transit marginally at N=10000 (1.38×). At PyORBIT's typical per-eval sizes (esp ≤ 5 k, trend, small kepler/transit), `gpu: false` remains the honest default: faster and equally correct, since parity holds on every surface (above).

### Known issue: the post-fit `pyorbit_results` step hung (NOT fit performance)

The post-fit best-effort `pyorbit_results` (plotting/results) step — 16-process pool — **hung ~118 minutes** until the driver's 2-hour timeout (**rc=124**). This is a **post-fit step issue, not evidence about fit performance** (both fits completed — see the P4 leg). Reviewer-certified sequence:

1. the first parallel step went fork→spawn via the repo guard in `pyorbit/samplers/pyorbit_getresults.py` (as designed);
2. a **spawned worker died re-initializing JAX on CUDA in the child**:
   `unable to create StreamExecutor for CUDA:0 ... CUDA_ERROR_OUT_OF_MEMORY: out of memory`, then
   `RuntimeError: Unable to initialize backend 'cuda': INTERNAL: no supported devices found for platform CUDA`
   (`pyorbit_results_gpu.log:3258-3307`);
3. the parent then blocked inside the context pool; the **defined serial fallback** (`pyorbit/samplers/pyorbit_getresults.py:154-161`) **never printed**; 6 leaked pool semaphores at death;
4. the GPU leg ran first, so the **CPU leg's results step never ran**; no `timings.json`/report artifacts were produced.

**Classification:** an **in-repo unbounded pool-wait / missing fallback trigger** atop an **environmental worker-CUDA failure**. The OOM micro-cause is unproven (no GPU-memory telemetry in that window; a foreign load — e.g. ollama re-loading its 17.7 GB model — is neither proven nor excluded). Track as a **known follow-up issue**; do **not** read it as evidence about fit performance.

**Caveats.**

- Every number in this section is from the 2026-11-18 clean run, or **explicitly attributed as inherited** (the F-batch ~1.83–1.87× reference comes from the certified `development/emcee_v2_ab.py` harness, not from this run). No contested-round figures are carried forward.
- The P4 wall times are for the per-walker emcee path class on a launch-overhead-dominated workload; they are not a general `gpu: true` vs `gpu: false` statement, and never a speedup claim for the GPU.
- Fit legs were single-process (68 walkers, 500 steps); the post-fit results step used a 16-process pool (see the known issue above).

**Repro / status of the evidence.** The clean run was driven by `development/gpu_clean_run.sh` (stages P0 contention audit → P1 precision → P2 perf surfaces → P3 emcee A/B → P4 Kepler-63b A/B → P5 stats report; artifacts default to `/tmp/opencode/pyorbit_clean_run`, with per-leg logs `P0…P4_*.log`). The precision re-confirmation above is the clean run's suite + b1/b3/b4/b5/b6 record; the known-issue worker traceback is at `pyorbit_results_gpu.log:3258-3307`.

## 9. Upstream compatibility

**Fork base:** `LucaMalavolta/PyORBIT`, branch `development`, **v12.0.0b0**; fork-point commit `0f25dc6` (2026-10-22), pulled into this fork by merge commit `9aa88ea` (parents `[2e1cf31, 0f25dc6]`).

**Adoption status (as of 2026-11-16): complete.** **All** changes upstream made since the fork point — the **6 commits to tip `cd3a707` (2026-11-09)** — are adopted into this fork:

| Upstream commit(s) | Change | File(s) |
|---|---|---|
| `092245c`, `2e56f3e` | `tc_tols` feature | `keywords_definitions.py`, `common/planets.py`, `models/dynamical_modelling.py` |
| `c742710` | `input_parser.py` robustness + the **`SharedPolynomialTrend.print_info(self)`** bug fix | `subroutines/input_parser.py`, `models/polynomial_trend.py` |
| `878ca2d` | `sinusoid.py` hygiene | `models/sinusoid.py` |
| `827fd1a` | plotting dynamical results (two `update_parameter_values` call sites) | `samplers/pyorbit_getresults.py` |
| `cd3a707` | `results_analysis.py` GP branch | `subroutines/results_analysis.py` |

**Verification evidence.**

- **4 of the 7 touched model/common files are byte-identical** to upstream: `keywords_definitions.py`, `sinusoid.py`, `input_parser.py`, `results_analysis.py`. The remaining three (`common/planets.py`, `models/dynamical_modelling.py`, `models/polynomial_trend.py`) verified **hunk-identical** — every upstream code line is present and the only residuals are out-of-hunk comment/blank-line lines (plus the fork's intentional GPU branches).
- The `print_info(self)` fix demonstrated **red→green**; trend GPU parity re-run **PASS @ 0.0** (60/60 bit-identical); `import pyorbit` OK in the CPU env; **zero conflicts** on the merge.
- **`gpu: false` zero-regression preserved:** the six upstream commits touch **no GPU-dispatch surface** — `pyorbit/subroutines/kepler_exo.py` (the CPU Kepler solver) is **byte-identical** between this fork and upstream.
- `pyorbit/samplers/pyorbit_getresults.py` also differs from upstream by exactly the fork's intentional JAX/spawn guard (see [Troubleshooting](#troubleshooting)); the upstream two-line fix in that file is adopted.

**Intentional deltas retained from upstream** (features *added* by this fork, not regressions): all GPU additions (gate, allowlist, JAX ports, dynesty worker pool), the widened pins + `[gpu]` extra in `pyproject.toml`, the 14× `scipy.matrix`→`numpy` compat fixes in the GP/spleaf models, and the `pyorbit_getresults` JAX guard.
