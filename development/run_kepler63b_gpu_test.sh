#!/usr/bin/env bash
# =============================================================================
# run_kepler63b_gpu_test.sh
#
# PyORBIT GPU-vs-CPU performance test on the Kepler-63b dataset.
#
# What it does:
#   0. Fail-fast checks: conda envs, GPU is free (ollama must be OFF),
#      and jax in the GPU env actually SEES the GPU (never a silent CPU run).
#   1. Idempotently creates the conda env `pyorbit-gpu` (clone of `pyorbit`
#      + a CUDA jax/jaxlib + tinygp) if it does not exist yet.
#   2. Creates a FRESH run dir  <TEST_FOLDER>/gpu_test_<UTC-stamp>/  with two
#      derived configs (the original TEST_FOLDER files are never modified):
#         config.gpu.yaml   parameters.gpu: true   -> output: gpu
#         config.cpu.yaml   parameters.gpu: false  -> output: cpu
#      Fresh output dirs also guarantee no checkpoint resume from old runs.
#   3. GPU run (env pyorbit-gpu, gpu:true):  PyDE + emcee, full log, plus a
#      per-second nvidia-smi utilization/memory sampler (gpu_usage.log).
#   4. CPU run (SAME env, gpu:false):       PyDE + emcee, full log.
#      Same env for both = a clean A/B of the GPU feature itself (identical
#      jax/numpy/tinygp; only the device changes).
#   5. Collects everything needed to rebuild the performance report:
#      timings.json, report_raw.md, both full logs, nvidia-smi snapshots,
#      best-effort pyorbit_results artifacts (ACF/corner/plots).
#
# Usage:
#   ./run_kepler63b_gpu_test.sh [TEST_FOLDER] [CONFIG_YAML]
#     TEST_FOLDER  (default: /home/zingales/pyorbit-tests/Kepler-63b)
#     CONFIG_YAML  (default: the single top-level *.yaml in TEST_FOLDER)
#
# BEFORE: stop ollama so the RTX PRO 4500 is not holding VRAM.
# AFTER:  restart ollama, then come back to the ollama conversation and ask to
#         rebuild the performance report from the printed run dir.
# =============================================================================
set -euo pipefail

# --- constants ---------------------------------------------------------------
TEST_FOLDER="${1:-/home/zingales/pyorbit-tests/Kepler-63b}"
CONFIG_SRC="${2:-}"
CPU_ENV="pyorbit"
GPU_ENV="pyorbit-gpu"
REPO_DIR="${PYORBIT_REPO_DIR:-/home/zingales/repos/PyORBIT-gpu}"
MIN_FREE_MIB=8192          # "full GPU" headroom required before the GPU run
SMON_PID=""

cleanup() { [ -n "${SMON_PID}" ] && kill "${SMON_PID}" 2>/dev/null || true; }
trap cleanup EXIT

die()  { echo; echo "ERROR: $*" >&2; exit 1; }
info() { echo; echo "==================================================================="; echo "==> $*"; echo "==================================================================="; }

case "${1:-}" in
  -h|--help) sed -n '2,35p' "$0"; exit 0 ;;
esac

# --- 0. fail-fast checks ------------------------------------------------------
command -v nvidia-smi >/dev/null 2>&1 || die "nvidia-smi not found (NVIDIA driver missing?)"
[ -d "$TEST_FOLDER" ] || die "test folder not found: $TEST_FOLDER"
[ -d "$REPO_DIR" ]    || die "PyORBIT repo not found: $REPO_DIR"
[ -f "$REPO_DIR/PyORBIT_Run.py" ] || die "PyORBIT_Run.py not found in $REPO_DIR"

find_conda() {
  if command -v conda >/dev/null 2>&1; then command -v conda; return 0; fi
  local p
  for p in "${HOME}/anaconda3/bin/conda" "${HOME}/miniconda3/bin/conda" \
           "${HOME}/miniforge3/bin/conda" /opt/conda/bin/conda; do
    [ -x "$p" ] && { echo "$p"; return 0; }
  done
  return 1
}
CONDA="$(find_conda)" || die "conda not found (checked PATH and common locations)"

env_exists() { "$CONDA" env list 2>/dev/null | awk '{print $1}' | grep -qx "$1"; }

run_env() { # env command...
  local env="$1"; shift
  "$CONDA" run --no-capture-output -n "$env" "$@"
}

env_exists "$CPU_ENV" || die "conda env '$CPU_ENV' not found — create it first (see memory: pyorbit-test-env)"

# GPU free? (the user stops ollama manually; we only verify and fail fast)
FREE_MIB="$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -n1 | tr -d '[:space:]' || true)"
case "${FREE_MIB:-}" in ''|*[!0-9]*) die "could not read free GPU memory via nvidia-smi";; esac
if [ "$FREE_MIB" -lt "$MIN_FREE_MIB" ]; then
  die "only ${FREE_MIB} MiB of GPU memory is free (need >= ${MIN_FREE_MIB}).
ollama is almost certainly still holding a model in VRAM.
Stop it first, e.g.:
    ollama stop <model>      # or: pkill -f ollama
then re-run this script."
fi
if pgrep -x ollama >/dev/null 2>&1; then
  echo "WARNING: an 'ollama' process is still running — if it re-loads a model it will fight PyORBIT for VRAM."
fi

# --- 1. GPU conda env (idempotent) -------------------------------------------
if ! env_exists "$GPU_ENV"; then
  info "creating conda env '$GPU_ENV' (clone of '$CPU_ENV')…"
  "$CONDA" create --yes --clone "$CPU_ENV" -n "$GPU_ENV"
  info "installing CUDA jax/jaxlib + tinygp into '$GPU_ENV' (needs network)…"
  run_env "$GPU_ENV" pip install --upgrade "jax[cuda12-pip]" "jaxlib[cuda12-pip]" tinygp \
    || run_env "$GPU_ENV" pip install --upgrade "jax[cuda13-pip]" "jaxlib[cuda13-pip]" tinygp \
    || run_env "$GPU_ENV" pip install --upgrade "jax[cuda12]" "jaxlib[cuda12]" tinygp \
    || die "could not install a CUDA jax into '$GPU_ENV' (see pip output above)"
fi

info "verifying jax in '$GPU_ENV' sees the GPU (fail-loud gate)"
run_env "$GPU_ENV" python -c 'import jax; ds=list(jax.devices()); assert any(d.platform in ("cuda","rocm") for d in ds) or len(jax.devices("gpu"))>0, "no GPU device visible to jax: %r" % (ds,); print("  jax devices:", ds)' \
  || die "jax in '$GPU_ENV' cannot see the GPU.
Fix: conda run -n $GPU_ENV pip install --upgrade 'jax[cuda12-pip]' 'jaxlib[cuda12-pip]'"

# --- 2. discover config + create fresh run dir -------------------------------
if [ -n "$CONFIG_SRC" ]; then
  [ -f "$CONFIG_SRC" ] || die "config not found: $CONFIG_SRC"
else
  cands=()
  for f in "$TEST_FOLDER"/*.yaml "$TEST_FOLDER"/*.yml; do
    [ -f "$f" ] && cands+=("$f")
  done
  if [ "${#cands[@]}" -eq 0 ]; then
    die "no top-level *.yaml/*.yml in $TEST_FOLDER — pass the config path as \$2"
  fi
  if [ "${#cands[@]}" -gt 1 ]; then
    { echo "ERROR: multiple top-level YAML configs in $TEST_FOLDER — pass the one to use as \$2:" >&2
      printf '   %s\n' "${cands[@]}" >&2; }
    exit 1
  fi
  CONFIG_SRC="${cands[0]}"
fi
echo "config: $CONFIG_SRC"

STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
RUN_DIR="$TEST_FOLDER/gpu_test_${STAMP}"
mkdir -p "$RUN_DIR"
cp -f "$CONFIG_SRC" "$RUN_DIR/config.orig.yaml"

REPORT_SRC=""
for f in "$TEST_FOLDER"/*.md; do
  if [ -f "$f" ]; then cp -f "$f" "$RUN_DIR/existing_report.md"; REPORT_SRC="$f"; break; fi
done
if [ -n "$REPORT_SRC" ]; then
  echo "existing report copied: $REPORT_SRC -> existing_report.md"
else
  echo "(note: no top-level *.md report in $TEST_FOLDER — the report will be formatted fresh)"
fi

make_config() { # src dst gpu(true|false) output_name
  local src="$1" dst="$2" gpu="$3" outname="$4"
  run_env "$GPU_ENV" python -c '
import sys, yaml
src, dst, gpu, out = sys.argv[1:5]
cfg = yaml.safe_load(open(src)) or {}
p = dict(cfg.get("parameters") or {})
p["gpu"] = (gpu == "true")
if gpu == "true":
    p.setdefault("gpu_device", 0)
    p["gpu_verbose"] = True
cfg["parameters"] = p
cfg["output"] = out
with open(dst, "w") as fh:
    fh.write(yaml.safe_dump(cfg, sort_keys=False))
print("  wrote %s  (gpu: %s, output: %s)" % (dst, p["gpu"], out))
' "$src" "$dst" "$gpu" "$outname"
}

info "writing derived configs into $RUN_DIR"
# The runs execute with CWD=$TEST_FOLDER (where the config's data files live),
# so the output dir must be qualified with the run-dir name to stay inside it.
make_config "$CONFIG_SRC" "$RUN_DIR/config.gpu.yaml" true  "gpu_test_${STAMP}/gpu"
make_config "$CONFIG_SRC" "$RUN_DIR/config.cpu.yaml" false "gpu_test_${STAMP}/cpu"

# --- 3. GPU run ----------------------------------------------------------------
info "GPU RUN (env $GPU_ENV, sampler: emcee = PyDE + MCMC), full GPU"
# CWD = $TEST_FOLDER: the config's data-file paths must resolve exactly as they
# did for the user's previous runs in this folder. All test artifacts of this
# script are still written explicitly into $RUN_DIR.
cd "$TEST_FOLDER"
nvidia-smi > "$RUN_DIR/nvidia-smi_before.txt" 2>&1 || true
nvidia-smi dmon -s um -o D > "$RUN_DIR/gpu_usage.log" 2>&1 &
SMON_PID=$!

T_GPU_0="$(date +%s.%N)"
RC_GPU=0
run_env "$GPU_ENV" python "$REPO_DIR/PyORBIT_Run.py" emcee "$RUN_DIR/config.gpu.yaml" > "$RUN_DIR/gpu_run.log" 2>&1 \
  || RC_GPU=$?
T_GPU_1="$(date +%s.%N)"

kill "$SMON_PID" 2>/dev/null || true
wait "$SMON_PID" 2>/dev/null || true
SMON_PID=""
nvidia-smi > "$RUN_DIR/nvidia-smi_after.txt" 2>&1 || true

WALL_GPU="$(awk -v a="$T_GPU_0" -v b="$T_GPU_1" 'BEGIN{printf "%.1f", b-a}')"
PYDE_GPU="$(awk '/PyDE completed, it took/ {s=$0; sub(/^.*it took */, "", s); sub(/ *seconds.*/, "", s); print s; exit}' "$RUN_DIR/gpu_run.log" || true)"
echo "GPU run: rc=$RC_GPU  wall=${WALL_GPU}s  pyde=${PYDE_GPU:-n/a}s   (log: $RUN_DIR/gpu_run.log)"
if [ "$RC_GPU" -ne 0 ]; then
  echo "WARNING: GPU run FAILED — see gpu_run.log (a GPUSupportedModelError is the expected, loud failure)."
  echo "         Continuing with the CPU run so the report can still be built."
fi

# --- 4. CPU run ----------------------------------------------------------------
info "CPU RUN (env $GPU_ENV, gpu:false) — identical stack, no GPU"
T_CPU_0="$(date +%s.%N)"
RC_CPU=0
run_env "$GPU_ENV" python "$REPO_DIR/PyORBIT_Run.py" emcee "$RUN_DIR/config.cpu.yaml" > "$RUN_DIR/cpu_run.log" 2>&1 \
  || RC_CPU=$?
T_CPU_1="$(date +%s.%N)"

WALL_CPU="$(awk -v a="$T_CPU_0" -v b="$T_CPU_1" 'BEGIN{printf "%.1f", b-a}')"
PYDE_CPU="$(awk '/PyDE completed, it took/ {s=$0; sub(/^.*it took */, "", s); sub(/ *seconds.*/, "", s); print s; exit}' "$RUN_DIR/cpu_run.log" || true)"
echo "CPU run: rc=$RC_CPU  wall=${WALL_CPU}s  pyde=${PYDE_CPU:-n/a}s   (log: $RUN_DIR/cpu_run.log)"

# --- 5. best-effort result artifacts (plots/ACF) -------------------------------
echo
echo "==> best-effort: pyorbit_results artifacts (ACF/corner/plots)"
run_env "$GPU_ENV" pyorbit_results emcee "$RUN_DIR/config.gpu.yaml" -all > "$RUN_DIR/pyorbit_results_gpu.log" 2>&1 \
  && echo "    pyorbit_results (gpu): ok" \
  || echo "    pyorbit_results (gpu): FAILED (non-fatal) — see pyorbit_results_gpu.log"
run_env "$GPU_ENV" pyorbit_results emcee "$RUN_DIR/config.cpu.yaml" -all > "$RUN_DIR/pyorbit_results_cpu.log" 2>&1 \
  && echo "    pyorbit_results (cpu): ok" \
  || echo "    pyorbit_results (cpu): FAILED (non-fatal) — see pyorbit_results_cpu.log"

# --- 6. collect: timings.json + report_raw.md ----------------------------------
cat > "$RUN_DIR/_collect.py" <<'PYEOF'
import datetime
import json
import os
import re
import sys

run_dir = sys.argv[1]

def read_log(path):
    try:
        with open(path, "r", errors="replace") as fh:
            return fh.read()
    except OSError:
        return ""

def pyde_seconds(text):
    m = re.search(r"PyDE completed, it took\s+([\d.]+) seconds", text)
    return float(m.group(1)) if m else None

# posterior table rows (chain case):  name  idx  median  +1sigma  -1sigma  (15-84 p)  [b0, b1]
ROW = re.compile(
    r"^\s*(\S+)\s+\d+\s+"
    r"(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s+"
    r"(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s+"
    r"(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s+\(15-84 p\)",
    re.MULTILINE)

def last_posterior(text):
    """{name: median} of the LAST posterior table printed in the log."""
    out = {}
    for m in ROW.finditer(text):
        out[m.group(1)] = float(m.group(2))
    return out

def usage_stats(path):
    peak, mem_sum, n = 0, 0, 0
    try:
        with open(path) as fh:
            for line in fh:
                parts = line.split()
                if not parts or not parts[0].isdigit():
                    continue
                # nvidia-smi dmon columns: id  bus  sm  mem  enc  dec  pm
                if len(parts) >= 4:
                    try:
                        sm, mem = int(parts[2]), int(parts[3])
                    except ValueError:
                        continue
                    peak = max(peak, sm)
                    mem_sum += mem
                    n += 1
    except OSError:
        pass
    return {"gpu_util_peak_pct": peak,
            "gpu_mem_mean_mib": round(mem_sum / n, 1) if n else None,
            "samples": n}

def nvidia_snap(path):
    out = []
    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if line:
                    out.append(line)
    except OSError:
        pass
    return out

def env_int(name, default):
    try:
        return int(os.environ.get(name, str(default)))
    except ValueError:
        return default

def env_float(name, default):
    try:
        return float(os.environ.get(name, str(default)))
    except ValueError:
        return default

g = read_log(os.path.join(run_dir, "gpu_run.log"))
c = read_log(os.path.join(run_dir, "cpu_run.log"))

gpu_med = last_posterior(g)
cpu_med = last_posterior(c)
parity = {}
for name in sorted(set(gpu_med) & set(cpu_med)):
    a, b = gpu_med[name], cpu_med[name]
    parity[name] = {"gpu": a, "cpu": b,
                    "abs_diff": abs(a - b),
                    "rel_diff": (abs(a - b) / abs(a)) if a else None}
max_abs = max((p["abs_diff"] for p in parity.values()), default=None)

rc_gpu, rc_cpu = env_int("RC_GPU", 1), env_int("RC_CPU", 1)
timings = {
    "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
    "run_dir": run_dir,
    "config_src": os.environ.get("CONFIG_SRC", ""),
    "gpu": {
        "status": "ok" if rc_gpu == 0 else "FAILED",
        "rc": rc_gpu,
        "wall_s": env_float("WALL_GPU", 0.0),
        "pyde_s": pyde_seconds(g),
        "nvidia_before": nvidia_snap(os.path.join(run_dir, "nvidia-smi_before.txt")),
        "nvidia_after": nvidia_snap(os.path.join(run_dir, "nvidia-smi_after.txt")),
    },
    "cpu": {
        "status": "ok" if rc_cpu == 0 else "FAILED",
        "rc": rc_cpu,
        "wall_s": env_float("WALL_CPU", 0.0),
        "pyde_s": pyde_seconds(c),
    },
    "gpu_usage": usage_stats(os.path.join(run_dir, "gpu_usage.log")),
    "parity": {
        "n_params_compared": len(parity),
        "max_abs_median_diff": max_abs,
        "per_param": parity,
    },
}
if timings["gpu"]["status"] == "ok" and timings["cpu"]["status"] == "ok" \
        and timings["gpu"]["wall_s"] > 0:
    timings["speedup_cpu_over_gpu"] = round(timings["cpu"]["wall_s"] / timings["gpu"]["wall_s"], 3)

with open(os.path.join(run_dir, "timings.json"), "w") as fh:
    json.dump(timings, fh, indent=2, sort_keys=True)

# ---- report_raw.md -----------------------------------------------------------
def fmt_parity(p):
    if not p:
        return "_no posterior tables found in the logs_"
    lines = ["| parameter | GPU median | CPU median | abs diff |",
             "|---|---|---|---|"]
    for name in sorted(p):
        lines.append("| {n} | {g:.6g} | {c:.6g} | {d:.3g} |".format(
            n=name, g=p[name]["gpu"], c=p[name]["cpu"], d=p[name]["abs_diff"]))
    return "\n".join(lines)

u = timings["gpu_usage"]
md = [
    "# PyORBIT GPU-vs-CPU raw data (Kepler-63b)",
    "",
    "- generated: %s" % timings["generated_utc"],
    "- run dir: `%s`" % run_dir,
    "- config: `%s`" % timings["config_src"],
    "",
    "## Status",
    "",
    "- GPU run: **%s** (rc=%s)" % (timings["gpu"]["status"], timings["gpu"]["rc"]),
    "- CPU run: **%s** (rc=%s)" % (timings["cpu"]["status"], timings["cpu"]["rc"]),
    "",
    "## Timings (wall clock, seconds)",
    "",
    "| phase | GPU | CPU |",
    "|---|---|---|",
    "| PyDE | %s | %s |" % (timings["gpu"]["pyde_s"], timings["cpu"]["pyde_s"]),
    "| total (PyDE + emcee) | %s | %s |" % (timings["gpu"]["wall_s"], timings["cpu"]["wall_s"]),
]
if "speedup_cpu_over_gpu" in timings:
    md.append("| **speedup (CPU wall / GPU wall)** | **%sx** | |" % timings["speedup_cpu_over_gpu"])
md += [
    "",
    "## GPU utilization (nvidia-smi dmon, 1 Hz samples during the GPU run)",
    "",
    "- peak SM utilization: **%s %%**" % u["gpu_util_peak_pct"],
    "- mean GPU memory used: **%s MiB** over %d samples" % (u["gpu_mem_mean_mib"], u["samples"]),
    "",
    "## Posterior parity (medians, 15-84% reported by PyORBIT)",
    "",
    fmt_parity(parity),
    "",
    "max |median diff| (GPU vs CPU): **%s**" % timings["parity"]["max_abs_median_diff"],
    "",
]
with open(os.path.join(run_dir, "report_raw.md"), "w") as fh:
    fh.write("\n".join(md))

print("wrote timings.json + report_raw.md")
PYEOF

export RC_GPU RC_CPU WALL_GPU WALL_CPU CONFIG_SRC
run_env "$GPU_ENV" python "$RUN_DIR/_collect.py" "$RUN_DIR"

# --- done ----------------------------------------------------------------------
echo
info "DONE"
echo "  run dir     : $RUN_DIR"
echo "  GPU run     : rc=$RC_GPU  wall=${WALL_GPU}s"
echo "  CPU run     : rc=$RC_CPU  wall=${WALL_CPU}s"
echo "  logs        : gpu_run.log  cpu_run.log  gpu_usage.log"
echo "  report data : timings.json  report_raw.md  nvidia-smi_before/after.txt"
echo
echo "  Next: restart ollama, come back to the ollama conversation"
echo "        and ask to rebuild the performance report from $RUN_DIR"
