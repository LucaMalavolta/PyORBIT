#!/usr/bin/env bash
# =============================================================================
# gpu_clean_run.sh — one-shot, unattended driver for the PyORBIT GPU-vs-CPU
# honest precision + performance statistics.
#
# Design rules (fixed for this deliverable):
#   * `set -u` only — deliberately NOT `set -e`: every leg may fail without
#     aborting the run; every leg's outcome is recorded in $OUT/legs.tsv and
#     its full output in a $OUT/Px_*.log file.
#   * Never kills/stops foreign processes (e.g. ollama). Per-leg timeouts use
#     coreutils `timeout`, which runs the child in its OWN process group and
#     signals *that group* on expiry (escalating to KILL) — i.e. only our own
#     leg tree can ever be killed.
#   * All logs/artifacts live under OUT (default
#     /tmp/opencode/pyorbit_clean_run; override with OUT=... bash ...).
#   * Stages:
#       P0  environment / GPU-contention audit        -> P0_env.log + flags.env
#       P1  precision (parity suite + b1/b3/b4/b5/b6) -> P1_*.log
#       P2  performance surfaces (reps>=11, N ceilings) -> P2_perf.log + P2_logs/
#       P3  emcee A/B (standalone leg)                -> P3_emcee.log
#       P4  Kepler-63b A/B (author's script, verbatim)-> P4_kepler63b.log + P4_artifacts/
#       P5  stats report                              -> HONEST_STATS.md + stats.json
#       P6  compact console summary (starts with `== CLEAN RUN COMPLETE ==`)
#   * MAX_STAGE (default 6, env-overridable) gates P1–P4 for smoke tests;
#     P5/P6 always run. MAX_STAGE=0 => P0 + P5 + P6 only.
#   * Contention policy: ollama reachable OR < 20 GiB free GPU memory =>
#     warning flag(s); the P4 Kepler-63b A/B is skipped under LOW_VRAM (logged),
#     and P2 runs its light N-ceiling list.
#
# Usage:
#   bash development/gpu_clean_run.sh
#   OUT=/path/to/out bash development/gpu_clean_run.sh
#   MAX_STAGE=0 bash development/gpu_clean_run.sh      # plumbing smoke only
# =============================================================================
set -u
export LC_ALL=C

# --- constants / config -------------------------------------------------------
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${OUT:-/tmp/opencode/pyorbit_clean_run}"
MAX_STAGE="${MAX_STAGE:-6}"
case "$MAX_STAGE" in (*[!0-9]*|"") MAX_STAGE=6 ;; esac

PY_ENV="pyorbit"
GPU_ENV="pyorbit-gpu"
TEST_FOLDER="/home/zingales/pyorbit-tests/Kepler-63b"
P4_CONFIG="$TEST_FOLDER/TOI1554_GPU.yaml"
REPS="${REPS:-11}"
WARMUP="${WARMUP:-3}"
case "$REPS"   in (*[!0-9]*|"") REPS=11 ;; esac
case "$WARMUP" in (*[!0-9]*|"") WARMUP=3 ;; esac
[ "$WARMUP" -lt 3 ] && WARMUP=3

# N ceilings (full lists = at/above the pre-existing defaults 10k kepler /
# 5k tinygp). Light mode (LOW_VRAM) reverts to the pre-existing lists.
KEPLER_N="${KEPLER_N:-100,1000,10000,100000}"
TINYGP_N="${TINYGP_N:-1000,5000,10000,20000}"

# timeouts (s): only our own leg process group is ever killed
TMO_LIGHT=900      # P1 legs (b6 gets 1500: it runs two JAX sub-legs)
TMO_B6=1500
TMO_P2=3600        # P2 full performance surfaces
TMO_P3=1800        # P3 emcee A/B
TMO_P4=7200        # P4 both legs inside the author's single call (capped config)
LOW_VRAM_GIB=20    # below this: light legs + P4 skip

UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)
CONTENDED=0
OLLAMA_REACHABLE=0
FREE_MIB=""
FREE_GIB="0"
LOW_VRAM=0
GPU_MODEL="(nvidia-smi unavailable)"
GPU_TOTAL_GIB="?"
P4_SKIPPED=0
P4_SKIP_REASON=""
CONDA=""

# --- boot ---------------------------------------------------------------------
mkdir -p "$OUT" || { echo "FATAL: cannot create OUT dir: $OUT" >&2; exit 2; }

LEGLEDGER="$OUT/legs.tsv"
if [ ! -f "$LEGLEDGER" ]; then
  printf 'name\trc\twall_s\ttimed_out\tstate\tnote\tlog\n' > "$LEGLEDGER"
fi

log_line() { echo "[clean-run] $*"; }

# Record a leg that was deliberately not run (reason in `note`).
record_skip() { # name log note
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$1" "-" "-" "0" "skipped" "$3" "$2" >> "$LEGLEDGER"
  log_line "leg $1: SKIPPED ($3)"
}

# Run one leg: NAME TIMEOUT_S LOGFILE NOTE CMD...
# Output -> $OUT/LOGFILE (stdout+stderr). Result row appended to legs.tsv.
# `timeout` (coreutils) puts the child in its own process group and signals
# that group on expiry (`-k` escalates to KILL after 10 s), so only the leg's
# own tree can be killed — never foreign processes.
run_leg() { # name tmo logfile note cmd...
  local name="$1" tmo="$2" logf="$3" note="${4:-}"
  shift 4
  local t0 t1 rc state to_flag elapsed
  log_line "leg $name: timeout ${tmo}s -> ${OUT}/$logf"
  t0=$(date +%s)
  ( cd "$REPO" && timeout -k 10 "$tmo" "$@" ) > "$OUT/$logf" 2>&1
  rc=$?
  t1=$(date +%s)
  elapsed=$((t1 - t0))
  to_flag=0
  if [ "$rc" -eq 0 ]; then
    state="ok"
  elif [ "$rc" -eq 124 ]; then
    state="timed_out"; to_flag=1
  else
    state="failed"
  fi
  [ -n "$note" ] || note="(rc=$rc)"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$name" "$rc" "$elapsed" "$to_flag" "$state" "$note" "$logf" >> "$LEGLEDGER"
  log_line "leg $name: state=$state rc=$rc wall=${elapsed}s (log ${OUT}/$logf)"
  return 0
}

# flags.env key=VAL (later lines override earlier on re-read)
update_flag() { echo "$1=$2" >> "$OUT/flags.env"; }

# conda: source a profile if present (never modifies envs; run scripts are read-only w.r.t. package installs,
# except the author's P4 script which is idempotent and pre-satisfied here)
for _p in "$HOME/anaconda3/etc/profile.d/conda.sh" \
          "$HOME/miniconda3/etc/profile.d/conda.sh" \
          "$HOME/miniforge3/etc/profile.d/conda.sh" \
          /opt/conda/etc/profile.d/conda.sh; do
  if [ -f "$_p" ]; then
    # shellcheck disable=SC1091
    . "$_p" 2>/dev/null || true
    break
  fi
done
if command -v conda >/dev/null 2>&1; then CONDA="$(command -v conda)"; fi

gpu_free_mib() {
  nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits 2>/dev/null \
    | head -n1 | tr -d '[:space:]'
}

pkg_versions() { # env
  "$CONDA" run --no-capture-output -n "$1" python -c '
import importlib.metadata as md
for name in ("pyorbit","jax","jaxlib","cupy","numpy","scipy","tinygp",
             "batman","spleaf","jaxoplanet"):
    try:
        print("  {0}={1}".format(name, md.version(name)))
    except Exception:
        print("  {0}=not-installed".format(name))
' 2>&1
}

# =============================================================================
log_line "PyORBIT GPU clean-run driver"
log_line "REPO=$REPO"
log_line "OUT=$OUT  MAX_STAGE=$MAX_STAGE"

# --- P0: environment / contamination audit -----------------------------------
log_line "P0: environment audit"
ollama_body="n/a"
ollama_code="n/a"
if command -v curl >/dev/null 2>&1; then
  ollama_body="$(curl -s --max-time 3 http://127.0.0.1:11434/api/tags 2>/dev/null | head -c 300)"
  ollama_code="$(curl -s --max-time 3 -o /dev/null -w '%{http_code}' http://127.0.0.1:11434/api/tags 2>/dev/null || echo err)"
else
  ollama_code="curl-missing"
fi
case "$ollama_code" in
  200) OLLAMA_REACHABLE=1 ;;
  *)   OLLAMA_REACHABLE=0 ;;
esac

if command -v nvidia-smi >/dev/null 2>&1; then
  GPU_MODEL="$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -n1 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
  GPU_TOTAL_MIB="$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits 2>/dev/null | head -n1 | tr -d '[:space:]')"
  GPU_TOTAL_GIB="$(awk -v m="${GPU_TOTAL_MIB:-0}" 'BEGIN{printf "%.1f", m/1024}')"
  FREE_MIB="$(gpu_free_mib)" || FREE_MIB=""
else
  GPU_MODEL="(nvidia-smi unavailable)"
  GPU_TOTAL_GIB="?"
fi
case "$FREE_MIB" in
  ''|*[!0-9]*) FREE_MIB=""; FREE_GIB="0" ;;
  *) FREE_GIB="$(awk -v m="$FREE_MIB" 'BEGIN{printf "%.1f", m/1024}')" ;;
esac
[ "$OLLAMA_REACHABLE" -eq 1 ] && CONTENDED=1
awk -v g="${FREE_GIB:-0}" -v t=$LOW_VRAM_GIB 'BEGIN{exit !(g+0 < t+0)}' && LOW_VRAM=1

if [ "$CONTENDED" -eq 1 ]; then
  log_line "WARN CONTENDED=1: ollama HTTP endpoint reachable (127.0.0.1:11434) — GPU shares VRAM; numbers are not clean-room. Stop ollama for the clean run and re-run."
fi
if [ "$LOW_VRAM" -eq 1 ]; then
  log_line "WARN LOW_VRAM=1: only ${FREE_GIB} GiB of ${GPU_TOTAL_GIB} GiB free (< ${LOW_VRAM_GIB} GiB) — running light legs; P4 Kepler-63b A/B will be SKIPPED."
fi
if [ -z "${KEPLER_N:-}" ]; then KEPLER_N="100,1000,10000,100000"; fi
if [ -z "${TINYGP_N:-}" ]; then TINYGP_N="1000,5000,10000,20000"; fi
if [ "$LOW_VRAM" -eq 1 ]; then
  KEPLER_N="100,1000,10000"
  TINYGP_N="1000,5000"
fi

# P0_env.log (the `--- commands ---` section must remain the LAST section:
# the report script captures it verbatim for provenance)
{
  echo "# P0 environment / contamination audit — $UTC"
  echo "REPO: $REPO"
  echo "OUT: $OUT"
  echo "MAX_STAGE: $MAX_STAGE"
  echo "GPU MODEL: $GPU_MODEL"
  echo "GPU total/free: ${GPU_TOTAL_GIB} GiB total / ${FREE_GIB} GiB free"
  echo "OLLAMA_REACHABLE: $OLLAMA_REACHABLE (CONTENDED=$CONTENDED)"
  echo "LOW_VRAM: $LOW_VRAM (threshold ${LOW_VRAM_GIB} GiB)"
  echo "conda: ${CONDA:-NOT FOUND}"
  echo
  echo "--- ollama probe ---"
  echo "  GET http://127.0.0.1:11434/api/tags  http_code=$ollama_code"
  echo "  body(300B): $(echo "$ollama_body" | head -n3 | tr '\n' ' ')"
  echo
  if command -v nvidia-smi >/dev/null 2>&1; then
    echo "--- nvidia-smi (full) ---"
    nvidia-smi 2>&1
    echo
    echo "--- nvidia-smi query ---"
    nvidia-smi --query-gpu=name,memory.total,memory.free,memory.used,utilization.gpu,utilization.memory --format=csv 2>&1
    echo
    echo "--- GPU process list (foreign processes are NEVER killed) ---"
    nvidia-smi --query-compute-apps=pid,used_memory,process_name --format=csv 2>&1
  else
    echo "--- nvidia-smi: NOT AVAILABLE ---"
  fi
  echo
  if [ -n "$CONDA" ]; then
    echo "--- package versions: env '$PY_ENV' (CPU) ---"
    pkg_versions "$PY_ENV"
    echo
    echo "--- package versions: env '$GPU_ENV' (GPU) ---"
    pkg_versions "$GPU_ENV"
  else
    echo "--- package versions: conda NOT FOUND ---"
  fi
  echo
  echo "--- commands ---"
  echo "[P1] $CONDA run --no-capture-output -n $PY_ENV python development/gpu_cpu_parity.py"
  echo "[P1] $CONDA run --no-capture-output -n $GPU_ENV python development/gpu_cpu_parity.py"
  echo "[P1] $CONDA run --no-capture-output -n $GPU_ENV python development/gpu_parity_b1_transit.py"
  echo "[P1] $CONDA run --no-capture-output -n $GPU_ENV python development/gpu_parity_b3_esp.py"
  echo "[P1] $CONDA run --no-capture-output -n $GPU_ENV python development/gpu_parity_b4_apodized.py"
  echo "[P1] $CONDA run --no-capture-output -n $GPU_ENV python development/gpu_parity_b5_trend.py"
  echo "[P1] $CONDA run --no-capture-output -n $GPU_ENV python development/gpu_parity_b6_tinygp.py"
  echo "[P1] $CONDA run --no-capture-output -n $GPU_ENV python development/emcee_v2_f1.py"
  echo "[P2] $CONDA run --no-capture-output -n $GPU_ENV python development/gpu_perf_surfaces.py --mode all --reps $REPS --warmup $WARMUP --logdir \$OUT/P2_logs --tinygp-n $TINYGP_N --kepler-n $KEPLER_N"
  echo "[P3] $CONDA run --no-capture-output -n $GPU_ENV python development/gpu_perf_surfaces.py --mode emcee --reps $REPS --warmup $WARMUP --logdir \$OUT/P3_logs"
  echo "[P4] PYORBIT_REPO_DIR=$REPO bash development/run_kepler63b_gpu_test.sh $TEST_FOLDER $P4_CONFIG"
  echo "[P5/P6] python3 development/gpu_perf_stats_report.py --out \$OUT [--print-summary]"
} > "$OUT/P0_env.log"

# flags.env (authoritative run metadata for the report aggregator)
{
  echo "UTC=$UTC"
  echo "CONTENDED=$CONTENDED"
  echo "OLLAMA_REACHABLE=$OLLAMA_REACHABLE"
  echo "FREE_MIB=${FREE_MIB:-0}"
  echo "FREE_GIB=${FREE_GIB:-0}"
  echo "LOW_VRAM=$LOW_VRAM"
  echo "GPU_MODEL=$GPU_MODEL"
  echo "GPU_TOTAL_GIB=$GPU_TOTAL_GIB"
  echo "PYORBIT_ENV=$PY_ENV"
  echo "GPU_ENV=$GPU_ENV"
  echo "REPS=$REPS"
  echo "WARMUP=$WARMUP"
  echo "KEPLER_N=$KEPLER_N"
  echo "TINYGP_N=$TINYGP_N"
  echo "P4_CONFIG=$P4_CONFIG"
  echo "P4_SKIPPED=$P4_SKIPPED"
  echo "P4_SKIP_REASON=$P4_SKIP_REASON"
  echo "MAX_STAGE=$MAX_STAGE"
  echo "CONDA=${CONDA:-}"
} > "$OUT/flags.env"
log_line "P0 done (CONTENDED=$CONTENDED LOW_VRAM=$LOW_VRAM free=${FREE_GIB} GiB)"

# --- P1: precision -------------------------------------------------------------
if [ "$MAX_STAGE" -ge 1 ]; then
  log_line "P1: precision suite + probes"
  if [ -n "$CONDA" ]; then
    run_leg parity_pyorbit     $TMO_LIGHT "P1_parity_pyorbit.log"    "parity suite, CPU env" \
      "$CONDA" run --no-capture-output -n "$PY_ENV"  python development/gpu_cpu_parity.py
    run_leg parity_pyorbit_gpu $TMO_LIGHT "P1_parity_pyorbit_gpu.log" "parity suite, GPU env" \
      "$CONDA" run --no-capture-output -n "$GPU_ENV" python development/gpu_cpu_parity.py
    run_leg b1_transit $TMO_LIGHT "P1_b1_transit.log"    "b1 transit" \
      "$CONDA" run --no-capture-output -n "$GPU_ENV" python development/gpu_parity_b1_transit.py
    run_leg b3_esp     $TMO_LIGHT "P1_b3_esp.log"        "b3 esp" \
      "$CONDA" run --no-capture-output -n "$GPU_ENV" python development/gpu_parity_b3_esp.py
    run_leg b4_apodized $TMO_LIGHT "P1_b4_apodized.log"  "b4 apodized" \
      "$CONDA" run --no-capture-output -n "$GPU_ENV" python development/gpu_parity_b4_apodized.py
    run_leg b5_trend   $TMO_LIGHT "P1_b5_trend.log"      "b5 trend" \
      "$CONDA" run --no-capture-output -n "$GPU_ENV" python development/gpu_parity_b5_trend.py
    run_leg b6_tinygp  $TMO_B6    "P1_b6_tinygp.log"     "b6 tinygp (two internal JAX sub-legs)" \
      "$CONDA" run --no-capture-output -n "$GPU_ENV" python development/gpu_parity_b6_tinygp.py
    run_leg f1_blob_parity $TMO_LIGHT "P1_f1_blob_parity.log" "f1 blob-parity + F1-class probe (obs-1)" \
      "$CONDA" run --no-capture-output -n "$GPU_ENV" python development/emcee_v2_f1.py
  else
    for spec in \
      "parity_pyorbit|P1_parity_pyorbit.log" \
      "parity_pyorbit_gpu|P1_parity_pyorbit_gpu.log" \
      "b1_transit|P1_b1_transit.log" \
      "b3_esp|P1_b3_esp.log" \
      "b4_apodized|P1_b4_apodized.log" \
      "b5_trend|P1_b5_trend.log" \
      "b6_tinygp|P1_b6_tinygp.log" \
      "f1_blob_parity|P1_f1_blob_parity.log"; do
      record_skip "${spec%%|*}" "${spec##*|}" "conda not found on this host"
    done
  fi
else
  log_line "P1: SKIPPED (MAX_STAGE=$MAX_STAGE < 1)"
  for spec in \
    "parity_pyorbit|P1_parity_pyorbit.log" \
    "parity_pyorbit_gpu|P1_parity_pyorbit_gpu.log" \
    "b1_transit|P1_b1_transit.log" \
    "b3_esp|P1_b3_esp.log" \
    "b4_apodized|P1_b4_apodized.log" \
    "b5_trend|P1_b5_trend.log" \
    "b6_tinygp|P1_b6_tinygp.log" \
    "f1_blob_parity|P1_f1_blob_parity.log"; do
    record_skip "${spec%%|*}" "${spec##*|}" "MAX_STAGE=$MAX_STAGE < 1"
  done
fi

# --- P2: performance surfaces --------------------------------------------------
if [ "$MAX_STAGE" -ge 2 ]; then
  log_line "P2: performance surfaces (reps=$REPS warmup=$WARMUP tinygp-n=$TINYGP_N kepler-n=$KEPLER_N)"
  if [ -n "$CONDA" ]; then
    rm -rf "$OUT/P2_logs"
    run_leg perf_surfaces $TMO_P2 "P2_perf.log" "gpu_perf_surfaces --mode all" \
      "$CONDA" run --no-capture-output -n "$GPU_ENV" \
        python development/gpu_perf_surfaces.py \
          --mode all --reps "$REPS" --warmup "$WARMUP" \
          --logdir "$OUT/P2_logs" --tinygp-n "$TINYGP_N" --kepler-n "$KEPLER_N"
  else
    record_skip perf_surfaces P2_perf.log "conda not found on this host"
  fi
else
  log_line "P2: SKIPPED (MAX_STAGE=$MAX_STAGE < 2)"
  record_skip perf_surfaces P2_perf.log "MAX_STAGE=$MAX_STAGE < 2"
fi

# --- P3: emcee A/B ---------------------------------------------------------------
if [ "$MAX_STAGE" -ge 3 ]; then
  log_line "P3: emcee A/B (standalone leg)"
  if [ -n "$CONDA" ]; then
    run_leg emcee_ab $TMO_P3 "P3_emcee.log" "gpu_perf_surfaces --mode emcee" \
      "$CONDA" run --no-capture-output -n "$GPU_ENV" \
        python development/gpu_perf_surfaces.py \
          --mode emcee --reps "$REPS" --warmup "$WARMUP" --logdir "$OUT/P3_logs"
  else
    record_skip emcee_ab P3_emcee.log "conda not found on this host"
  fi
else
  log_line "P3: SKIPPED (MAX_STAGE=$MAX_STAGE < 3)"
  record_skip emcee_ab P3_emcee.log "MAX_STAGE=$MAX_STAGE < 3"
fi

# --- P4: Kepler-63b A/B (author's script, verbatim) ------------------------------
if [ "$MAX_STAGE" -ge 4 ]; then
  log_line "P4: Kepler-63b A/B (author script, verbatim)"
  if [ "$LOW_VRAM" -eq 1 ]; then
    P4_SKIPPED=1
    P4_SKIP_REASON="LOW_VRAM: free ${FREE_GIB} GiB < ${LOW_VRAM_GIB} GiB — the P4 A/B would thrash the shared card; light legs only (see author script MIN_FREE_MIB=8192 fail-fast)"
    record_skip kepler63b_ab P4_kepler63b.log "$P4_SKIP_REASON"
    update_flag P4_SKIPPED "$P4_SKIPPED"
    update_flag P4_SKIP_REASON "$P4_SKIP_REASON"
  else
    # pre-leg VRAM guard (re-check: P0 can be minutes old). Never touch
    # foreign processes — the author script's own fail-fast is preserved.
    CUR_MIB="$(gpu_free_mib)" || CUR_MIB=""
    case "$CUR_MIB" in
      *[!0-9]*|"") ;;
      *)
        if [ "$CUR_MIB" -lt 8192 ]; then
          P4_SKIPPED=1
          P4_SKIP_REASON="pre-P4 guard: only ${CUR_MIB} MiB free (< 8192 MiB author threshold)"
          record_skip kepler63b_ab P4_kepler63b.log "$P4_SKIP_REASON"
          update_flag P4_SKIPPED "$P4_SKIPPED"
          update_flag P4_SKIP_REASON "$P4_SKIP_REASON"
        fi
        ;;
    esac
  fi

  if [ "$P4_SKIPPED" -eq 0 ]; then
    if [ -n "$CONDA" ] && [ -f "$P4_CONFIG" ] && [ -d "$TEST_FOLDER" ]; then
      run_leg kepler63b_ab $TMO_P4 "P4_kepler63b.log" \
        "run_kepler63b_gpu_test.sh (both legs inside)" \
        env PYORBIT_REPO_DIR="$REPO" bash development/run_kepler63b_gpu_test.sh \
          "$TEST_FOLDER" "$P4_CONFIG"
    else
      P4_SKIPPED=1
      P4_SKIP_REASON="missing conda or P4 inputs ($P4_CONFIG / $TEST_FOLDER)"
      record_skip kepler63b_ab P4_kepler63b.log "$P4_SKIP_REASON"
      update_flag P4_SKIPPED "$P4_SKIPPED"
      update_flag P4_SKIP_REASON "$P4_SKIP_REASON"
    fi
  fi

  # collect the author-run artifacts if a run dir was created
  RUN_DIR_LINE="$(sed -n 's/^[[:space:]]*run dir[[:space:]]*:[[:space:]]*//p' "$OUT/P4_kepler63b.log" 2>/dev/null | tail -n1)"
  if [ -n "$RUN_DIR_LINE" ] && [ -d "$RUN_DIR_LINE" ]; then
    mkdir -p "$OUT/P4_artifacts"
    for f in timings.json report_raw.md nvidia-smi_before.txt nvidia-smi_after.txt gpu_usage.log; do
      [ -f "$RUN_DIR_LINE/$f" ] && cp -f "$RUN_DIR_LINE/$f" "$OUT/P4_artifacts/" 2>/dev/null || true
    done
    for f in gpu_run.log cpu_run.log; do
      if [ -f "$RUN_DIR_LINE/$f" ]; then
        { echo "## tail of $f ($RUN_DIR_LINE/$f)"; tail -n 60 "$RUN_DIR_LINE/$f"; } \
          > "$OUT/P4_artifacts/${f%.log}.tail.log" 2>/dev/null || true
      fi
    done
    log_line "P4 artifacts collected from run dir: $RUN_DIR_LINE"
  fi
else
  log_line "P4: SKIPPED (MAX_STAGE=$MAX_STAGE < 4)"
  record_skip kepler63b_ab P4_kepler63b.log "MAX_STAGE=$MAX_STAGE < 4"
fi

# --- P5 + P6: stats report + compact console summary ----------------------------
log_line "P5: stats report (HONEST_STATS.md + stats.json)"
( cd "$REPO" && exec python3 development/gpu_perf_stats_report.py --out "$OUT" ) \
  > "$OUT/P5_report.log" 2>&1
RC5=$?
if [ "$RC5" -ne 0 ]; then
  log_line "WARN: P5 report rc=$RC5 (see P5_report.log); continuing to P6"
fi

log_line "P6: compact console summary"
if ! ( cd "$REPO" && python3 development/gpu_perf_stats_report.py --out "$OUT" --print-summary 2>&1 | tee "$OUT/P6_summary.log" ); then
  log_line "ERROR: P6 summary failed (see P6_summary.log)"
  echo "== CLEAN RUN INCOMPLETE (P6 failed) =="
  echo "OUT=$OUT"
  exit 1
fi

log_line "done. all artifacts under: $OUT"
exit 0
