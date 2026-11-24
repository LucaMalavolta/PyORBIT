#!/usr/bin/env bash
# =============================================================================
# check_gpu_usage.sh — live monitor for the RTX PRO 4500 Blackwell (32 GiB).
#
#   Run this in a SECOND terminal while an A/B test is going (e.g.
#   run_nested_test_phot_WIDE_Trend.sh) to WATCH the GPU actually work —
#   SM / memory utilization, VRAM in use, power draw, clocks, temperature,
#   and the per-process VRAM breakdown (so you can see the dynesty 8-worker
#   pool AND/OR ollama holding VRAM at a glance).
#
#   It ONLY READS from nvidia-smi. It launches nothing and does no GPU work —
#   so it is always safe to run alongside the benchmark (it adds no VRAM).
#
# USAGE
#   bash check_gpu_usage.sh                 # live, refresh every 1 s  (Ctrl-C stops)
#   bash check_gpu_usage.sh -n 2            # live, refresh every 2 s
#   bash check_gpu_usage.sh --once          # ONE snapshot, then exit (pipe/tee-safe)
#   bash check_gpu_usage.sh --log gpu.csv   # also append one CSV row per sample
#                                           # (combine with live mode to record)
#   bash check_gpu_usage.sh -h              # this help
#
#   Examples:
#     # watch it live while the photometric dynesty leg runs:
#     bash check_gpu_usage.sh -n 1 --log ~/pyorbit-tests/Kepler-63b/gpu_live.csv
#     # a single snapshot into a file:
#     bash check_gpu_usage.sh --once | tee gpu_snapshot.txt
#
#   CSV columns (one row per sample):
#     timestamp_utc, sm_util_pct, mem_util_pct, mem_used_mib, mem_free_mib,
#     mem_total_mib, power_w, temp_c, sm_mhz
#
#   WHY THE PER-PROCESS LIST: the GPU leg of the A/B test spawns 8 JAX worker
#   processes (each opens its own CUDA context + reserves gpu_vram_fraction of
#   VRAM). If the monitor shows only `llama-server` / `ollama` holding tens of
#   GiB and NO pyorbit python processes, the pool never started (or ollama is
#   still hogging VRAM) — that is the "GPU sits at ~1% power" symptom.
# =============================================================================
set -uo pipefail

# --- defaults / argument parsing --------------------------------------------
INTERVAL="1"
MODE="live"
LOG=""

usage() { sed -n '2,30p' "$0"; }

while [ $# -gt 0 ]; do
  case "$1" in
    -h|--help)          usage; exit 0 ;;
    --once|-1)          MODE="once" ;;
    -n|--interval)      INTERVAL="${2:-}"; shift ;;
    -n=*|--interval=*)  INTERVAL="${1#*=}" ;;
    --log)              LOG="${2:-}"; shift ;;
    --log=*)            LOG="${1#*=}" ;;
    --)                 shift; break ;;
    -*)                 echo "unknown option: $1 (try --help)" >&2; exit 2 ;;
    *)                  echo "unexpected arg: $1 (try --help)"  >&2; exit 2 ;;
  esac
  shift
done

command -v nvidia-smi >/dev/null 2>&1 || { echo "ERROR: nvidia-smi not found (NVIDIA driver missing?)" >&2; exit 1; }
case "$INTERVAL" in ''|*[!0-9.]*) echo "ERROR: interval must be a number, got '$INTERVAL'" >&2; exit 2 ;; esac

# --- small helpers -----------------------------------------------------------
# trim leading/trailing whitespace (nvidia-smi CSV pads fields with a space)
trim() { local s="$1"; s="${s#"${s%%[![:space:]]*}"}"; s="${s%"${s##*[![:space:]]}"}"; printf '%s' "$s"; }

# resolve a pid to a short human command (falls back to nvidia-smi's name)
short_cmd() {
  local pid="$1" fallback="$2" c
  if [ -r "/proc/$pid/cmdline" ]; then
    c=$(tr '\0' ' ' < "/proc/$pid/cmdline" 2>/dev/null | sed -E 's/[[:space:]]+/ /g; s/^ +//; s/ +$//')
    if [ -n "$c" ]; then
      [ ${#c} -gt 78 ] && c="${c:0:75}..."
      printf '%s' "$c"; return 0
    fi
  fi
  printf '%s' "$fallback"
}

# --- gather device metrics into globals --------------------------------------
# One combined query (noheader, nounits) keeps this a single NVML call for the
# numbers; a second for the string labels (name/driver/pcie).
gather() {
  local line
  if ! line=$(nvidia-smi --query-gpu=utilization.gpu,utilization.memory,memory.used,memory.free,memory.total,power.draw,power.limit,temperature.gpu,clocks.current.sm,clocks.max.sm,clocks.current.memory,clocks.max.memory --format=csv,noheader,nounits 2>/dev/null); then
    GATHER_OK=0; return 0
  fi
  GATHER_OK=1
  local v; IFS=', ' read -r -a v <<< "$line"
  SMU=${v[0]:-?};  MEMU=${v[1]:-?}
  USED=${v[2]:-?}; FREE=${v[3]:-?}; TOTAL=${v[4]:-?}
  PDRAW=${v[5]:-?}; PLIMIT=${v[6]:-?}; TEMP=${v[7]:-?}
  SMCUR=${v[8]:-?}; SMX=${v[9]:-?}; MEMCUR=${v[10]:-?}; MEMX=${v[11]:-?}
  local labels l
  labels=$(nvidia-smi --query-gpu=name,driver_version,pcie.link.gen.current,pcie.link.width.current --format=csv,noheader 2>/dev/null) || labels="?, ?, ?, ?"
  IFS=',' read -r -a l <<< "$labels"
  NAME=$(trim "${l[0]:-?}"); DRIVER=$(trim "${l[1]:-?}")
  PCIEGEN=$(trim "${l[2]:-?}"); PCIEX=$(trim "${l[3]:-?}")
}

# --- gather per-process VRAM into PROCS_LINES (and detect ollama) ------------
PROCS_LINES=""
PROCS_EMPTY=1
OLLAMA_MEM=0
gather_procs() {
  PROCS_LINES=""; PROCS_EMPTY=1; OLLAMA_MEM=0
  local out
  out=$(nvidia-smi --query-compute-apps=pid,used_memory,process_name --format=csv,noheader 2>/dev/null) || true
  [ -n "$out" ] || return 0
  local line pid mem name cmd
  while IFS=',' read -r pid mem name; do
    pid=$(trim "$pid"); [ -n "$pid" ] || continue
    mem=$(trim "$mem"); mem=${mem%% *}            # strip trailing " MiB"
    name=$(trim "$name")
    cmd=$(short_cmd "$pid" "$name")
    PROCS_LINES+=$(printf '  %9s  %9s MiB  %s\n' "$pid" "$mem" "$cmd")$'\n'
    PROCS_EMPTY=0
    case "$(printf '%s' "$cmd")" in
      *ollama*|*llama-server*) OLLAMA_MEM=$(( OLLAMA_MEM + ${mem:-0} )) ;;
    esac
  done <<< "$out"
}

# --- render one snapshot ------------------------------------------------------
render() {
  local pct free_gib
  pct=$(awk -v u="${USED:-0}"  -v t="${TOTAL:-0}" 'BEGIN{printf "%.1f", (t>0)?100*u/t:0}')
  free_gib=$(awk -v f="${FREE:-0}"  'BEGIN{printf "%.1f", f/1024}')
  echo "============================================================================="
  echo "  GPU check   $(date -u +%Y-%m-%dT%H:%M:%SZ)   $(hostname)   (sample #$SAMPLE_N)"
  echo "  ${NAME}    driver ${DRIVER}    PCIe Gen${PCIEGEN} x${PCIEX}"
  echo "============================================================================="
  printf '  SM util   : %5s %%      mem util : %5s %%\n' "${SMU}" "${MEMU}"
  printf '  VRAM      : %6s / %s MiB used  (%s%%)      free %s MiB  (%s GiB)\n' \
         "${USED}" "${TOTAL}" "${pct}" "${FREE}" "${free_gib}"
  printf '  Power     : %7s / %s W      temp %s C\n' "${PDRAW}" "${PLIMIT}" "${TEMP}"
  printf '  Clocks    : SM %s / %s MHz      Mem %s / %s MHz\n' "${SMCUR}" "${SMX}" "${MEMCUR}" "${MEMX}"
  echo "-------------------------------------------------------------------------------"
  if [ "${PROCS_EMPTY}" = "1" ]; then
    echo "  GPU processes : (none — nothing is using the GPU)"
  else
    echo "  GPU processes (per-process VRAM):"
    printf '   PID         VRAM(MiB)  COMMAND\n'
    printf '%s' "${PROCS_LINES}"
  fi
  if [ "${OLLAMA_MEM}" -gt 0 ]; then
    echo "  ! ollama / llama-server is holding ~${OLLAMA_MEM} MiB of VRAM."
    echo "    Stop it before the 8-worker pool:   pkill -f llama-server   (or: ollama stop)"
  fi
  echo "============================================================================="
}

# --- CSV row (for --log) ------------------------------------------------------
csv_row() {
  local ts
  ts=$(date -u +%Y-%m-%dT%H:%M:%SZ)
  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$ts" "${SMU}" "${MEMU}" "${USED}" "${FREE}" "${TOTAL}" "${PDRAW}" "${TEMP}" "${SMCUR}"
}
csv_ensure_header() {
  [ -n "$LOG" ] || return 0
  [ -s "$LOG" ] || echo "timestamp_utc,sm_util_pct,mem_util_pct,mem_used_mib,mem_free_mib,mem_total_mib,power_w,temp_c,sm_mhz" >> "$LOG"
}

# --- one full sample ----------------------------------------------------------
run_once() {
  gather
  if [ "${GATHER_OK:-0}" != "1" ]; then
    echo "  (nvidia-smi failed this tick — retrying next sample)"
    return 1
  fi
  gather_procs
  render
  csv_ensure_header
  [ -n "$LOG" ] && csv_row >> "$LOG"
  return 0
}

# --- main loop ----------------------------------------------------------------
SAMPLE_N=0
if [ "$MODE" = "once" ]; then
  SAMPLE_N=1
  run_once
  [ -n "$LOG" ] && echo "  (logged to $LOG)"
  exit 0
fi

trap 'echo; echo "(stopped after $SAMPLE_N samples @ ${INTERVAL}s)"; [ -n "$LOG" ] && echo "  log written: $LOG"; exit 0' INT TERM
echo "Monitoring — Ctrl-C to stop.  (interval ${INTERVAL}s${LOG:+, logging to $LOG})"
while :; do
  [ -t 1 ] && printf '\033[2J\033[H'   # clear screen only on a real terminal
  SAMPLE_N=$((SAMPLE_N + 1))
  run_once
  sleep "$INTERVAL"
done
