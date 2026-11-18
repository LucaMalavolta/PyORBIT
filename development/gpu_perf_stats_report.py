#!/usr/bin/env python3
"""PyORBIT-GPU clean-run stats aggregator.

Parses the artifacts of ``development/gpu_clean_run.sh`` under ``--out DIR``
and writes two files into that directory:

  * ``HONEST_STATS.md`` -- human-readable report, caveats first, every figure
    traceable to a specific log in OUT (nothing invented; missing data is
    reported as not-run with the reason and log path).
  * ``stats.json``      -- machine-readable: per-rep wall-time arrays,
    per-leg stats (median/mean/std/min/max/IQR), parity blocks, contention
    snapshots, leg ledger, and file hashes for traceability.

Inputs expected in OUT (produced by gpu_clean_run.sh; all optional -- a
missing input degrades to a "not run" row rather than a crash):
  flags.env, legs.tsv, P0_env.log,
  P1_parity_pyorbit.log, P1_parity_pyorbit_gpu.log,
  P1_b1_transit.log, P1_b3_esp.log, P1_b4_apodized.log, P1_b5_trend.log,
  P1_b6_tinygp.log,
  P2_perf.log, P2_logs/{bench.log,results.json,contention.log,emcee_T2.log},
  P3_emcee.log, P3_logs/bench.log,
  P4_kepler63b.log, P4_artifacts/{timings.json,report_raw.md},
  P5_stats.log.

Usage
-----
  python3 development/gpu_perf_stats_report.py --out /tmp/opencode/pyorbit_clean_run
  python3 development/gpu_perf_stats_report.py --out DIR --print-summary

Stdlib only.
"""

import argparse
import datetime
import hashlib
import json
import math
import os
import re
import statistics
import sys

# --------------------------------------------------------------------------
# contest-session reference values (NOT recomputed; used only to state whether
# the clean-env numbers match) -- taken 2026-11-17 from the verified logs in
# /tmp/opencode/pyorbit_gpu_verify/ (GPU was contended then: ~26.7 GiB of the
# 32.6 GiB card held by ollama).
# --------------------------------------------------------------------------
CONTESTED_SRC = "/tmp/opencode/pyorbit_gpu_verify/*.log (2026-11-17)"
CONTESTED_REFERENCE = {
    "parity-suite (pyorbit CPU env)": {
        "result": "PASS",
        "summary": "5 passed, 1 skipped",
        "worst_delta": None,
        "band": None,
    },
    "parity-suite (pyorbit-gpu)": {
        "result": "PASS",
        "summary": "8 passed, 0 skipped",
        "worst_delta": None,
        "band": None,
    },
    "b1 transit": {"result": "PASS", "worst_delta": 7.794456e-09, "band": 1e-03},
    "b3 esp": {"result": "PASS", "worst_delta": 1.742301e-07, "band": 1e-06},
    "b4 apodized": {"result": "PASS", "worst_delta": 7.897982e-13, "band": 1e-10},
    "b5 trend": {"result": "PASS", "worst_delta": 0.0, "band": 1e-10},
    "b6 tinygp": {"result": "PASS", "worst_delta": 1.935405e-09, "band": 1e-08},
}

# (file in OUT, display name)
PRECISION_SOURCES = (
    ("P1_parity_pyorbit.log", "parity-suite (pyorbit CPU env)"),
    ("P1_parity_pyorbit_gpu.log", "parity-suite (pyorbit-gpu)"),
    ("P1_b1_transit.log", "b1 transit"),
    ("P1_b3_esp.log", "b3 esp"),
    ("P1_b4_apodized.log", "b4 apodized"),
    ("P1_b5_trend.log", "b5 trend"),
    ("P1_b6_tinygp.log", "b6 tinygp"),
)

RESULT_RE = re.compile(r"^RESULT:\s*(PASS|SUSPICIOUS|BUG|FAIL)\b(.*)$")
NUM = r"([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)"
DELTA_PATTERNS = (
    ("max|diff|", re.compile(r"max\|diff\|=\s*" + NUM)),
    ("max rel", re.compile(r"max rel=\s*" + NUM)),
    ("|delta| log_probability", re.compile(r"\|delta\| log_probability=\s*" + NUM)),
)
BAND_RE = re.compile(r"(?:<=|<)\s*" + NUM)

EMCEE_LEG_RE = re.compile(
    r"RESULT emcee leg=(\w+) n=(\d+) ndim=(\d+) nwalkers=(\d+) "
    r"steps_run=(\d+) median_steps_per_s=" + NUM + r" walls_s=([\d. ]+)"
)
EMCEE_SP_RE = re.compile(r"RESULT emcee speedup_gpu_vs_cpu_steps_per_s=" + NUM)
K63_LEG_RE_TPL = r"{leg} run: rc=(\d+)  wall=([0-9.]+)s  pyde=([0-9.]+|n/a)s"


def _read(path):
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            return fh.read()
    except OSError:
        return None


def _hash(path):
    try:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()[:16]
    except OSError:
        return None


def _n(s):
    """Parse a numeric field; None when absent/invalid."""
    try:
        v = float(s)
        if math.isnan(v) or math.isinf(v):
            return None
        return v
    except (TypeError, ValueError):
        return None


def _f(x, spec="{0:.4g}"):
    if x is None:
        return "-"
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if v == 0:
        return "0"
    if abs(v) < 1e-6 or abs(v) >= 1e5:
        return "{0:.3e}".format(v)
    return spec.format(v)


def _mdcell(s):
    """Escape a value for a markdown table cell (pipes break the table)."""
    return str(s).replace("|", "\\|")


# --------------------------------------------------------------- precision --
def parse_precision_log(path, name):
    """Extract verdict/worst-delta/band from a parity suite or probe log."""
    text = _read(path)
    row = {
        "surface": name,
        "log": os.path.basename(path),
        "result": None,
        "pass": None,
        "summary": None,
        "worst_delta": None,
        "band": None,
        "pass_count": None,
        "skip_count": None,
        "fail_count": None,
        "reference": CONTESTED_REFERENCE.get(name),
        "matches_reference": None,
        "missing": text is None,
    }
    if text is None:
        return row

    last_result = None
    for line in text.splitlines():
        m = RESULT_RE.match(line.strip())
        if m:
            last_result = (m.group(1), m.group(2).strip())
    if last_result:
        result, rest = last_result
        row["result"] = result
        row["pass"] = result == "PASS"
        row["summary"] = rest
        mcount = re.search(r"(\d+) passed, (\d+) skipped", rest)
        if mcount:
            row["pass_count"] = int(mcount.group(1))
            row["skip_count"] = int(mcount.group(2))
        nfail = len(re.findall(r"\[FAIL\]", text))
        if nfail:
            row["fail_count"] = nfail
        for _label, pat in DELTA_PATTERNS:
            m = pat.search(rest)
            if m:
                row["worst_delta"] = _n(m.group(1))
                break
        mb = BAND_RE.search(rest)
        if mb:
            row["band"] = _n(mb.group(1))
    else:
        # no RESULT line: a crashed/timed-out leg -- say so explicitly
        row["result"] = "NOT_RUN"
        row["pass"] = False

    ref = row["reference"]
    if ref is None:
        row["matches_reference"] = None
    elif row["missing"] or row["result"] is None:
        row["matches_reference"] = None
    elif "summary" in ref and ref["summary"]:
        m = re.search(r"(\d+) passed, (\d+) skipped", row["summary"] or "")
        row["matches_reference"] = bool(m and m.group(0) == ref["summary"])
    elif ref.get("worst_delta") is not None and row["worst_delta"] is not None:
        old, new = ref["worst_delta"], row["worst_delta"]
        tol = max(1e-15, abs(old) * 1e-6)
        row["matches_reference"] = abs(new - old) <= tol
    return row


# --------------------------------------------------------------- performance -
def parse_label(label):
    """'kepler N=100  [CPU]' / 'tinygp gp N=1000' / 'mc cpu n=600' -> parts."""
    l = (label or "").strip()
    jitref = l.endswith("[GPU-jit-reference]")
    body = l
    mbr = re.search(r"\[\s*(CPU|GPU)\s*\]", body)
    leg = mbr.group(1).lower() if mbr else None
    if jitref:
        body = body[: -len("[GPU-jit-reference]")]
        leg = "gpu-jitref"
    toks = body.split()
    surface = toks[0] if toks else "?"
    mm = re.search(r"[nN]\s*=\s*([0-9]+)", body)
    n = int(mm.group(1)) if mm else None
    if leg is None:
        for t in toks[1:3]:
            if t.lower() in ("cpu", "gpu"):
                leg = t.lower()
                break
    return surface, n, leg, jitref


def _stats_block(rec):
    times = rec.get("times_s") or []
    block = {
        "n_reps": rec.get("reps"),
        "warmup": rec.get("warmup"),
        "cold_ms": (rec.get("cold_s") or 0.0) * 1000.0,
        "median_ms": (rec.get("med_s") or 0.0) * 1000.0,
        "mean_ms": None,
        "std_ms": None,
        "min_ms": (rec.get("min_s") or 0.0) * 1000.0,
        "max_ms": (rec.get("max_s") or 0.0) * 1000.0,
        "iqr_ms": None,
        "parity_guard": None,
        "times_ms": [t * 1000.0 for t in times],
    }
    if len(times) >= 2:
        block["mean_ms"] = statistics.fmean(times) * 1000.0
        block["std_ms"] = statistics.pstdev(times) * 1000.0
        if len(times) >= 4:
            q = statistics.quantiles(times, n=4)
            block["iqr_ms"] = (q[2] - q[0]) * 1000.0
    elif times:
        block["mean_ms"] = times[0] * 1000.0
    note = rec.get("parity_guard")
    if note:
        # notes like "parity max|diff|=1.2e-14",
        # "max|diff| (perturbed geometry) = 1.2e-3",
        # "|dloglike|=1.0e-08 (B3 parity PASS<=1e-6)"
        mm = re.search(r"(?:max\|diff\||max rel|dloglike)[^=]*=\s*" + NUM, note)
        if mm:
            block["parity_guard"] = _n(mm.group(1))
    return block


def load_perf_records(out_dir):
    """Bench records (last session) from P2_logs/bench.log, fallback
    results.json; returns (records_by_label, source_relpath)."""
    bench_path = os.path.join(out_dir, "P2_logs", "bench.log")
    text = _read(bench_path)
    by_label = {}
    src = None
    if text is not None:
        # last session only (bench.log is append-only across runs)
        idx = text.rfind("### === session")
        chunk = text[idx:] if idx >= 0 else text
        for line in chunk.splitlines():
            line = line.strip()
            if not line.startswith("{"):
                continue
            payload, note = (line.split("    #", 1) + [""])[:2]
            try:
                rec = json.loads(payload)
            except ValueError:
                continue
            lab = rec.get("label")
            if not lab:
                continue
            rec["parity_guard"] = note.strip()
            by_label[lab] = rec
        src = "P2_logs/bench.log"
    if not by_label:
        res_path = os.path.join(out_dir, "P2_logs", "results.json")
        text = _read(res_path)
        if text is not None:
            try:
                for rec in json.loads(text):
                    lab = rec.get("label")
                    if lab:
                        by_label[lab] = rec
            except ValueError:
                pass
            src = "P2_logs/results.json"
    return by_label, src


def build_performance(out_dir):
    records, src = load_perf_records(out_dir)
    groups = {}
    for lab, rec in records.items():
        surface, n, leg, _jitref = parse_label(lab)
        if leg is None:
            leg = "unknown"
        groups.setdefault((surface, n), {}).setdefault(leg, rec)

    rows = []
    for (surface, n), legs in sorted(
        groups.items(),
        key=lambda kv: (kv[0][0], kv[0][1] if kv[0][1] is not None else -1),
    ):
        leg_stats = {k: _stats_block(v) for k, v in legs.items()}
        speedup = None
        parity_note = None
        cpu = legs.get("cpu")
        gpu = legs.get("gpu")
        if cpu is not None and gpu is not None:
            cg = gpu.get("med_s") or 0.0
            cc = cpu.get("med_s") or 0.0
            if cg > 0:
                speedup = cc / cg
            elif cc > 0:
                speedup = None
                parity_note = "gpu median 0 -- speedup undefined"
        for k, v in legs.items():
            sp = v.get("speedup_gpu_vs_cpu")
            if sp is not None and math.isfinite(sp) and speedup is None:
                speedup = sp
        rows.append(
            {
                "surface": surface,
                "N": n,
                "legs": leg_stats,
                "speedup_gpu_vs_cpu": speedup,
                "speedup_note": parity_note,
                "n_legs": len(legs),
            }
        )
    return rows, src, len(records)


def load_emcee(out_dir):
    """emcee A/B RESULT lines from P3 (primary) and P2 T2 (secondary)."""

    def extract(path):
        text = _read(path)
        if text is None:
            return None
        legs = {}
        sp = None
        for line in text.splitlines():
            m = EMCEE_LEG_RE.search(line)
            if m:
                walls = [float(x) for x in m.group(7).split()]
                legs[m.group(1)] = {
                    "n_points": int(m.group(2)),
                    "ndim": int(m.group(3)),
                    "nwalkers": int(m.group(4)),
                    "steps_run": int(m.group(5)),
                    "median_steps_per_s": _n(m.group(6)),
                    "walls_s": walls,
                }
            m = EMCEE_SP_RE.search(line)
            if m:
                sp = _n(m.group(1))
        if not legs and sp is None:
            return None
        return {"legs": legs, "speedup_gpu_vs_cpu": sp, "log": os.path.basename(path)}

    p3 = extract(os.path.join(out_dir, "P3_emcee.log"))
    p2 = extract(os.path.join(out_dir, "P2_logs", "emcee_T2.log"))
    out = {"P3_standalone": p3, "P2_T2_in_orchestrator": p2, "primary": p3 or p2}
    if out["primary"]:
        pr = out["primary"]
        cpu = pr["legs"].get("cpu")
        gpu = pr["legs"].get("gpu")
        if cpu and gpu and cpu["median_steps_per_s"] and gpu["median_steps_per_s"]:
            pr["derived_speedup_gpu_vs_cpu"] = (
                gpu["median_steps_per_s"] / cpu["median_steps_per_s"]
            )
    return out


def load_p4(out_dir):
    row = {"log": "P4_kepler63b.log"}
    tim = os.path.join(out_dir, "P4_artifacts", "timings.json")
    text = _read(tim)
    if text is not None:
        try:
            d = json.loads(text)
            row.update(
                {
                    "source": "P4_artifacts/timings.json",
                    "gpu": {
                        k: d["gpu"].get(k) for k in ("status", "rc", "wall_s", "pyde_s")
                    },
                    "cpu": {
                        k: d["cpu"].get(k) for k in ("status", "rc", "wall_s", "pyde_s")
                    },
                    "speedup_cpu_over_gpu": d.get("speedup_cpu_over_gpu"),
                    "parity_max_abs_median_diff": (d.get("parity") or {}).get(
                        "max_abs_median_diff"
                    ),
                    "n_params_compared": (d.get("parity") or {}).get(
                        "n_params_compared"
                    ),
                    "gpu_usage": d.get("gpu_usage"),
                    "run_dir": d.get("run_dir"),
                }
            )
            return _p4_artifact_truth(out_dir, row)
        except ValueError:
            row["timings_json_note"] = "present but unparseable"
    text = _read(os.path.join(out_dir, "P4_kepler63b.log"))
    if text is not None:
        g = re.search(K63_LEG_RE_TPL.format(leg="GPU"), text)
        c = re.search(K63_LEG_RE_TPL.format(leg="CPU"), text)
        if g:
            row["gpu"] = {
                "rc": int(g.group(1)),
                "wall_s": _n(g.group(2)),
                "pyde_s": _n(g.group(3)),
                "status": "ok" if g.group(1) == "0" else "FAILED",
            }
        if c:
            row["cpu"] = {
                "rc": int(c.group(1)),
                "wall_s": _n(c.group(2)),
                "pyde_s": _n(c.group(3)),
                "status": "ok" if c.group(1) == "0" else "FAILED",
            }
        if "gpu" not in row and "cpu" not in row:
            row["result"] = "NOT_RUN_OR_ABORTED"
            tail = "\n".join(text.splitlines()[-12:])
            row["log_tail"] = tail
    skip_reason = _p4_skip_reason(out_dir)
    if skip_reason:
        row["result"] = "SKIPPED"
        row["skip_reason"] = skip_reason
    return _p4_artifact_truth(out_dir, row)


def _p4_artifact_truth(out_dir, row):
    """Set the artifact-copy truth bits on ``row`` (consumed by HONEST_STATS §4).

    The driver (``gpu_clean_run.sh`` P4 section) creates
    ``$OUT/P4_artifacts/`` and fills it ONLY when the author script reaches
    its DONE block — the block that announces a ``run dir     : <path>``
    line in the leg log and where the leg copyable files live. So a
    "artifacts copied" claim in the report is true ONLY when BOTH:
    the directory exists in OUT AND is non-empty, AND a run dir was
    announced (in ``timings.json`` or the log line above). When either is
    missing (typical: the leg hit its timeout budget before the DONE
    block), the report must say NOT copied — never claim a copy and never
    emit a stale ``run dir: ?``.
    """
    if row.get("run_dir") is None:
        log_text = _read(os.path.join(out_dir, "P4_kepler63b.log"))
        if log_text:
            m = re.search(r"^\s*run dir\s*:\s*(\S+)", log_text, re.MULTILINE)
            if m:
                row["run_dir"] = m.group(1)
    art_dir = os.path.join(out_dir, "P4_artifacts")
    present = os.path.isdir(art_dir)
    nonempty = False
    if present:
        try:
            nonempty = bool(os.listdir(art_dir))
        except OSError:
            nonempty = False
    row["artifacts_dir_present"] = present
    row["run_dir_announced"] = row.get("run_dir") is not None
    row["artifacts_copied"] = bool(nonempty and row.get("run_dir"))
    return row


def _p4_skip_reason(out_dir):
    flags = parse_flags(os.path.join(out_dir, "flags.env"))
    if flags.get("P4_SKIPPED") == "1":
        return flags.get("P4_SKIP_REASON") or "(see flags.env)"
    for line in _read(os.path.join(out_dir, "legs.tsv") or "").splitlines():
        parts = line.split("\t")
        # legs.tsv columns: name rc wall_s timed_out state note log
        if (
            parts
            and parts[0] == "kepler63b_ab"
            and len(parts) >= 6
            and parts[4] == "skipped"
        ):
            return parts[5]
    return None


# ------------------------------------------------------------------- flags --
def parse_flags(path):
    out = {}
    text = _read(path)
    if text is None:
        return out
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        out[k.strip()] = v.strip().strip('"')
    return out


def load_legs(out_dir):
    legs = []
    text = _read(os.path.join(out_dir, "legs.tsv"))
    if text is None:
        return legs
    lines = [l for l in text.splitlines() if l.strip()]
    if not lines:
        return legs
    header = lines[0].split("\t")
    for line in lines[1:]:
        parts = line.split("\t")
        row = dict(zip(header, parts))
        legs.append(row)
    return legs


def load_contention(out_dir):
    snaps = []
    for rel in ("P2_logs/contention.log",):
        text = _read(os.path.join(out_dir, rel))
        if text is None:
            continue
        for line in text.splitlines():
            if line.strip():
                snaps.append({"log": rel, "line": line.strip()})
    return snaps[-40:]


def load_p0_meta(out_dir):
    meta = {}
    flags = parse_flags(os.path.join(out_dir, "flags.env"))
    meta["flags"] = flags
    text = _read(os.path.join(out_dir, "P0_env.log"))
    if text is None:
        meta["p0_log_missing"] = True
        return meta
    mm = re.search(r"GPU MODEL:\s*(.+)", text)
    if mm:
        meta["gpu_model"] = mm.group(1).strip()
    return meta


# ------------------------------------------------------------------ proven --
def provenance(out_dir):
    prov = {}
    rels = [
        "flags.env",
        "legs.tsv",
        "P0_env.log",
        "P1_parity_pyorbit.log",
        "P1_parity_pyorbit_gpu.log",
        "P1_b1_transit.log",
        "P1_b3_esp.log",
        "P1_b4_apodized.log",
        "P1_b5_trend.log",
        "P1_b6_tinygp.log",
        "P2_perf.log",
        "P2_logs/bench.log",
        "P2_logs/results.json",
        "P2_logs/contention.log",
        "P2_logs/emcee_T2.log",
        "P3_emcee.log",
        "P4_kepler63b.log",
        "P4_artifacts/timings.json",
        "P4_artifacts/report_raw.md",
        "HONEST_STATS.md",
        "stats.json",
    ]
    for rel in rels:
        p = os.path.join(out_dir, rel)
        if not os.path.isfile(p):
            continue
        try:
            size = os.path.getsize(p)
        except OSError:
            continue
        entry = {"size_bytes": size, "sha256_16": _hash(p)}
        try:
            with open(p, encoding="utf-8", errors="replace") as fh:
                entry["lines"] = sum(1 for _ in fh)
        except OSError:
            pass
        prov[rel] = entry
    return prov


# ------------------------------------------------------------------ report --
def _ref_cell(row):
    ref = row.get("reference")
    if not ref:
        return "-"
    if ref.get("summary"):
        return ref["summary"]
    if ref.get("worst_delta") is not None:
        d = ref["worst_delta"]
        band = ref.get("band")
        return "{0} (band {1})".format(_f(d), _f(band))
    return "-"


def _match_cell(row):
    if row.get("missing"):
        return "not run"
    if row.get("result") == "NOT_RUN":
        return "not run"
    m = row.get("matches_reference")
    if m is None:
        return "-"
    return "matches contested ref" if m else "DIFFERS from contested ref"


def build_markdown(ctx):
    L = []
    a = L.append
    flags = ctx["meta"].get("flags", {})
    contended = flags.get("CONTENDED") == "1"
    low_vram = flags.get("LOW_VRAM") == "1"

    a("# PyORBIT GPU — CLEAN RUN — honest performance & precision statistics")
    a("")
    a("Generated: {0}".format(ctx["generated_utc"]))
    a("Run dir (OUT): `{0}`".format(ctx["out_dir"]))
    if ctx["meta"].get("gpu_model"):
        a("GPU: {0}".format(ctx["meta"]["gpu_model"]))
    a("")
    # ---- caveats first ----
    a("## 0. Environment & contention (read first)")
    a("")
    if contended:
        a(
            "> **CAVEAT: the run is marked CONTENDED** — ollama's "
            "local API (http://127.0.0.1:11434/api/tags) was reachable at P0. "
            "GPU wall times below may be inflated or skewed; treat speedups as "
            "lower/upper bounds, not point estimates. Precision (value parity) "
            "is unaffected by contention, only performance is."
        )
    else:
        a(
            "Contention check: ollama API **not** reachable at P0 → card assumed "
            "uncontended (`CONTENDED=0`). No foreign process was touched; this "
            "script only probes, never kills."
        )
    a("")
    a(
        "- Free VRAM at P0: {0} GiB of {1} GiB total".format(
            _n(flags.get("FREE_GIB")), _n(flags.get("GPU_TOTAL_GIB"))
        )
    )
    a(
        "- `LOW_VRAM={0}` ({1})".format(
            flags.get("LOW_VRAM", "n/a"),
            "the P4 Kepler-63b A/B is skipped by policy when free < 20 GiB"
            if low_vram
            else "not set",
        )
    )
    if flags.get("PY_ORBIT_ENV") or flags.get("PYORBIT_ENV"):
        a(
            "- Conda envs: CPU=`{0}`, GPU=`{1}` (versions in `P0_env.log`)".format(
                flags.get("PYORBIT_ENV", "pyorbit"), flags.get("GPU_ENV", "pyorbit-gpu")
            )
        )
    a("- Reference session: contested numbers from {0}".format(CONTESTED_SRC))
    a("")

    # ---- precision ----
    a("## 1. Precision (clean-env re-confirmation)")
    a("")
    a(
        "| surface | verdict | worst-Δ | band (tolerance) | contested ref | match? | log |"
    )
    a("|---|---|---|---|---|---|---|")
    for row in ctx["precision"]:
        verdict = row["result"] if row["result"] else "-"
        if row.get("summary") and row["surface"].startswith("parity-suite"):
            verdict = "{0} ({1})".format(verdict, row["summary"])
        a(
            "| {0} | {1} | {2} | {3} | {4} | {5} | `{6}` |".format(
                _mdcell(row["surface"]),
                _mdcell(verdict),
                _mdcell(_f(row["worst_delta"])),
                _mdcell(_f(row["band"])),
                _mdcell(_ref_cell(row)),
                _mdcell(_match_cell(row)),
                _mdcell(row["log"]),
            )
        )
    a("")
    any_missing = [r["surface"] for r in ctx["precision"] if r.get("missing")]
    if any_missing:
        a(
            "Not-run in this clean run: {0} (leg never produced output; see leg "
            "ledger in §7).".format(", ".join(any_missing))
        )
        a("")

    # ---- performance ----
    a("## 2. Performance (clean env, GPU vs CPU)")
    a("")
    if ctx["performance"]:
        a(
            "| surface | N | leg | reps | warmup | cold ms | median ms | "
            "mean±std ms | min–max ms | IQR ms |"
        )
        a("|---|---|---|---|---|---|---|---|---|---|")
        for row in ctx["performance"]:
            legs = row["legs"]
            first = True
            for leg_name in ("cpu", "gpu", "gpu-jitref"):
                if leg_name not in legs:
                    continue
                s = legs[leg_name]
                mean_std = (
                    "{0}±{1}".format(
                        _f(s["mean_ms"], "{0:.2f}"), _f(s["std_ms"], "{0:.2f}")
                    )
                    if s["mean_ms"] is not None
                    else "-"
                )
                a(
                    "| {0} | {1} | {2} | {3} | {4} | {5} | {6} | {7} | {8} | {9} |".format(
                        row["surface"] if first else "",
                        row["N"] if first else "",
                        leg_name,
                        s["n_reps"],
                        s["warmup"],
                        _f(s["cold_ms"], "{0:.2f}"),
                        _f(s["median_ms"], "{0:.2f}"),
                        mean_std,
                        "{0}–{1}".format(
                            _f(s["min_ms"], "{0:.2f}"), _f(s["max_ms"], "{0:.2f}")
                        ),
                        _f(s["iqr_ms"], "{0:.2f}"),
                    )
                )
                first = False
            sp = row.get("speedup_gpu_vs_cpu")
            sp_txt = "{0:.2f}×".format(sp) if sp else "-"
            a(
                "| | | **speedup (CPU med / GPU med)** | | | | | **{0}**{1} | | |".format(
                    sp_txt,
                    " ({0})".format(row["speedup_note"])
                    if row.get("speedup_note")
                    else "",
                )
            )
        a("")
        a(
            "Source: `{0}` (last session); per-rep arrays in `stats.json` "
            "(`performance[*].legs[*].times_ms`).".format(ctx["performance_src"])
        )
        guards = []
        for r in ctx["performance"]:
            for _leg, _stats in r["legs"].items():
                g = _stats.get("parity_guard")
                if g is not None:
                    guards.append((r["surface"], r["N"], g))
        if guards:
            mx = max(g for _, _, g in guards)
            a("")
            a(
                "In-run value-parity guards (enforced per leg, from `bench.log` "
                "sanity notes): worst |Δ| across all measured legs = **{0}** "
                "(all finite; a leg that diverged in value would have failed the "
                "leg's own sanity check and been marked `sanity-fail`).".format(_f(mx))
            )
    else:
        a(
            "**No performance data.** Either P2 was skipped or it produced no "
            "bench records before failing; see the leg ledger in §7 and "
            "`P2_perf.log` / `P2_logs/` for the exact reason. Nothing is "
            "invented here."
        )
    a("")

    # ---- emcee ----
    a("## 3. emcee A/B (P3 standalone leg, clean env)")
    a("")
    em = ctx["emcee"]
    prim = em.get("primary")
    if prim:
        for leg_name in ("cpu", "gpu"):
            s = prim["legs"].get(leg_name)
            if not s:
                continue
            a(
                "- **{0} leg**: median {1} steps/s over {2} timed runs "
                "({3} steps × {4} walkers, ndim={5}; walls: {6} s)".format(
                    leg_name,
                    _f(s["median_steps_per_s"], "{0:.2f}"),
                    len(s["walls_s"]),
                    s["steps_run"] // max(len(s["walls_s"]), 1),
                    s["nwalkers"],
                    s["ndim"],
                    " / ".join(_f(w, "{0:.2f}") for w in s["walls_s"]),
                )
            )
        sp = prim.get("speedup_gpu_vs_cpu") or prim.get("derived_speedup_gpu_vs_cpu")
        if sp:
            a(
                "- **speedup (GPU steps/s / CPU steps/s): {0:.3f}×** — honest "
                "note: emcee's per-eval work is small here, so the GPU's fixed "
                "launch overhead dominates; a speedup near 1× (or < 1×) at this "
                "scale is expected, not a failure.".format(sp)
            )
        a("- Source: `{0}`".format(prim["log"]))
        if (
            em.get("P3_standalone")
            and em.get("P2_T2_in_orchestrator")
            and em["P2_T2_in_orchestrator"] is not em["P3_standalone"]
        ):
            p2x = em["P2_T2_in_orchestrator"]
            a(
                "- (P2's in-orchestrator T2 emcee also produced: speedup={0})".format(
                    _f(p2x.get("speedup_gpu_vs_cpu"))
                )
            )
    else:
        a(
            "**Not run** — no emcee RESULT lines found in `P3_emcee.log` or "
            "`P2_logs/emcee_T2.log`. See §7 leg ledger for the recorded rc/timeout."
        )
    a("")

    # ---- P4 ----
    a("## 4. Kepler-63b A/B (P4, author's script run verbatim)")
    a("")
    p4 = ctx["p4"]
    if p4.get("result") == "SKIPPED":
        a("**SKIPPED** — {0}".format(p4.get("skip_reason", "no reason recorded")))
        a("")
        a(
            "(The pre-leg free-VRAM guard is satisfied by the author's own "
            "fail-fast: `MIN_FREE_MIB=8192` in `run_kepler63b_gpu_test.sh`, "
            "which hard-exits if the card is not free. We never disable it and "
            "never kill foreign processes.)"
        )
        a("")
    elif p4.get("result") == "NOT_RUN_OR_ABORTED" or (
        not p4.get("gpu") and not p4.get("cpu") and "source" not in p4
    ):
        a(
            "**Not completed** (aborted or timed out inside the author's "
            "script). Captured tail of `P4_kepler63b.log`:"
        )
        a("")
        a("```")
        for line in (p4.get("log_tail") or "(no log)").splitlines():
            a(line)
        a("```")
        a("")
    else:
        g, c = p4.get("gpu") or {}, p4.get("cpu") or {}
        a("| leg | status | rc | wall s | PyDE s |")
        a("|---|---|---|---|---|")
        for name, leg in (("GPU (gpu:true)", g), ("CPU (gpu:false)", c)):
            a(
                "| {0} | {1} | {2} | {3} | {4} |".format(
                    name,
                    leg.get("status", "-"),
                    leg.get("rc", "-"),
                    _f(leg.get("wall_s"), "{0:.1f}"),
                    _f(leg.get("pyde_s"), "{0:.1f}"),
                )
            )
        sp = p4.get("speedup_cpu_over_gpu")
        if sp:
            a("")
            a(
                "**CPU wall / GPU wall = {0:.2f}×** (wall-clock ratio, as "
                "computed by the author's `_collect.py`; >1 means the CPU leg "
                "was *slower* overall, i.e. the GPU won end-to-end; <1 the "
                "opposite).".format(sp)
            )
        if p4.get("parity_max_abs_median_diff") is not None:
            a("")
            a(
                "Posterior parity: max |median diff| GPU vs CPU = **{0}** "
                "over {1} parameters (sampler chains differ by seed/proposal "
                "path; medians compared, not samples).".format(
                    _f(p4["parity_max_abs_median_diff"]),
                    p4.get("n_params_compared", "?"),
                )
            )
        a("")
        cfg = ctx["meta"].get("flags", {}).get("P4_CONFIG", "")
        a(
            "Config source: `{0}` (the clean-run orchestrator defaults to the "
            "author's *capped benchmark* config — PyDE ngen=100, emcee "
            "nsteps=500 — because the production config with nsteps=100000 "
            "runs for hours and would break the ≤60 min budget).".format(
                cfg or "see P4_kepler63b.log"
            )
        )
        if p4.get("artifacts_copied"):
            a(
                "Full artifacts copied to `P4_artifacts/` (run dir: {0}).".format(
                    p4.get("run_dir")
                )
            )
        else:
            a(
                "Artifacts NOT copied (leg timed out / DONE block not "
                "reached; run dir never announced in the leg log — no "
                "copy to claim for this run; see §7 leg ledger for the "
                "leg's rc/timeout)."
            )
    a("")

    # ---- contention evidence ----
    a("## 5. Contention evidence (nvidia-smi)")
    a("")
    snaps = ctx["contention"]
    if snaps:
        a(
            "Last nvidia-smi snapshots recorded by the P2 orchestrator (`{0}`):".format(
                snaps[0]["log"]
            )
        )
        a("")
        for s in snaps[-12:]:
            a("> {0}".format(s["line"]))
        a("")
        a("Full P0 audit (nvidia-smi + ollama probe + pip versions): `P0_env.log`.")
    else:
        a(
            "No `contention.log` snapshots found (P2 orchestrator never ran or "
            "logdir differs). P0 audit: see `P0_env.log`."
        )
    a("")

    # ---- provenance ----
    a("## 6. Reproducibility / provenance")
    a("")
    a(
        "All outputs of this report derive ONLY from files in the run dir "
        "(sha256-16 below, sizes and line counts in `stats.json.provenance`):"
    )
    a("")
    for rel, entry in sorted(ctx["provenance"].items()):
        a(
            "- `{0}` — {1} B, {2} lines, sha256 `{3}`".format(
                rel,
                entry["size_bytes"],
                entry.get("lines", "?"),
                entry.get("sha256_16", "?"),
            )
        )
    a("")
    a("Commands that produced these (see `gpu_clean_run.sh`):")
    a("")
    for cmd in ctx["commands"]:
        a("- `{0}`".format(cmd))
    a("")

    # ---- leg ledger ----
    a("## 7. Leg ledger")
    a("")
    legs = ctx["legs"]
    if legs:
        a("| leg | rc | wall s | timed_out | state | note | log |")
        a("|---|---|---|---|---|---|---|")
        for leg in legs:
            a(
                "| {0} | {1} | {2} | {3} | {4} | {5} | `{6}` |".format(
                    _mdcell(leg.get("name", "?")),
                    _mdcell(leg.get("rc", "-")),
                    _mdcell(leg.get("wall_s", "-")),
                    _mdcell(leg.get("timed_out", "-")),
                    _mdcell(leg.get("state", "-")),
                    _mdcell(leg.get("note", "-") or "-"),
                    _mdcell(leg.get("log", "-")),
                )
            )
    else:
        a("No leg ledger recorded.")
    a("")
    a("---")
    a(
        "_Generated by `development/gpu_perf_stats_report.py` (stdlib only). "
        "If any row above says not-run/skipped, that is a fact, not a "
        "placeholder — the corresponding log path is given._"
    )
    a("")
    return "\n".join(L)


def build_summary(ctx):
    """Compact pasteable console summary (P6 of the clean run)."""
    flags = ctx["meta"].get("flags", {})
    lines = []
    a = lines.append
    a("== CLEAN RUN COMPLETE ==")
    cont = flags.get("CONTENDED") == "1"
    lv = flags.get("LOW_VRAM") == "1"
    a("UTC: {0} | OUT: {1}".format(ctx["generated_utc"], ctx["out_dir"]))
    a(
        "CONTENDED: {0} (ollama reachable at P0: {1}) | LOW_VRAM: {2}"
        " | free VRAM at P0: {3} GiB".format(
            "YES *** WARNING ***" if cont else "no",
            flags.get("OLLAMA_REACHABLE", "?"),
            "yes" if lv else "no",
            flags.get("FREE_GIB", "?"),
        )
    )
    if cont:
        a(
            "WARNING: card was contended at P0 — GPU wall times are NOT honest "
            "point estimates; precision results are still valid."
        )
    a("-- precision (clean env) --")
    for row in ctx["precision"]:
        if row.get("missing"):
            a("  {0:<28s} NOT RUN (leg did not complete)".format(row["surface"]))
            continue
        delta = (
            _f(row["worst_delta"])
            if row["worst_delta"] is not None
            else (row.get("summary") or "-")
        )
        a(
            "  {0:<28s} {1:<9s} worst-Δ {2}  [contested ref: {3}, {4}]".format(
                row["surface"],
                row["result"] or "?",
                delta,
                _ref_cell(row),
                _match_cell(row),
            )
        )
    a("-- performance (median ms; speedup = CPU med / GPU med) --")
    n_perf = 0
    for row in ctx["performance"]:
        cpu = row["legs"].get("cpu")
        gpu = row["legs"].get("gpu")
        if not cpu or not gpu:
            continue
        sp = row.get("speedup_gpu_vs_cpu")
        a(
            "  {0:<8s} N={1:<6d} CPU {2:>10s} ms | GPU {3:>10s} ms | "
            "speedup {4:<6s} | reps {5}/{6}".format(
                row["surface"],
                row["N"],
                _f(cpu["median_ms"], "{0:.1f}"),
                _f(gpu["median_ms"], "{0:.1f}"),
                "{0:.2f}x".format(sp) if sp else "-",
                cpu["n_reps"],
                gpu["n_reps"],
            )
        )
        n_perf += 1
    if n_perf == 0 and not ctx["performance"]:
        a("  (no perf data — see HONEST_STATS.md §2 and leg ledger)")
    em = (ctx["emcee"] or {}).get("primary")
    if em:
        sp = em.get("speedup_gpu_vs_cpu") or em.get("derived_speedup_gpu_vs_cpu")
        cpu = em["legs"].get("cpu") or {}
        gpu = em["legs"].get("gpu") or {}
        a("-- emcee A/B (P3) --")
        a(
            "  CPU {0} steps/s | GPU {1} steps/s | speedup {2} (n={3})".format(
                _f(cpu.get("median_steps_per_s"), "{0:.2f}"),
                _f(gpu.get("median_steps_per_s"), "{0:.2f}"),
                _f(sp, "{0:.3f}") if sp else "-",
                (cpu.get("n_points") or gpu.get("n_points") or "?"),
            )
        )
    else:
        a("-- emcee A/B (P3) --  NOT RUN (see leg ledger)")
    p4 = ctx["p4"]
    a("-- Kepler-63b A/B (P4, author script verbatim) --")
    if p4.get("result") == "SKIPPED":
        a("  SKIPPED: {0}".format(p4.get("skip_reason", "?")))
    elif p4.get("gpu") and p4.get("cpu"):
        g, c = p4["gpu"], p4["cpu"]
        sp = p4.get("speedup_cpu_over_gpu")
        a(
            "  GPU wall {0} s (pyde {1}) | CPU wall {2} s (pyde {3}) | "
            "CPU/GPU wall {4}".format(
                _f(g.get("wall_s"), "{0:.1f}"),
                _f(g.get("pyde_s"), "{0:.1f}"),
                _f(c.get("wall_s"), "{0:.1f}"),
                _f(c.get("pyde_s"), "{0:.1f}"),
                _f(sp, "{0:.2f}") + "x" if sp else "-",
            )
        )
    else:
        a(
            "  NOT COMPLETED: {0}".format(
                p4.get("skip_reason") or p4.get("result") or "see P4_kepler63b.log"
            )
        )
    timeouts = [l for l in ctx["legs"] if str(l.get("timed_out")) == "1"]
    failed = [l for l in ctx["legs"] if str(l.get("rc")) not in ("0", "", "-")]
    if timeouts:
        a(
            "TIMEOUTS hit: {0}".format(
                ", ".join(
                    "{0}({1}s)".format(l.get("name"), l.get("wall_s")) for l in timeouts
                )
            )
        )
    if failed and not timeouts:
        a(
            "Non-zero rc legs: {0}".format(
                ", ".join(
                    "{0}(rc={1})".format(l.get("name"), l.get("rc")) for l in failed
                )
            )
        )
    skipped = [
        l
        for l in ctx["legs"]
        if str(l.get("state")) == "skipped" and str(l.get("timed_out")) != "1"
    ]
    if skipped:
        a(
            "Skipped legs: {0}".format(
                ", ".join(
                    "{0}({1})".format(l.get("name"), (l.get("note") or "")[:40])
                    for l in skipped
                )
            )
        )
    a("OUT={0}  (HONEST_STATS.md, stats.json)".format(ctx["out_dir"]))
    return "\n".join(lines)


# --------------------------------------------------------------------- main --
def collect_commands(out_dir):
    """Best-effort: recover the exact per-leg commands from the orchestrator
    's P0_env.log / legs.tsv so the report cites the real invocations."""
    cmds = [
        "bash development/gpu_clean_run.sh   (driver; per-leg commands below)",
    ]
    text = _read(os.path.join(out_dir, "P0_env.log"))
    if text is None:
        return cmds
    m = re.search(r"--- commands ---\n(.*?)(?:\n---|\Z)", text, re.DOTALL)
    if m:
        for line in m.group(1).splitlines():
            line = line.strip()
            if line:
                cmds.append(line)
    return cmds


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True, help="clean-run OUT directory")
    ap.add_argument(
        "--print-summary",
        action="store_true",
        help="also print the compact console summary",
    )
    args = ap.parse_args(argv)

    out_dir = os.path.abspath(args.out)
    if not os.path.isdir(out_dir):
        print("ERROR: OUT directory not found: {0}".format(out_dir), file=sys.stderr)
        return 2

    ctx = {
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime(
            "%Y-%m-%d %H:%M:%S UTC"
        ),
        "out_dir": out_dir,
        "meta": load_p0_meta(out_dir),
        "precision": [
            parse_precision_log(os.path.join(out_dir, rel), name)
            for rel, name in PRECISION_SOURCES
        ],
        "emcee": load_emcee(out_dir),
        "p4": load_p4(out_dir),
        "contention": load_contention(out_dir),
        "provenance": provenance(out_dir),
        "commands": collect_commands(out_dir),
        "legs": load_legs(out_dir),
    }
    perf, perf_src, n_recs = build_performance(out_dir)
    ctx["performance"] = perf
    ctx["performance_src"] = perf_src or "(none found)"
    ctx["n_bench_records"] = n_recs

    md_path = os.path.join(out_dir, "HONEST_STATS.md")
    js_path = os.path.join(out_dir, "stats.json")
    with open(md_path, "w", encoding="utf-8") as fh:
        fh.write(build_markdown(ctx))
    ctx["provenance"] = provenance(out_dir)  # include the two outputs' hashes
    with open(js_path, "w", encoding="utf-8") as fh:
        json.dump(ctx, fh, indent=2, sort_keys=True, default=str)

    print("wrote {0}".format(md_path))
    print("wrote {0}".format(js_path))
    if args.print_summary:
        print()
        print(build_summary(ctx))
    return 0


if __name__ == "__main__":
    sys.exit(main())
