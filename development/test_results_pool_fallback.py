#!/usr/bin/env python3
"""Deterministic repro test for the pyorbit_results pool-hang fix.

Pins the production failure mode of ``pyorbit.samplers.pyorbit_getresults``:
a spawned results worker that dies while rebuilding its state (in the
production run the child died unpicking its first payload with

    unable to create StreamExecutor for CUDA:0 ... CUDA_ERROR_OUT_OF_MEMORY

). In the affected code the parent then blocked forever inside the pool,
leaked pool semaphores, and never emitted the documented
``WARNING: parallel ... failed (...); retrying serially`` fallback.

This script simulates worker death in two ways (no GPU, no OOM, fully
deterministic):
  1. the child dies while unpickling the initargs / running the initializer
    (``DeadOnUnpickle``) -- the reported hang mode;
  2. the child dies mid-task (``DyingMidTask``, ``os._exit``) -- the
    in-flight-worker mode.

Acceptance assertions per case:
  (a) the parent does NOT block unboundedly (bounded by the module budgets,
      plus a hard SIGALRM test-level timeout so CI can never hang);
  (b) the explicit ``WARNING: parallel ... failed ... retrying serially``
      fallback fires;
  (c) the step still completes with EXACT results (serial path) and leaves
      no leaked worker processes or pool threads behind.

A happy-path case proves the parallel pool still works and yields the same
exact numbers.

Run standalone (this repo has no pytest configured; the script is plain
python and is also pytest-collectable if pytest is ever added):

    python3 development/test_results_pool_fallback.py
"""

import contextlib
import io
import os
import signal
import sys
import threading
import time

sys.path.insert(
    0,
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
)

from pyorbit.samplers import pyorbit_getresults as results_mod  # noqa: E402

# ---------------------------------------------------------------------------
# Model-container stand-ins (module-level so they are picklable by reference
# in spawned children)
# ---------------------------------------------------------------------------


class FakeMC:
    """Harmless picklable stand-in for the model container."""

    def nautilus_call(self, theta):
        return -10.0 * (1 + int(theta))


class DeadOnUnpickle:
    """Parent-side it works fine; a child dies the moment it rebuilds it.

    Mirrors the production child death (CUDA OOM while re-initializing in
    the child) by failing inside ``__setstate__``, i.e. during unpickle of
    the worker's initargs -- before the worker does any work.

    NOTE: the instance must carry non-empty state. This platform's pickle
    omits the state-rebuild step for empty-``__dict__`` instances (verified
    in-process: an empty instance round-trips without calling
    ``__setstate__`` at all), so a stateless death hook would never fire in
    the child and the "death" case would silently become a happy path.
    """

    def __init__(self):
        self._state = 1  # non-empty __dict__ so unpickle runs __setstate__

    def nautilus_call(self, theta):
        return -10.0 * (1 + int(theta))

    def __setstate__(self, state):
        raise RuntimeError(
            "simulated results worker death: unable to re-initialize in child (OOM)"
        )


class SlowInitMC:
    """Parent-side fine; a child hangs forever rebuilding it.

    ``__setstate__`` blocks past any sane budget, exercising the preflight
    timeout + bounded terminate/kill branch (a child that OOM-hangs inside
    JAX/CUDA re-init in the reported scenario).
    """

    def __init__(self):
        self._state = 1  # non-empty state so unpickle runs __setstate__

    def nautilus_call(self, theta):
        return -10.0 * (1 + int(theta))

    def __setstate__(self, state):
        import time

        time.sleep(999)  # longer than any test budget; SIGTERM kills it


class DyingMidTask:
    """Worker dies mid-task (``os._exit``); parent-side use is fine."""

    def nautilus_call(self, theta):
        import multiprocessing

        if multiprocessing.current_process().name != "MainProcess":
            os._exit(9)
        return -10.0 * (1 + int(theta))


class TrapTermMC:
    """Worker traps SIGTERM (child-side) and then hangs mid-task (B-F1).

    Pins the blocked teardown defect: a worker that ignores SIGTERM must
    still be SIGKILLed by the bounded reap, so ``pool.terminate()`` -- whose
    internal teardown in this stdlib runs synchronously on the caller thread
    and ends in a timeoutless worker join after SIGTERing the workers --
    cannot block the parent.
    Parent-side use is fine (the real MainProcess just gets its number).
    """

    def __init__(self):
        self._state = 1  # non-empty state so unpickle runs __setstate__

    def nautilus_call(self, theta):
        import multiprocessing
        import time

        if multiprocessing.current_process().name != "MainProcess":
            time.sleep(999)  # hang mid-task while SIGTERM is ignored
        return -10.0 * (1 + int(theta))

    def __setstate__(self, state):
        import signal

        signal.signal(signal.SIGTERM, signal.SIG_IGN)


# ---------------------------------------------------------------------------
# Test scaffolding
# ---------------------------------------------------------------------------


class HangGuard:
    """Test-level timeout (SIGALRM): CI can never hang, whatever happens."""

    def __init__(self, seconds):
        self.seconds = seconds

    def __enter__(self):
        if hasattr(signal, "SIGALRM"):
            signal.signal(signal.SIGALRM, self._on_alarm)
            signal.alarm(self.seconds)
        return self

    def __exit__(self, *exc):
        if hasattr(signal, "SIGALRM"):
            signal.alarm(0)
        return False

    @staticmethod
    def _on_alarm(signum, frame):
        raise TimeoutError(
            "TEST HUNG: parent blocked unboundedly inside the results pool "
            "(acceptance (a) violated)"
        )


def own_child_pids():
    """PIDs of this process's real children (empty <=> no leaked workers).

    Excludes the stdlib multiprocessing *resource tracker* daemon: the first
    spawned process starts it as a long-lived child of this process (its
    cmdline is ``...multiprocessing.resource_tracker import main;main(N)``),
    and it legitimately outlives pool teardown. Real pool workers never
    match that cmdline, so a genuinely leaked worker is still flagged.
    """
    me = os.getpid()
    pids = []
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            with open(os.path.join("/proc", entry, "stat"), "r") as fh:
                data = fh.read()
        except OSError:
            continue
        try:
            ppid = int(data.rsplit(")", 1)[1].split()[1])
        except (IndexError, ValueError):
            continue
        if ppid != me:
            continue
        try:
            with open(os.path.join("/proc", entry, "cmdline"), "r") as fh:
                cmd = fh.read()
        except OSError:
            cmd = ""
        if "resource_tracker" in cmd:
            continue  # legitimate stdlib daemon, not a leaked worker
        pids.append(int(entry))
    return pids


def expected_values(mc, items):
    return [mc.nautilus_call(x) for x in items]


def run_case(mc, items, workers=2, mp_method="spawn"):
    """Run one ``_results_parallel_map`` call with small, CI-safe budgets."""
    original = (
        results_mod._RESULTS_WORKER_INIT_TIMEOUT,
        results_mod._RESULTS_ITEM_TIMEOUT,
    )
    results_mod._RESULTS_WORKER_INIT_TIMEOUT = 15.0
    results_mod._RESULTS_ITEM_TIMEOUT = 5.0
    baseline_threads = len(threading.enumerate())
    start = time.monotonic()
    try:
        with HangGuard(90), contextlib.redirect_stdout(io.StringIO()) as out:
            got = results_mod._results_parallel_map(
                mc, "nautilus", items, workers, mp_method
            )
    finally:
        (
            results_mod._RESULTS_WORKER_INIT_TIMEOUT,
            results_mod._RESULTS_ITEM_TIMEOUT,
        ) = original
    elapsed = time.monotonic() - start
    text = out.getvalue()
    warnings = [line for line in text.splitlines() if "WARNING: parallel" in line]
    # give teardown a moment to settle, then assert no leaks
    time.sleep(2.0)
    leaked_pids = own_child_pids()
    leaked_threads = len(threading.enumerate()) - baseline_threads
    return {
        "got": list(got),
        "warnings": warnings,
        "elapsed": elapsed,
        "leaked_pids": leaked_pids,
        "leaked_threads": leaked_threads,
    }


def check_case(name, mc, items, expect_fallback):
    result = run_case(mc, items)
    expected = expected_values(mc, items)

    # (c) exact results: serial fallback must keep identical numbers
    assert result["got"] == expected, (
        "{0}: results differ. got={1!r} expected={2!r}".format(
            name, result["got"], expected
        )
    )

    # (b) explicit fallback warning fires (or is absent on the happy path)
    if expect_fallback:
        assert len(result["warnings"]) == 1, (
            "{0}: expected exactly one fallback warning, got: {1!r}".format(
                name, result["warnings"]
            )
        )
        assert "failed" in result["warnings"][0], result["warnings"]
        assert "retrying serially" in result["warnings"][0], result["warnings"]
    else:
        assert not result["warnings"], "{0}: unexpected fallback warning: {1!r}".format(
            name, result["warnings"]
        )

    # (a) bounded wait: well under the hard test-level watchdog (90 s)
    assert result["elapsed"] < 60.0, (
        "{0}: took {1:.1f}s; the wait is not bounded".format(name, result["elapsed"])
    )

    # (c) no leaked worker processes / pool threads
    assert result["leaked_pids"] == [], "{0}: leaked worker processes: {1}".format(
        name, result["leaked_pids"]
    )
    assert result["leaked_threads"] <= 0, "{0}: leaked {1} threads".format(
        name, result["leaked_threads"]
    )
    return result


# ---------------------------------------------------------------------------
# Test cases (pytest-collectable test_* functions)
# ---------------------------------------------------------------------------

ITEMS = list(range(5))


def test_parallel_happy_path_spawn():
    """The normal parallel pool still works and yields exact results."""
    result = check_case("happy-path", FakeMC(), ITEMS, expect_fallback=False)
    print("    happy-path: {0:.1f}s".format(result["elapsed"]))


def test_worker_death_in_child_initializer_falls_back():
    """Reported hang mode: child dies during initargs unpickle/initializer.

    (a) parent must not block unboundedly, (b) fallback warning fires,
    (c) exact serial results, no leaked workers/threads.
    """
    result = check_case(
        "initializer-death", DeadOnUnpickle(), ITEMS, expect_fallback=True
    )
    print(
        "    initializer-death: {0:.1f}s, warning: {1!r}".format(
            result["elapsed"], result["warnings"]
        )
    )


def test_worker_death_mid_task_falls_back():
    """In-flight worker death: the unfilled slot must time out, not hang."""
    result = check_case("mid-task-death", DyingMidTask(), ITEMS, expect_fallback=True)
    print(
        "    mid-task-death: {0:.1f}s, warning: {1!r}".format(
            result["elapsed"], result["warnings"]
        )
    )


def test_worker_hang_in_initializer_falls_back():
    """Preflight hang: the child blocks rebuilding state past the init
    budget; it must be bounded-terminated/killed, the warning fired, and
    the exact serial results still returned.
    """
    result = check_case("initializer-hang", SlowInitMC(), ITEMS, expect_fallback=True)
    print(
        "    initializer-hang: {0:.1f}s, warning: {1!r}".format(
            result["elapsed"], result["warnings"]
        )
    )


def test_sigterm_trapping_worker_falls_back_bounded():
    """B-F1: a SIGTERM-trapping worker that hangs mid-task must still be
    reaped (SIGKILL) during teardown; the parent returns within budgets
    with the single fallback warning and exact serial results.
    """
    result = check_case(
        "sigterm-trap-mid-task", TrapTermMC(), ITEMS, expect_fallback=True
    )
    print(
        "    sigterm-trap-mid-task: {0:.1f}s, warning: {1!r}".format(
            result["elapsed"], result["warnings"]
        )
    )


def _summary():
    lines = [
        "SUMMARY",
        "  (a) no unbounded parent wait (watchdog 90 s; budgets 15 s / 5 s)",
        "  (b) 'WARNING: parallel ... failed ... retrying serially' fires",
        "  (c) exact results, no leaked worker processes/threads",
    ]
    return "\n".join(lines)


def main():
    cases = [
        ("parallel happy path (spawn)", test_parallel_happy_path_spawn),
        (
            "worker death in child initializer -> serial fallback",
            test_worker_death_in_child_initializer_falls_back,
        ),
        (
            "worker death mid-task -> bounded timeout -> serial fallback",
            test_worker_death_mid_task_falls_back,
        ),
        (
            "worker hang in child initializer -> bounded kill -> serial fallback",
            test_worker_hang_in_initializer_falls_back,
        ),
        (
            "worker traps SIGTERM + hangs mid-task -> bounded SIGKILL reap -> "
            "serial fallback",
            test_sigterm_trapping_worker_falls_back_bounded,
        ),
    ]
    failures = 0
    for name, case in cases:
        print("CASE: {0}".format(name))
        try:
            case()
        except Exception as exc:  # noqa: BLE001 - test harness
            failures += 1
            print("  FAIL: {0}: {1}".format(type(exc).__name__, exc))
        else:
            print("  PASS")
        print()
    print(_summary())
    if failures:
        print("RESULT: {0}/{1} cases FAILED".format(failures, len(cases)))
        return 1
    print("RESULT: all {0} cases PASSED".format(len(cases)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
