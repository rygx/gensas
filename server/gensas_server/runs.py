"""Runs and checks: state, worker processes, limits and clean-up."""

import multiprocessing
import os
import shutil
import tempfile
import threading
import time
import uuid
from collections import OrderedDict

from . import worker

QUEUED, RUNNING, COMPLETED, TERMINATED, FAILED, CANCELLED = "queued", "running", "completed", "terminated", "failed", "cancelled"
ACTIVE = (QUEUED, RUNNING)
HAS_RESULT = (COMPLETED, TERMINATED)


class Busy(Exception):
    """No capacity: too many runs, or no worker free for a check in time."""


class Run:
    def __init__(self, run_id, job, directory):
        self.id = run_id
        self.job = job
        self.directory = directory
        self.state = QUEUED
        self.created = time.time()
        self.started = None
        self.finished = None
        self.progress_time = 0.0
        self.progress_end = 0.0
        self.error = None
        self.diagnostics = []
        self.warnings = []
        self.termination = None
        self.summary = None
        self.cancel_requested = False
        self.process = None

    def to_dict(self):
        elapsed = None
        if self.started is not None:
            elapsed = (self.finished if self.finished is not None else time.time()) - self.started
        return {
            "id": self.id,
            "kind": "modelica",
            "state": self.state,
            "createdAt": self.created,
            "startedAt": self.started,
            "finishedAt": self.finished,
            "progress": {"time": self.progress_time, "end": self.progress_end},
            "error": self.error,
            "diagnostics": list(self.diagnostics),
            "warnings": list(self.warnings),
            "elapsedSeconds": elapsed,
            "termination": self.termination,
            "summary": self.summary,
        }


class RunManager:
    def __init__(self, config):
        self.config = config
        self._ctx = multiprocessing.get_context("spawn")  # safe with threads and with the shared library
        self._lock = threading.RLock()
        self._runs = OrderedDict()
        self._slots = threading.Semaphore(config.max_concurrent_runs)
        self._check_slots = threading.Semaphore(config.max_concurrent_checks)
        if config.work_dir:
            os.makedirs(config.work_dir, exist_ok=True)
            self._root = tempfile.mkdtemp(prefix="gensas-", dir=config.work_dir)
        else:
            self._root = tempfile.mkdtemp(prefix="gensas-")
        self._closed = False

    # ---- runs -------------------------------------------------------------------------------

    def create(self, job):
        with self._lock:
            self._collect()
            active = sum(1 for r in self._runs.values() if r.state in ACTIVE)
            if active >= self.config.max_concurrent_runs + self.config.max_queued_runs:
                raise Busy("too many runs are waiting or running; try again later")
            run_id = uuid.uuid4().hex[:12]
            directory = os.path.join(self._root, run_id)
            os.makedirs(directory)
            run = Run(run_id, job, directory)
            self._runs[run_id] = run
        threading.Thread(target=self._drive, args=(run,), name="run-" + run_id, daemon=True).start()
        return run

    def get(self, run_id):
        with self._lock:
            return self._runs.get(run_id)

    def snapshot(self, run):
        with self._lock:
            return run.to_dict()

    def list(self):
        with self._lock:
            return [r.to_dict() for r in reversed(self._runs.values())]

    def cancel(self, run):
        """Stops a queued or running run. Returns False if it had already finished."""
        with self._lock:
            if run.state not in ACTIVE:
                return False
            run.cancel_requested = True
            process = run.process
        if process is not None and process.is_alive():
            process.kill()
        return True

    def delete(self, run):
        with self._lock:
            self._runs.pop(run.id, None)
        shutil.rmtree(run.directory, ignore_errors=True)

    def _collect(self):
        """Removes finished runs that are too old or too many."""
        now = time.time()
        finished = [r for r in self._runs.values() if r.state not in ACTIVE]
        drop = [r for r in finished if r.finished is not None and now - r.finished > self.config.retention]
        rest = [r for r in finished if r not in drop]
        drop += rest[: max(0, len(rest) - self.config.max_kept_runs)]
        for r in drop:
            self._runs.pop(r.id, None)
            shutil.rmtree(r.directory, ignore_errors=True)

    def _finish(self, run, state, error=None):
        with self._lock:
            run.state = state
            run.finished = time.time()
            if error is not None:
                run.error = error

    def _drive(self, run):
        """Thread of one run: waits for a worker slot, starts the process and follows its messages."""
        while not self._slots.acquire(timeout=0.1):
            if run.cancel_requested or self._closed:
                self._finish(run, CANCELLED, "cancelled before it started")
                return
        try:
            if run.cancel_requested or self._closed:
                self._finish(run, CANCELLED, "cancelled before it started")
                return
            self._execute(run)
        finally:
            self._slots.release()

    def _execute(self, run):
        parent, child = self._ctx.Pipe(duplex=False)
        process = self._ctx.Process(target=worker.run_worker, args=(run.job, run.directory, child), daemon=True)
        with self._lock:
            run.state = RUNNING
            run.started = time.time()
            run.process = process
        process.start()
        child.close()
        deadline = run.started + self.config.run_timeout
        outcome = None
        try:
            while outcome is None:
                if parent.poll(0.05):
                    try:
                        message = parent.recv()
                    except (EOFError, OSError):
                        break
                    kind = message[0]
                    if kind == "progress":
                        with self._lock:
                            run.progress_time, run.progress_end = message[1], message[2]
                    else:
                        outcome = message
                elif not process.is_alive():
                    # the process ended; take a message that arrived just before
                    if parent.poll(0.05):
                        continue
                    break
                if outcome is None and time.time() > deadline:
                    process.kill()
                    process.join()
                    self._finish(run, FAILED, "the run exceeded the time limit of %g seconds" % self.config.run_timeout)
                    return
        finally:
            if process.is_alive():
                process.kill()
            process.join()
            parent.close()
        if run.cancel_requested:
            self._finish(run, CANCELLED, "cancelled")
        elif outcome is None:
            code = process.exitcode
            why = "the solver process crashed (signal %d)" % -code if code is not None and code < 0 else "the solver process ended without a result (exit code %s)" % code
            self._finish(run, FAILED, why)
        elif outcome[0] == "completed":
            info = outcome[1]
            with self._lock:
                run.summary = info["summary"]
                run.warnings = info["warnings"]
                run.diagnostics = info["diagnostics"]
                run.termination = info["termination"]
                run.progress_time = run.summary["endTime"]
                if run.progress_end <= 0.0:
                    run.progress_end = run.summary["endTime"]
            self._finish(run, TERMINATED if info["termination"] else COMPLETED)
        else:
            with self._lock:
                run.diagnostics = outcome[2]
                run.warnings = outcome[3]
            self._finish(run, FAILED, outcome[1])

    # ---- checks -----------------------------------------------------------------------------

    def check(self, job):
        """Checks a model in a worker process. Blocks up to the check timeout; raises Busy if no worker
        is free in time or the check does not finish."""
        if not self._check_slots.acquire(timeout=self.config.check_timeout):
            raise Busy("the server is busy checking other models; try again")
        try:
            parent, child = self._ctx.Pipe(duplex=False)
            process = self._ctx.Process(target=worker.check_worker, args=(job, child), daemon=True)
            process.start()
            child.close()
            try:
                if parent.poll(self.config.check_timeout):
                    message = parent.recv()
                    if message[0] == "checked":
                        return message[1]
                    raise Busy("the check ended without a result")
                raise Busy("the check took longer than %g seconds" % self.config.check_timeout)
            except (EOFError, OSError):
                crash = {"severity": "error", "line": 0, "column": 0, "source": "server",
                         "message": "the checker crashed on this model; please report it with the model text"}
                return {"ok": False, "diagnostics": [crash], "info": None}
            finally:
                if process.is_alive():
                    process.kill()
                process.join()
                parent.close()
        finally:
            self._check_slots.release()

    # ---- shutdown ---------------------------------------------------------------------------

    def close(self):
        self._closed = True
        with self._lock:
            runs = list(self._runs.values())
        for run in runs:
            self.cancel(run)
        for run in runs:
            deadline = time.time() + 2.0
            while run.state in ACTIVE and time.time() < deadline:
                time.sleep(0.02)
        shutil.rmtree(self._root, ignore_errors=True)
