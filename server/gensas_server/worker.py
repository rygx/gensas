"""Code that runs in the worker processes.

Every check and every run happens in its own process, so a crash of the solver (or a run that never
ends) cannot take the server down: the server only reads the messages the worker sends over a pipe and
kills the process when it must. Messages are tuples whose first item names the kind:

    ("progress", time, end)
    ("completed", info)       info: summary, warnings, diagnostics, termination
    ("failed", message, diagnostics, warnings)
    ("checked", response)     response of POST /check
"""

import json
import os
import time as _time

PROGRESS_INTERVAL = 0.1  # seconds between progress messages


def _diagnostics(items):
    return [d.as_dict() if hasattr(d, "as_dict") else d for d in items]


def run_worker(job, out_dir, conn):
    """Runs one job and stores the result in `out_dir` (result.npz, meta.json, result.csv, result.mat)."""
    import gensas
    import numpy as np

    if os.environ.get("GENSAS_SERVER_TEST_CRASH"):  # lets the tests check that a crashing worker is survived
        os.abort()
    sim = None
    try:
        sim = gensas.Simulation(job)
        last = [0.0]

        def progress(t, end):
            now = _time.monotonic()
            if now - last[0] >= PROGRESS_INTERVAL:
                last[0] = now
                conn.send(("progress", t, end))
            return False

        try:
            result = sim.run(progress)
        except gensas.GenSASError as err:
            conn.send(("failed", str(err), _diagnostics(getattr(err, "diagnostics", ())), list(sim.warnings)))
            return
        names = list(result.variables)
        discrete = list(result.discrete)
        arrays = {"time": np.asarray(result.time, dtype=float)}
        for i, name in enumerate(names):
            arrays["v%d" % i] = np.asarray(result.variables[name], dtype=float)
        for i, name in enumerate(discrete):
            arrays["d%d" % i] = np.asarray(result.discrete[name], dtype=float)
        np.savez(os.path.join(out_dir, "result.npz"), **arrays)
        meta = {"variables": names, "discrete": discrete, "events": result.events, "termination": result.termination}
        with open(os.path.join(out_dir, "meta.json"), "w", encoding="utf-8") as f:
            json.dump(meta, f)
        sim.write("csv", os.path.join(out_dir, "result.csv"))
        sim.write("mat", os.path.join(out_dir, "result.mat"))
        time_axis = arrays["time"]
        conn.send(
            (
                "completed",
                {
                    "summary": {
                        "samples": int(len(time_axis)),
                        "variables": names,
                        "discrete": discrete,
                        "events": len(result.events),
                        "endTime": float(time_axis[-1]) if len(time_axis) else 0.0,
                    },
                    "warnings": list(sim.warnings),
                    "diagnostics": _diagnostics(sim.diagnostics),
                    "termination": result.termination,
                },
            )
        )
    except gensas.GenSASError as err:  # the job itself was refused
        conn.send(("failed", str(err), _diagnostics(getattr(err, "diagnostics", ())), []))
    finally:
        if sim is not None:
            sim.close()
        conn.close()


def check_worker(job, conn):
    """Checks a job's model without running it and sends the response of POST /check."""
    import gensas

    sim = None
    try:
        try:
            sim = gensas.Simulation(job)
            found = sim.check()
        except gensas.GenSASError as err:
            found = list(getattr(err, "diagnostics", ())) or [gensas.Diagnostic("error", 0, 0, str(err), "job")]
        diagnostics = _diagnostics(found)
        info = sim.model_info() if sim is not None else None
        conn.send(("checked", {"ok": not any(d["severity"] == "error" for d in diagnostics), "diagnostics": diagnostics, "info": info}))
    finally:
        if sim is not None:
            sim.close()
        conn.close()
