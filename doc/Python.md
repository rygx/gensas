# Python interface

`python/gensas` is a small `ctypes` wrapper around the [C API](CApi.md). It needs only the Python standard library (results are NumPy arrays if NumPy is installed, lists of floats otherwise) and the shared library built from `//capi`.

```python
import sys
sys.path.insert(0, "python")          # or copy python/gensas next to your script
import gensas

result = gensas.run({
    "gensas": {"format": "job", "version": 1}, "kind": "modelica",
    "model": {"file": "resources/mofile/test_bouncing_ball.mo"},
    "solver": {"endTime": 1.2},
})
result.time                 # sample times
result.variables["h"]       # continuous variables by name
result.discrete             # discrete variables (hybrid models)
result.events               # [{"time": ..., "assignments": {...}, "reinits": {"v": ...}}, ...]
result.termination          # None, or {"time": ..., "message": ...} if terminate() ended the run early
```

* A job is the dictionary form of a [job file](JobFormat.md), a JSON string, or the path of a job file (then relative paths resolve against that file; for a dictionary or string use `base_dir=`).
* `gensas.Simulation(job)` gives control: `run(progress)` returns the result, `write("json" | "csv" | "mat", path)` writes the last result, `warnings` lists adjusted options, `close()` (or a `with` block) frees the handle.
* `progress(time, end_time)` is called after each solution segment; return `True` to cancel (`gensas.Cancelled` is raised). An exception raised in the callback stops the run and is re-raised.
* `sim.check()` parses and checks the model without running it and returns a list of `Diagnostic` objects (`severity` "error" or "warning", `line`, `column`, `message`, `source`; `as_dict()` for JSON); an empty list means the model is fine. `sim.diagnostics` gives the same for the last `run()` or `check()`, and a `GenSASError` raised for a bad model carries them in its `diagnostics` attribute.
* `sim.model_info()` returns what the last `check()` or `run()` found out about the model as a dict (`name`, `description`, `hybrid`, `variables` with `kind` "state"/"algebraic"/"discrete", `unit`, `description`, `start`, `input`, `line`, and `parameters` with `value`, `constant`, ...), or `None` if the model could not be read. An editor can list variables and offer parameters for changing from it.
* Errors raise `gensas.GenSASError` with the library's message (including the parser's messages for a bad model).
* Several simulations may run in different threads at the same time.

The library is found through the `GENSAS_LIBRARY` environment variable, next to the package, or in `bazel-bin/capi`. Build it with `bazel build //capi:libgensas.so` (`//capi:gensas.dll` on Windows). Tests: `bazel test //python:gensas_python_test` (Linux and macOS), or `GENSAS_LIBRARY=<path> python3 -m unittest discover -s python/tests` (this is what the Windows CI job runs against `gensas.dll`).
