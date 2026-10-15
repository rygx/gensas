"""Python interface to GenSAS (ModelicaSAS) through the C API in capi/gensas.h.

    import gensas
    result = gensas.run({
        "gensas": {"format": "job", "version": 1}, "kind": "modelica",
        "model": {"source": "model D Real x(start=1); equation der(x) = -x; end D;"},
        "solver": {"endTime": 2},
    })
    result.time, result.variables["x"]

A job is the dictionary form of a job file (doc/JobFormat.md), a JSON string, or the path of a job
file. The shared library is found through the GENSAS_LIBRARY environment variable, next to this
package, or in bazel-bin/capi of a source checkout (build it with `bazel build //capi:libgensas.so`,
or //capi:gensas.dll on Windows). Several Simulation objects may run on different threads at once.
Only the standard library is required; results are NumPy arrays when NumPy is installed and lists
of floats otherwise.
"""

import ctypes
import json
import os
import sys

try:
    import numpy as _np
except ImportError:  # NumPy is optional
    _np = None

__all__ = ["Simulation", "Result", "run", "GenSASError", "Cancelled", "library_path"]

_OK, _ERROR, _CANCELLED = 0, 1, 2
_ERR_LEN = 1024


class GenSASError(Exception):
    """A job, model or solver error reported by the library."""


class Cancelled(GenSASError):
    """The progress callback asked to stop the run."""


def _candidate_paths():
    env = os.environ.get("GENSAS_LIBRARY")
    if env:
        yield env
        return
    names = {"win32": ["gensas.dll"], "darwin": ["libgensas.dylib", "libgensas.so"]}.get(sys.platform, ["libgensas.so"])
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(os.path.dirname(here))
    # next to the package, in a Bazel runfiles tree (capi/ below the working directory), or in bazel-bin
    for folder in (here, os.path.join(os.getcwd(), "capi"), os.path.join(root, "bazel-bin", "capi")):
        for name in names:
            yield os.path.join(folder, name)


def library_path():
    """Path of the shared library that will be (or was) loaded."""
    for path in _candidate_paths():
        if os.path.exists(path):
            return path
    raise GenSASError("cannot find the GenSAS shared library; build //capi:libgensas.so or set GENSAS_LIBRARY")


_lib = None


def _load():
    global _lib
    if _lib is not None:
        return _lib
    path = library_path()
    if hasattr(os, "add_dll_directory"):  # Windows: let the DLL find its neighbours
        os.add_dll_directory(os.path.dirname(os.path.abspath(path)))
    lib = ctypes.CDLL(path)
    c_dp = ctypes.POINTER(ctypes.c_double)
    progress = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_void_p)
    sim = ctypes.c_void_p
    sig = {
        "gensas_api_version": (ctypes.c_int, []),
        "gensas_version": (ctypes.c_char_p, []),
        "gensas_create": (sim, [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_size_t]),
        "gensas_destroy": (None, [sim]),
        "gensas_run": (ctypes.c_int, [sim, progress, ctypes.c_void_p]),
        "gensas_error": (ctypes.c_char_p, [sim]),
        "gensas_warning_count": (ctypes.c_size_t, [sim]),
        "gensas_warning": (ctypes.c_char_p, [sim, ctypes.c_size_t]),
        "gensas_sample_count": (ctypes.c_size_t, [sim]),
        "gensas_time": (c_dp, [sim]),
        "gensas_variable_count": (ctypes.c_size_t, [sim]),
        "gensas_variable_name": (ctypes.c_char_p, [sim, ctypes.c_size_t]),
        "gensas_variable_values": (c_dp, [sim, ctypes.c_size_t]),
        "gensas_discrete_count": (ctypes.c_size_t, [sim]),
        "gensas_discrete_name": (ctypes.c_char_p, [sim, ctypes.c_size_t]),
        "gensas_discrete_values": (c_dp, [sim, ctypes.c_size_t]),
        "gensas_event_count": (ctypes.c_size_t, [sim]),
        "gensas_event_time": (ctypes.c_double, [sim, ctypes.c_size_t]),
        "gensas_event_assignment_count": (ctypes.c_size_t, [sim, ctypes.c_size_t]),
        "gensas_event_assignment": (ctypes.c_char_p, [sim, ctypes.c_size_t, ctypes.c_size_t, c_dp]),
        "gensas_event_reinit_count": (ctypes.c_size_t, [sim, ctypes.c_size_t]),
        "gensas_event_reinit": (ctypes.c_char_p, [sim, ctypes.c_size_t, ctypes.c_size_t, c_dp]),
        "gensas_write_result": (ctypes.c_int, [sim, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_size_t]),
    }
    for name, (restype, argtypes) in sig.items():
        fn = getattr(lib, name)
        fn.restype = restype
        fn.argtypes = argtypes
    lib._progress_type = progress
    if lib.gensas_api_version() != 1:
        raise GenSASError("the GenSAS library has API version %d, this module expects 1" % lib.gensas_api_version())
    _lib = lib
    return lib


def _array(pointer, n):
    if not pointer or n == 0:
        return _np.zeros(0) if _np is not None else []
    if _np is not None:
        return _np.ctypeslib.as_array(pointer, shape=(n,)).copy()  # the handle owns the C memory
    return pointer[:n]


def _text(raw):
    return raw.decode("utf-8", "replace") if raw is not None else None


class Result:
    """Samples of one run. `time` and every array in `variables` / `discrete` have the same length;
    a hybrid model has two samples at each event time (before and after the event)."""

    def __init__(self, time, variables, discrete, events):
        self.time = time
        self.variables = variables
        self.discrete = discrete
        self.events = events  # list of {"time", "assignments": {name: value}, "reinits": {name: value}}

    def __repr__(self):
        return "Result(%d samples, variables=%s, %d events)" % (len(self.time), list(self.variables), len(self.events))


class Simulation:
    """One simulation. `job` is a dict, a JSON string or the path of a job file; relative paths in
    a dict or string job resolve against `base_dir` (default: the current directory)."""

    def __init__(self, job, base_dir=None):
        lib = _load()
        if base_dir is not None:
            base_dir = os.fspath(base_dir)
        if isinstance(job, dict):
            text = json.dumps(job)
        elif isinstance(job, (str, os.PathLike)) and not str(job).lstrip().startswith("{"):
            path = os.fspath(job)
            with open(path, "r", encoding="utf-8") as f:
                text = f.read()
            if base_dir is None:
                base_dir = os.path.dirname(os.path.abspath(path))
        else:
            text = str(job)
        err = ctypes.create_string_buffer(_ERR_LEN)
        handle = lib.gensas_create(text.encode("utf-8"), base_dir.encode("utf-8") if base_dir else None, err, _ERR_LEN)
        if not handle:
            raise GenSASError(err.value.decode("utf-8", "replace"))
        self._lib = lib
        self._handle = ctypes.c_void_p(handle)
        self.warnings = [_text(lib.gensas_warning(self._handle, i)) for i in range(lib.gensas_warning_count(self._handle))]

    def run(self, progress=None):
        """Runs the simulation and returns a Result. `progress(time, end_time)` is called after each
        solution segment; return True from it to cancel (raises Cancelled)."""
        lib = self._lib
        self._check_open()
        failure = []
        callback = lib._progress_type()  # a NULL function pointer
        if progress is not None:
            def callback_impl(t, end, _user):
                try:
                    return 1 if progress(t, end) else 0
                except BaseException as exc:  # never let an exception cross the C boundary
                    failure.append(exc)
                    return 1
            callback = lib._progress_type(callback_impl)
        status = lib.gensas_run(self._handle, callback, None)
        if failure:
            raise failure[0]
        if status == _CANCELLED:
            raise Cancelled("the run was cancelled")
        if status != _OK:
            raise GenSASError(_text(lib.gensas_error(self._handle)))
        return self._result()

    def _result(self):
        lib, h = self._lib, self._handle
        n = lib.gensas_sample_count(h)
        time = _array(lib.gensas_time(h), n)
        variables = {_text(lib.gensas_variable_name(h, i)): _array(lib.gensas_variable_values(h, i), n) for i in range(lib.gensas_variable_count(h))}
        discrete = {_text(lib.gensas_discrete_name(h, i)): _array(lib.gensas_discrete_values(h, i), n) for i in range(lib.gensas_discrete_count(h))}
        events = []
        value = ctypes.c_double()
        for e in range(lib.gensas_event_count(h)):
            assignments, reinits = {}, {}
            for k in range(lib.gensas_event_assignment_count(h, e)):
                name = lib.gensas_event_assignment(h, e, k, ctypes.byref(value))
                assignments[_text(name)] = value.value
            for k in range(lib.gensas_event_reinit_count(h, e)):
                name = lib.gensas_event_reinit(h, e, k, ctypes.byref(value))
                reinits[_text(name)] = value.value
            events.append({"time": lib.gensas_event_time(h, e), "assignments": assignments, "reinits": reinits})
        return Result(time, variables, discrete, events)

    def write(self, format, path):
        """Writes the result of the last successful run as "json", "csv" or "mat"."""
        self._check_open()
        err = ctypes.create_string_buffer(_ERR_LEN)
        if self._lib.gensas_write_result(self._handle, format.encode("utf-8"), os.fspath(path).encode("utf-8"), err, _ERR_LEN) != 0:
            raise GenSASError(err.value.decode("utf-8", "replace"))

    def _check_open(self):
        if not getattr(self, "_handle", None):
            raise GenSASError("the simulation is closed")

    def close(self):
        if getattr(self, "_handle", None):
            self._lib.gensas_destroy(self._handle)
            self._handle = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass


def run(job, base_dir=None, progress=None):
    """Creates a Simulation, runs it and returns the Result."""
    with Simulation(job, base_dir) as sim:
        return sim.run(progress)
