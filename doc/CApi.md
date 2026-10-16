# C API

`capi/gensas.h` is a small C interface to ModelicaSAS for programs written in C, Python (`ctypes`), Julia, Rust and so on. It wraps the `Simulation` class (see [Architecture](Architecture.md)); power flow is not part of it.

A simulation is described by a job (JSON text, see [Job format](JobFormat.md)) and lives behind an opaque handle:

```c
#include "capi/gensas.h"

char err[256];
gensas_simulation *sim = gensas_create(job_json, "/path/of/job/dir", err, sizeof err);
if (sim == NULL) { /* err has the message */ }
int rc = gensas_run(sim, progress, user);      /* GENSAS_OK, GENSAS_ERROR or GENSAS_CANCELLED */
if (rc == GENSAS_OK) {
    size_t n = gensas_sample_count(sim);
    const double *t = gensas_time(sim);
    const double *x = gensas_variable_values(sim, 0);   /* gensas_variable_name(sim, 0) */
}
gensas_destroy(sim);
```

* **Progress and cancellation**: `progress(time, end_time, user)` is called after each solution segment; return non-zero to cancel (`gensas_run` then returns `GENSAS_CANCELLED`).
* **Results** (`gensas_time`, `gensas_variable_*`, `gensas_discrete_*`, `gensas_event_*`) belong to the handle and stay valid until `gensas_destroy` or the next `gensas_run` on it. Arrays are aligned with the time array; a hybrid model has two samples at each event time (before and after). `gensas_terminated`, `gensas_termination_time` and `gensas_termination_message` tell whether a `terminate()` in the model ended the run early (the run still succeeds; the results stop at that time). Failed `assert()` calls at warning level appear in `gensas_warning` after the run. Result files can be written with `gensas_write_result(sim, "json" | "csv" | "mat", path, ...)`.
* **Checking without running**: `gensas_check(sim)` parses and checks the model only (for editors). The errors and warnings about the model text are then available as diagnostics: `gensas_diagnostic_count`, and per index `gensas_diagnostic_severity` (`GENSAS_DIAG_ERROR` or `GENSAS_DIAG_WARNING`), `gensas_diagnostic_line`, `gensas_diagnostic_column` (1-based, 0 if unknown), `gensas_diagnostic_message` and `gensas_diagnostic_source` (`parser` or `checker`). `gensas_run` fills them too. Errors of the solver itself (a function leaving its domain, a stalled step) are not tied to a line and are only in `gensas_error`. A model that is not hybrid is also checked for a matching number of equations and unknowns.
* **Model info**: after `gensas_check`, `gensas_model_info_json(sim)` returns JSON text with the model `name`, `description`, `hybrid`, the `variables` (`name`, `kind` "state", "algebraic" or "discrete", `type`, `unit`, `description`, `start` or null, `input`, `line`) and the `parameters` (`name`, `type`, `unit`, `description`, `value`, `constant`, `line`). It is NULL until a check or run has read the model, and when the text could not be parsed.
* **Errors**: no function throws. `gensas_create` and `gensas_write_result` write a message into the buffer you pass (truncated and terminated if too short, or pass `NULL`/0); `gensas_run` keeps its message for `gensas_error`. Options moved into the supported range are listed by `gensas_warning`.
* **Threads**: handles are independent, so different handles may run on different threads at the same time (`capi/test/capi_test.cpp` does this). One handle must not be used from several threads at once. Solver messages are not printed unless the job asks for `verbose` output, which is meant for single runs.

Build: `bazel build //capi:libgensas.so` (Linux, macOS) or `//capi:gensas.dll` (Windows) produces the shared library (the C test is `//capi/test:capi_c_test`). `GENSAS_API_VERSION` and `gensas_api_version()` identify the interface version (2 added `gensas_check`, the diagnostics and `gensas_model_info_json`); a change that breaks existing callers increases it.
