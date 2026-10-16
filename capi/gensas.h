//
// Copyright (C) 2022, UChicago Argonne, LLC. All rights reserved.
//
// Software Name: Generic Semi-Analytical Simulation Tool (GenSAS)
// By: Argonne National Laboratory
// OPEN SOURCE LICENSE
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
//
// ******************************************************************************************************
// DISCLAIMER
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// ***************************************************************************************************
//
#ifndef GENSAS_H
#define GENSAS_H

/*
 * GenSAS C API, version 1: run Modelica simulations from C, Python (ctypes), Julia, Rust, ...
 *
 * A simulation is described by a job (JSON text, see doc/JobFormat.md) and lives behind an opaque
 * handle. Handles are independent: different handles may be used from different threads at the same
 * time, one handle must not be used from several threads at once. No function throws; failures are
 * reported through return values and gensas_error() / the err buffers.
 *
 * Typical use:
 *     gensas_simulation *sim = gensas_create(job_json, base_dir, err, sizeof err);
 *     if (sim && gensas_run(sim, progress, user) == GENSAS_OK) {
 *         size_t n = gensas_sample_count(sim);
 *         const double *t = gensas_time(sim);
 *         const double *x = gensas_variable_values(sim, 0);
 *     }
 *     gensas_destroy(sim);
 */

#include <stddef.h>

#if defined(_WIN32)
#if defined(GENSAS_BUILD_SHARED)
#define GENSAS_API __declspec(dllexport)
#else
#define GENSAS_API
#endif
#elif defined(__GNUC__)
#define GENSAS_API __attribute__((visibility("default")))
#else
#define GENSAS_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Version 2 adds gensas_check, gensas_diagnostic_* and gensas_model_info_json. */
#define GENSAS_API_VERSION 2

/* Return codes of gensas_run. */
#define GENSAS_OK 0
#define GENSAS_ERROR 1
#define GENSAS_CANCELLED 2

typedef struct gensas_simulation gensas_simulation;

/* Called after each solution segment with the time reached and the end time. Return non-zero to cancel. */
typedef int (*gensas_progress_fn)(double time, double end_time, void *user);

/* Version of this API (GENSAS_API_VERSION at build time) and the library version text. */
GENSAS_API int gensas_api_version(void);
GENSAS_API const char *gensas_version(void);

/* Creates a simulation from the text of a job file. Relative paths in the job resolve against
 * base_dir (NULL or "" for the current directory). Returns NULL on error and writes the message to
 * err (if err_len > 0). Warnings about adjusted options are available through gensas_warning. */
GENSAS_API gensas_simulation *gensas_create(const char *job_json, const char *base_dir, char *err, size_t err_len);

/* Frees a simulation and its results. NULL is allowed. */
GENSAS_API void gensas_destroy(gensas_simulation *sim);

/* Runs the simulation (progress may be NULL). Returns GENSAS_OK, GENSAS_ERROR (see gensas_error) or
 * GENSAS_CANCELLED. Running again replaces the previous result. */
GENSAS_API int gensas_run(gensas_simulation *sim, gensas_progress_fn progress, void *user);

/* Message of the last failed gensas_run; empty if none. Valid until the next call on this handle. */
GENSAS_API const char *gensas_error(const gensas_simulation *sim);

/* Parses and checks the model without running it, for editors. Returns GENSAS_OK or GENSAS_ERROR (see
 * gensas_error); the errors and warnings about the text are available through gensas_diagnostic_*. */
GENSAS_API int gensas_check(gensas_simulation *sim);

/* Errors and warnings about the model text from the last gensas_check or gensas_run: severity
 * (GENSAS_DIAG_ERROR or GENSAS_DIAG_WARNING), 1-based line and column (0 if unknown), the message, and
 * the stage that found it ("parser" or "checker"). Errors of the solver itself (a function leaving its
 * domain, ...) are not tied to a line and only available through gensas_error. Index out of range gives
 * 0 or NULL. Strings stay valid until the next call that runs or checks the same handle. */
/* What gensas_check found out about the model, as JSON text: name, description, hybrid, the variables
 * (name, kind "state"/"algebraic"/"discrete", type, unit, description, start or null, input, line) and the
 * parameters (name, type, unit, description, value, constant, line). NULL if the model could not be read.
 * Valid until the next call that runs or checks the same handle. */
GENSAS_API const char *gensas_model_info_json(const gensas_simulation *sim);

#define GENSAS_DIAG_ERROR 0
#define GENSAS_DIAG_WARNING 1
GENSAS_API size_t gensas_diagnostic_count(const gensas_simulation *sim);
GENSAS_API int gensas_diagnostic_severity(const gensas_simulation *sim, size_t index);
GENSAS_API int gensas_diagnostic_line(const gensas_simulation *sim, size_t index);
GENSAS_API int gensas_diagnostic_column(const gensas_simulation *sim, size_t index);
GENSAS_API const char *gensas_diagnostic_message(const gensas_simulation *sim, size_t index);
GENSAS_API const char *gensas_diagnostic_source(const gensas_simulation *sim, size_t index);

/* Warnings produced when the job was created (options moved into the supported range, ...) and by
 * the last run (assert() calls with AssertionLevel.warning that failed). */
GENSAS_API size_t gensas_warning_count(const gensas_simulation *sim);
GENSAS_API const char *gensas_warning(const gensas_simulation *sim, size_t index);

/*
 * Results of the last successful run. Pointers stay valid until gensas_destroy or the next
 * gensas_run on the same handle. All variable arrays have gensas_sample_count entries, aligned with
 * gensas_time. Times are non-decreasing; a hybrid model has two samples at each event time (before
 * and after the event). Counts are 0 and pointers NULL when there is no result.
 */
GENSAS_API size_t gensas_sample_count(const gensas_simulation *sim);
GENSAS_API const double *gensas_time(const gensas_simulation *sim);

/* Continuous variables. */
GENSAS_API size_t gensas_variable_count(const gensas_simulation *sim);
GENSAS_API const char *gensas_variable_name(const gensas_simulation *sim, size_t index);
GENSAS_API const double *gensas_variable_values(const gensas_simulation *sim, size_t index);

/* Discrete variables (hybrid models only). */
GENSAS_API size_t gensas_discrete_count(const gensas_simulation *sim);
GENSAS_API const char *gensas_discrete_name(const gensas_simulation *sim, size_t index);
GENSAS_API const double *gensas_discrete_values(const gensas_simulation *sim, size_t index);

/* Events (hybrid models only): time, assigned discrete variables and reinitialized continuous variables. */
GENSAS_API size_t gensas_event_count(const gensas_simulation *sim);
GENSAS_API double gensas_event_time(const gensas_simulation *sim, size_t event);
GENSAS_API size_t gensas_event_assignment_count(const gensas_simulation *sim, size_t event);
GENSAS_API const char *gensas_event_assignment(const gensas_simulation *sim, size_t event, size_t index, double *value);
GENSAS_API size_t gensas_event_reinit_count(const gensas_simulation *sim, size_t event);
GENSAS_API const char *gensas_event_reinit(const gensas_simulation *sim, size_t event, size_t index, double *value);

/*
 * Hybrid models: 1 if a terminate() in a when-equation ended the run early (the result then ends at
 * gensas_termination_time and the run still counts as successful), else 0. The message is the text
 * given to terminate(), NULL if the run was not terminated.
 */
GENSAS_API int gensas_terminated(const gensas_simulation *sim);
GENSAS_API double gensas_termination_time(const gensas_simulation *sim);
GENSAS_API const char *gensas_termination_message(const gensas_simulation *sim);

/* Writes the result as "json", "csv" or "mat". Returns 0 on success; otherwise non-zero with a message in err. */
GENSAS_API int gensas_write_result(const gensas_simulation *sim, const char *format, const char *path, char *err, size_t err_len);

#ifdef __cplusplus
}
#endif

#endif
