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
#include "capi/gensas.h"

#include "sas/Simulation.h"

#include <cstring>
#include <exception>
#include <memory>
#include <string>
#include <vector>

using namespace che::core;

struct gensas_simulation {
	std::unique_ptr<Simulation> sim;
	std::string error;
	bool hasResult = false;
};

namespace {

	void setError(char *err, size_t errLen, const std::string &msg) {
		if (err == nullptr || errLen == 0) {
			return;
		}
		size_t n = msg.size() < errLen - 1 ? msg.size() : errLen - 1;
		std::memcpy(err, msg.data(), n);
		err[n] = '\0';
	}

	const SasResult *resultOf(const gensas_simulation *h) {
		return (h != nullptr && h->hasResult) ? &h->sim->result() : nullptr;
	}

	const SasEventRecord *eventOf(const gensas_simulation *h, size_t i) {
		const SasResult *r = resultOf(h);
		return (r != nullptr && i < r->events.size()) ? &r->events[i] : nullptr;
	}

} // namespace

extern "C" {

int gensas_api_version(void) {
	return GENSAS_API_VERSION;
}

const char *gensas_version(void) {
	return "GenSAS C API 1";
}

gensas_simulation *gensas_create(const char *job_json, const char *base_dir, char *err, size_t err_len) {
	try {
		if (job_json == nullptr) {
			setError(err, err_len, "no job text given");
			return nullptr;
		}
		ModelicaJob job;
		std::string message;
		if (!parseModelicaJob(job_json, base_dir != nullptr ? base_dir : "", job, message)) {
			setError(err, err_len, message);
			return nullptr;
		}
		std::unique_ptr<gensas_simulation> h(new gensas_simulation());
		h->sim.reset(new Simulation(job));
		return h.release();
	} catch (const std::exception &e) {
		setError(err, err_len, e.what());
	} catch (...) {
		setError(err, err_len, "unknown error");
	}
	return nullptr;
}

void gensas_destroy(gensas_simulation *sim) {
	try {
		delete sim;
	} catch (...) {
	}
}

int gensas_run(gensas_simulation *sim, gensas_progress_fn progress, void *user) {
	if (sim == nullptr) {
		return GENSAS_ERROR;
	}
	try {
		sim->hasResult = false;
		sim->error.clear();
		ProgressFn fn;
		if (progress != nullptr) {
			fn = [progress, user](double t, double end) { return progress(t, end, user) != 0; };
		}
		if (sim->sim->run(fn)) {
			sim->hasResult = true;
			return GENSAS_OK;
		}
		sim->error = sim->sim->error();
		return sim->sim->cancelled() ? GENSAS_CANCELLED : GENSAS_ERROR;
	} catch (const std::exception &e) {
		sim->error = e.what();
	} catch (...) {
		sim->error = "unknown error";
	}
	return GENSAS_ERROR;
}

const char *gensas_error(const gensas_simulation *sim) {
	return sim != nullptr ? sim->error.c_str() : "";
}

size_t gensas_warning_count(const gensas_simulation *sim) {
	return sim != nullptr ? sim->sim->warnings().size() : 0;
}

const char *gensas_warning(const gensas_simulation *sim, size_t index) {
	if (sim == nullptr || index >= sim->sim->warnings().size()) {
		return nullptr;
	}
	return sim->sim->warnings()[index].c_str();
}

size_t gensas_sample_count(const gensas_simulation *sim) {
	const SasResult *r = resultOf(sim);
	return r != nullptr ? r->time.size() : 0;
}

const double *gensas_time(const gensas_simulation *sim) {
	const SasResult *r = resultOf(sim);
	return r != nullptr && !r->time.empty() ? r->time.data() : nullptr;
}

size_t gensas_variable_count(const gensas_simulation *sim) {
	const SasResult *r = resultOf(sim);
	return r != nullptr ? r->names.size() : 0;
}

const char *gensas_variable_name(const gensas_simulation *sim, size_t index) {
	const SasResult *r = resultOf(sim);
	return r != nullptr && index < r->names.size() ? r->names[index].c_str() : nullptr;
}

const double *gensas_variable_values(const gensas_simulation *sim, size_t index) {
	const SasResult *r = resultOf(sim);
	return r != nullptr && index < r->values.size() && !r->values[index].empty() ? r->values[index].data() : nullptr;
}

size_t gensas_discrete_count(const gensas_simulation *sim) {
	const SasResult *r = resultOf(sim);
	return r != nullptr ? r->discreteNames.size() : 0;
}

const char *gensas_discrete_name(const gensas_simulation *sim, size_t index) {
	const SasResult *r = resultOf(sim);
	return r != nullptr && index < r->discreteNames.size() ? r->discreteNames[index].c_str() : nullptr;
}

const double *gensas_discrete_values(const gensas_simulation *sim, size_t index) {
	const SasResult *r = resultOf(sim);
	return r != nullptr && index < r->discrete.size() && !r->discrete[index].empty() ? r->discrete[index].data() : nullptr;
}

int gensas_terminated(const gensas_simulation *sim) {
	const SasResult *r = resultOf(sim);
	return r != nullptr && r->terminated ? 1 : 0;
}

double gensas_termination_time(const gensas_simulation *sim) {
	const SasResult *r = resultOf(sim);
	return r != nullptr && r->terminated ? r->endTime : 0.0;
}

const char *gensas_termination_message(const gensas_simulation *sim) {
	const SasResult *r = resultOf(sim);
	return r != nullptr && r->terminated ? r->terminateMessage.c_str() : nullptr;
}

size_t gensas_event_count(const gensas_simulation *sim) {
	const SasResult *r = resultOf(sim);
	return r != nullptr ? r->events.size() : 0;
}

double gensas_event_time(const gensas_simulation *sim, size_t event) {
	const SasEventRecord *e = eventOf(sim, event);
	return e != nullptr ? e->time : 0.0;
}

size_t gensas_event_assignment_count(const gensas_simulation *sim, size_t event) {
	const SasEventRecord *e = eventOf(sim, event);
	return e != nullptr ? e->assignments.size() : 0;
}

const char *gensas_event_assignment(const gensas_simulation *sim, size_t event, size_t index, double *value) {
	const SasEventRecord *e = eventOf(sim, event);
	if (e == nullptr || index >= e->assignments.size()) {
		return nullptr;
	}
	if (value != nullptr) {
		*value = e->assignments[index].second;
	}
	return e->assignments[index].first.c_str();
}

size_t gensas_event_reinit_count(const gensas_simulation *sim, size_t event) {
	const SasEventRecord *e = eventOf(sim, event);
	return e != nullptr ? e->reinits.size() : 0;
}

const char *gensas_event_reinit(const gensas_simulation *sim, size_t event, size_t index, double *value) {
	const SasEventRecord *e = eventOf(sim, event);
	if (e == nullptr || index >= e->reinits.size()) {
		return nullptr;
	}
	if (value != nullptr) {
		*value = e->reinits[index].second;
	}
	return e->reinits[index].first.c_str();
}

int gensas_write_result(const gensas_simulation *sim, const char *format, const char *path, char *err, size_t err_len) {
	try {
		if (sim == nullptr || !sim->hasResult) {
			setError(err, err_len, "there is no result to write");
			return 1;
		}
		if (format == nullptr || path == nullptr) {
			setError(err, err_len, "format and path are required");
			return 1;
		}
		std::string message;
		if (!sim->sim->writeResult(format, path, message)) {
			setError(err, err_len, message);
			return 1;
		}
		return 0;
	} catch (const std::exception &e) {
		setError(err, err_len, e.what());
	} catch (...) {
		setError(err, err_len, "unknown error");
	}
	return 1;
}

} // extern "C"
