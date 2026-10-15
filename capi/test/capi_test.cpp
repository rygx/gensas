// C API from C++: hybrid results with events, warnings, error reporting and several handles running
// on different threads at the same time.
#include "capi/gensas.h"
#include "sas/test/test_util.h"

#include <cmath>
#include <string>
#include <thread>
#include <vector>

static std::string jobFor(const std::string &source, double endTime) {
	std::string escaped;
	for (char c : source) {
		if (c == '\n') {
			escaped += "\\n";
		} else if (c == '"' || c == '\\') {
			escaped += '\\';
			escaped += c;
		} else {
			escaped += c;
		}
	}
	return "{\"gensas\":{\"format\":\"job\",\"version\":1},\"kind\":\"modelica\",\"model\":{\"source\":\"" + escaped + "\"},\"solver\":{\"endTime\":" + std::to_string(endTime) + "}}";
}

static const char *kTerminate =
	"model M\n"
	"  Real x(start = 1);\n"
	"equation\n"
	"  der(x) = -x;\n"
	"  when x < 0.3 then\n"
	"    terminate(\"x fell below 0.3\");\n"
	"  end when;\n"
	"end M;\n";

static const char *kBall =
	"model Ball\n"
	"  Real h(start = 1);\n"
	"  Real v(start = 0);\n"
	"equation\n"
	"  der(h) = v;\n"
	"  der(v) = -9.81;\n"
	"  when h < 0 then\n"
	"    reinit(v, -0.5*pre(v));\n"
	"  end when;\n"
	"end Ball;\n";

static const char *kGain =
	"model M\n"
	"  discrete Real a(start = 1);\n"
	"  Real x(start = 1);\n"
	"equation\n"
	"  der(x) = -a*x;\n"
	"  when time >= 1 then\n"
	"    a = 3;\n"
	"  end when;\n"
	"end M;\n";

static void testHybridResult() {
	char err[256];
	gensas_simulation *sim = gensas_create(jobFor(kGain, 2.0).c_str(), nullptr, err, sizeof err);
	CHECK(sim != nullptr);
	if (sim == nullptr) {
		return;
	}
	CHECK_EQ(gensas_run(sim, nullptr, nullptr), GENSAS_OK);
	CHECK_EQ(gensas_discrete_count(sim), (size_t)1);
	CHECK_EQ(std::string(gensas_discrete_name(sim, 0)), std::string("a"));
	size_t n = gensas_sample_count(sim);
	CHECK(gensas_discrete_values(sim, 0) != nullptr);
	if (n > 0 && gensas_discrete_values(sim, 0) != nullptr) {
		CHECK_NEAR(gensas_discrete_values(sim, 0)[n - 1], 3.0, 0.0);
		CHECK_NEAR(gensas_discrete_values(sim, 0)[0], 1.0, 0.0);
	}
	CHECK_EQ(gensas_event_count(sim), (size_t)1);
	CHECK_NEAR(gensas_event_time(sim, 0), 1.0, 1e-9);
	CHECK_EQ(gensas_event_assignment_count(sim, 0), (size_t)1);
	double value = 0.0;
	const char *name = gensas_event_assignment(sim, 0, 0, &value);
	CHECK(name != nullptr && std::string(name) == "a");
	CHECK_NEAR(value, 3.0, 0.0);
	CHECK(gensas_event_assignment(sim, 0, 1, &value) == nullptr);
	CHECK_EQ(gensas_event_reinit_count(sim, 0), (size_t)0);
	CHECK_EQ(gensas_terminated(sim), 0);
	CHECK(gensas_termination_message(sim) == nullptr);
	gensas_destroy(sim);

	// terminate() ends the run early and successfully.
	sim = gensas_create(jobFor(kTerminate, 5.0).c_str(), nullptr, err, sizeof err);
	CHECK(sim != nullptr);
	if (sim != nullptr) {
		CHECK_EQ(gensas_run(sim, nullptr, nullptr), GENSAS_OK);
		CHECK_EQ(gensas_terminated(sim), 1);
		CHECK_NEAR(gensas_termination_time(sim), std::log(1.0 / 0.3), 1e-5);
		const char *msg = gensas_termination_message(sim);
		CHECK(msg != nullptr && std::string(msg) == "x fell below 0.3");
		size_t last = gensas_sample_count(sim);
		CHECK(last > 0 && std::fabs(gensas_time(sim)[last - 1] - gensas_termination_time(sim)) < 1e-9);
		gensas_destroy(sim);
	}

	sim = gensas_create(jobFor(kBall, 1.2).c_str(), nullptr, err, sizeof err);
	CHECK(sim != nullptr);
	if (sim != nullptr) {
		CHECK_EQ(gensas_run(sim, nullptr, nullptr), GENSAS_OK);
		CHECK_EQ(gensas_event_count(sim), (size_t)3);
		CHECK_EQ(gensas_event_reinit_count(sim, 0), (size_t)1);
		const char *r = gensas_event_reinit(sim, 0, 0, &value);
		CHECK(r != nullptr && std::string(r) == "v");
		gensas_destroy(sim);
	}
}

static void testErrorsAndWarnings() {
	char err[64];
	// Truncated message is still terminated.
	CHECK(gensas_create("{\"gensas\":{\"format\":\"job\",\"version\":1},\"kind\":\"modelica\",\"model\":{\"source\":\"x\"},\"colour\":1}", nullptr, err, 12) == nullptr);
	CHECK_EQ(std::string(err).size(), (size_t)11);
	// err may be absent.
	CHECK(gensas_create("nonsense", nullptr, nullptr, 0) == nullptr);

	gensas_simulation *sim = gensas_create(jobFor("model M\n  Real x;\nequation\n  der(x) = \n", 1.0).c_str(), nullptr, err, sizeof err);
	CHECK(sim != nullptr);
	if (sim != nullptr) {
		CHECK_EQ(gensas_run(sim, nullptr, nullptr), GENSAS_ERROR);
		CHECK(std::string(gensas_error(sim)).size() > 0);
		CHECK_EQ(gensas_sample_count(sim), (size_t)0);
		CHECK_EQ(gensas_write_result(sim, "json", "unused.json", err, sizeof err), 1);
		gensas_destroy(sim);
	}

	sim = gensas_create(jobFor("model D Real x(start=1); equation der(x) = -x; end D;", 5000.0).c_str(), nullptr, err, sizeof err);
	CHECK(sim != nullptr);
	if (sim != nullptr) {
		CHECK_EQ(gensas_warning_count(sim), (size_t)1);
		CHECK(gensas_warning(sim, 0) != nullptr);
		CHECK(gensas_warning(sim, 1) == nullptr);
		gensas_destroy(sim);
	}
	CHECK_EQ(gensas_run(nullptr, nullptr, nullptr), GENSAS_ERROR);
	CHECK_EQ(gensas_sample_count(nullptr), (size_t)0);
}

// Handles are independent: many threads, each with its own simulation, get the closed-form answer.
static void testConcurrentHandles() {
	const int kThreads = 8;
	std::vector<double> out(kThreads, NAN);
	std::vector<std::thread> threads;
	for (int t = 0; t < kThreads; t++) {
		threads.emplace_back([&out, t]() {
			std::string src = "model D\n  parameter Real k = " + std::to_string(0.5 + 0.25 * t) + ";\n  Real x(start = 1);\nequation\n  der(x) = -k*x;\nend D;\n";
			char err[128];
			for (int round = 0; round < 3; round++) {
				gensas_simulation *sim = gensas_create(jobFor(src, 2.0).c_str(), nullptr, err, sizeof err);
				if (sim == nullptr) {
					return;
				}
				if (gensas_run(sim, nullptr, nullptr) == GENSAS_OK) {
					size_t n = gensas_sample_count(sim);
					out[t] = gensas_variable_values(sim, 0)[n - 1];
				}
				gensas_destroy(sim);
			}
		});
	}
	for (std::thread &th : threads) {
		th.join();
	}
	for (int t = 0; t < kThreads; t++) {
		CHECK_NEAR(out[t], std::exp(-2.0 * (0.5 + 0.25 * t)), 1e-4);
	}
}

int main() {
	testHybridResult();
	testErrorsAndWarnings();
	testConcurrentHandles();
	return sastest::finish("capi_test");
}
