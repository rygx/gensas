// Modelica language features end to end: each model is run through the Simulation facade and
// compared with its closed-form solution.
#include "sas/Simulation.h"
#include "sas/test/test_util.h"

#include <cmath>
#include <string>

using namespace che::core;
using std::string;

static double lastValue(const SasResult &r, const string &name) {
	for (size_t i = 0; i < r.names.size(); i++) {
		if (r.names[i] == name) {
			return r.values[i].back();
		}
	}
	return NAN;
}

static ModelicaJob job(const char *src, double endTime) {
	ModelicaJob j;
	j.modelSource = src;
	j.options.endTime = endTime;
	return j;
}

// An equation side that is only a function call: der(x) = sqrt(y).
// The temporary standing for sqrt(y) must not take its initial value from der(x).
static void testBareFunctionCall() {
	const char *src =
		"model M\n"
		"  Real x(start = 1);\n"
		"  Real y(start = 1);\n"
		"equation\n"
		"  der(x) = sqrt(y);\n"
		"  der(y) = -0.1*y;\n"
		"end M;\n";
	Simulation sim(job(src, 2.0));
	CHECK(sim.run());
	CHECK(sim.error().empty());
	CHECK_NEAR(lastValue(sim.result(), "y"), std::exp(-0.2), 1e-4);
	CHECK_NEAR(lastValue(sim.result(), "x"), 1.0 + 20.0 * (1.0 - std::exp(-0.1)), 1e-4);
}

int main() {
	testBareFunctionCall();
	return sastest::finish("modelica_language_test");
}
