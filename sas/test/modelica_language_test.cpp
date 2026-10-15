// Modelica language features end to end: each model is run through the Simulation facade and
// compared with its closed-form solution.
#include "sas/Simulation.h"
#include "sas/test/test_util.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <string>

using namespace che::core;
using std::string;

static double lastValue(const SasResult &r, const string &name) {
	for (size_t i = 0; i < r.names.size(); i++) {
		if (r.names[i] == name) {
			if (r.values[i].empty()) {
				break;
			}
			return r.values[i].back();
		}
	}
	return std::numeric_limits<double>::quiet_NaN();
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

// z = f(u) for a state u that moves linearly: u(t) = u0 + rate*t. Compares z at the end with f(u(end)).
static void checkFunction(const char *name, double u0, double rate, double endTime, double (*f)(double)) {
	char src[512];
	std::snprintf(src, sizeof src,
				  "model M\n  Real u(start = %.17g);\n  Real z;\nequation\n  der(u) = %.17g;\n  z = %s(u);\nend M;\n",
				  u0, rate, name);
	Simulation sim(job(src, endTime));
	if (!sim.run()) {
		std::cerr << name << ": " << sim.error() << std::endl;
		CHECK(false);
		return;
	}
	// The last sample can fall slightly short of the end time; compare at the time actually reached.
	double uEnd = u0 + rate * sim.result().time.back();
	double z = lastValue(sim.result(), "z");
	if (std::fabs(z - f(uEnd)) > 1e-5 * std::fmax(1.0, std::fabs(f(uEnd)))) {
		std::cerr << name << ": z = " << z << ", expected " << f(uEnd) << std::endl;
		CHECK(false);
	}
}

static double fSin(double x) { return std::sin(x); }
static double fExp(double x) { return std::exp(x); }
static double fSqrt(double x) { return std::sqrt(x); }
static double fTan(double x) { return std::tan(x); }
static double fLog(double x) { return std::log(x); }
static double fLog10(double x) { return std::log10(x); }
static double fAtan(double x) { return std::atan(x); }
static double fAsin(double x) { return std::asin(x); }
static double fAcos(double x) { return std::acos(x); }
static double fSinh(double x) { return std::sinh(x); }
static double fCosh(double x) { return std::cosh(x); }
static double fTanh(double x) { return std::tanh(x); }

// Every elementary function, as an algebraic output of a moving state.
static void testElementaryFunctions() {
	checkFunction("sin", 0.2, 0.5, 3.0, fSin);
	checkFunction("exp", 0.2, 0.5, 3.0, fExp);
	checkFunction("sqrt", 0.5, 0.5, 3.0, fSqrt);
	checkFunction("tan", 0.1, 0.3, 2.0, fTan);
	checkFunction("log", 0.5, 0.4, 2.0, fLog);
	checkFunction("log", 4.0, -1.0, 3.0, fLog);
	checkFunction("log10", 0.5, 0.4, 2.0, fLog10);
	checkFunction("atan", -2.0, 1.5, 3.0, fAtan);
	checkFunction("asin", -0.5, 0.45, 2.0, fAsin);
	checkFunction("acos", -0.5, 0.45, 2.0, fAcos);
	checkFunction("sinh", -1.0, 1.0, 3.0, fSinh);
	checkFunction("cosh", -1.0, 1.0, 3.0, fCosh);
	checkFunction("tanh", -1.0, 1.0, 3.0, fTanh);
}

// Functions inside differential equations, with closed-form solutions.
static void testFunctionsInOdes() {
	// sinh(x) = sinh(1) exp(-t)
	const char *tanhOde =
		"model M\n  Real x(start = 1);\nequation\n  der(x) = -tanh(x);\nend M;\n";
	Simulation a(job(tanhOde, 1.5));
	CHECK(a.run());
	CHECK_NEAR(lastValue(a.result(), "x"), std::asinh(std::sinh(1.0) * std::exp(-a.result().time.back())), 1e-5);

	// x = exp(c exp(t)) with c = log(x0)
	const char *logOde =
		"model M\n  Real x(start = 1.6487212707001282);\nequation\n  der(x) = x*log(x);\nend M;\n";
	Simulation b(job(logOde, 1.0));
	CHECK(b.run());
	CHECK_NEAR(lastValue(b.result(), "x"), std::exp(0.5 * std::exp(b.result().time.back())), 1e-4);

	// tan(x) = tan(x0) exp(t) is not solvable for x, but x = atan(c exp(t)) solves der(x) = sin(x)*cos(x)
	const char *atanOde =
		"model M\n  Real x(start = 0.4636476090008061);\nequation\n  der(x) = sin(x)*cos(x);\nend M;\n";
	Simulation c(job(atanOde, 1.0));
	CHECK(c.run());
	CHECK_NEAR(lastValue(c.result(), "x"), std::atan(0.5 * std::exp(c.result().time.back())), 1e-5);

	// A function whose argument is itself a function: z = asin(sin(u)).
	const char *nested =
		"model M\n  Real u(start = 0.1);\n  Real z;\nequation\n  der(u) = 0.2;\n  z = asin(sin(u));\nend M;\n";
	Simulation d(job(nested, 2.0));
	CHECK(d.run());
	CHECK_NEAR(lastValue(d.result(), "z"), 0.5, 1e-6);
}

// Parameters and event conditions use the same functions.
static void testFunctionsInParametersAndEvents() {
	const char *src =
		"model M\n"
		"  parameter Real p = tanh(0.5) + log10(100);\n"
		"  discrete Real a(start = 0);\n"
		"  Real x(start = 2);\n"
		"equation\n"
		"  der(x) = -x;\n"
		"  when sinh(x) < 1 then\n"
		"    a = p;\n"
		"  end when;\n"
		"end M;\n";
	Simulation sim(job(src, 3.0));
	CHECK(sim.run());
	CHECK_EQ(sim.result().events.size(), (size_t)1);
	if (!sim.result().events.empty()) {
		const SasEventRecord &ev = sim.result().events[0];
		CHECK_NEAR(ev.time, std::log(2.0 / std::asinh(1.0)), 1e-6);
		CHECK_EQ(ev.assignments.size(), (size_t)1);
		if (!ev.assignments.empty()) {
			CHECK_NEAR(ev.assignments[0].second, std::tanh(0.5) + 2.0, 1e-9);
		}
	}
}

// noEvent(e) in a continuous equation is accepted; the switch is still located as an event.
static void testNoEvent() {
	const char *src =
		"model M\n"
		"  Real x(start = 1);\n"
		"equation\n"
		"  der(x) = noEvent(if x > 0.5 then -x else -2*x);\n"
		"end M;\n";
	Simulation sim(job(src, 2.0));
	CHECK(sim.run());
	double t = sim.result().time.back();
	CHECK_NEAR(lastValue(sim.result(), "x"), 0.5 * std::exp(-2.0 * (t - std::log(2.0))), 1e-4);
}

int main() {
	testBareFunctionCall();
	testNoEvent();
	testElementaryFunctions();
	testFunctionsInOdes();
	testFunctionsInParametersAndEvents();
	return sastest::finish("modelica_language_test");
}
