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

// Real y = expr; is an equation, wherever y and the variables in expr are declared.
static void testBindingEquations() {
	const char *src =
		"model M\n"
		"  Real y = 2*x;\n"
		"  Real x(start = 1);\n"
		"  Real s(start = 0) = sin(x);\n"
		"equation\n"
		"  der(x) = -y;\n"
		"end M;\n";
	Simulation sim(job(src, 1.0));
	CHECK(sim.run());
	CHECK(sim.error().empty());
	double t = sim.result().time.back();
	CHECK_NEAR(lastValue(sim.result(), "x"), std::exp(-2.0 * t), 1e-4);
	CHECK_NEAR(lastValue(sim.result(), "y"), 2.0 * std::exp(-2.0 * t), 1e-4);
	CHECK_NEAR(lastValue(sim.result(), "s"), std::sin(std::exp(-2.0 * t)), 1e-4);

	// Combined with events.
	const char *hybrid =
		"model M\n"
		"  Real x(start = 1);\n"
		"  Real y = 3*x;\n"
		"  discrete Real k(start = 1);\n"
		"equation\n"
		"  der(x) = -k*x;\n"
		"  when time >= 1 then\n"
		"    k = 2;\n"
		"  end when;\n"
		"end M;\n";
	Simulation h(job(hybrid, 2.0));
	CHECK(h.run());
	t = h.result().time.back();
	double x = std::exp(-1.0) * std::exp(-2.0 * (t - 1.0));
	CHECK_NEAR(lastValue(h.result(), "x"), x, 1e-4);
	CHECK_NEAR(lastValue(h.result(), "y"), 3.0 * x, 1e-4);

	// Discrete variables are assigned in when-equations.
	Simulation bad(job("model M\n  discrete Real d = 1;\n  Real x(start = 1);\nequation\n  der(x) = -x;\nend M;\n", 1.0));
	CHECK(!bad.run());
	CHECK(bad.error().find("discrete variable 'd'") != std::string::npos);
}

// A model written the way Modelica sources usually are.
static void testSourceLayout() {
	const char *src =
		"model M \"a model\"\n"
		"  /* block\n     comment */\n"
		"  parameter Real k = 2 \"gain\"; // trailing\n"
		"  Real x(start = 1) \"state\", z(start = 2);\n"
		"  output Real y;\n"
		"  input Real u = 0.5;\n"
		"  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}})));\n"
		"equation\n"
		"  der(x) = -k*x + u; /* inline */\n"
		"  der(z) = -z;\n"
		"  y = x + z;\n"
		"  annotation(experiment(StopTime = 1, Tolerance = 1e-6));\n"
		"end M;\n";
	Simulation sim(job(src, 1.0));
	CHECK(sim.run());
	double t = sim.result().time.back();
	double x = 0.25 + 0.75 * std::exp(-2.0 * t);
	CHECK_NEAR(lastValue(sim.result(), "x"), x, 1e-4);
	CHECK_NEAR(lastValue(sim.result(), "y"), x + 2.0 * std::exp(-t), 1e-4);
}

// Implicit algebraic equations are solved with Newton's method from the start values; a missing
// start value that makes the Jacobian singular is reported instead of giving a wrong result.
static void testImplicitAlgebraicEquations() {
	const char *src =
		"model M\n  Real x(start = 4);\n  Real y(start = 1);\nequation\n  der(x) = -x;\n  y*y = x;\nend M;\n";
	Simulation sim(job(src, 1.0));
	CHECK(sim.run());
	CHECK_NEAR(lastValue(sim.result(), "y"), std::sqrt(4.0 * std::exp(-sim.result().time.back())), 1e-4);

	const char *cubic =
		"model M\n  Real x(start = 4);\n  Real y(start = 1);\nequation\n  der(x) = -x;\n  x = y*y*y;\nend M;\n";
	Simulation c(job(cubic, 1.0));
	CHECK(c.run());
	CHECK_NEAR(lastValue(c.result(), "y"), std::cbrt(4.0 * std::exp(-c.result().time.back())), 1e-4);

	const char *noStart =
		"model M\n  Real x(start = 4);\n  Real y;\nequation\n  der(x) = -x;\n  y*y = x;\nend M;\n";
	Simulation bad(job(noStart, 1.0));
	CHECK(!bad.run());
	CHECK(bad.error().find("cannot be solved consistently") != std::string::npos);
}

// Arrays are sets of scalar variables; for-equations are unrolled.
static void testArrays() {
	// A cascade of first-order lags: x1 = exp(-t), x2 = t exp(-t), x3 = t^2/2 exp(-t), x4 = t^3/6 exp(-t).
	const char *cascade =
		"model M\n"
		"  parameter Integer n = 4;\n"
		"  Real x[n](start = {1, 0, 0, 0});\n"
		"equation\n"
		"  der(x[1]) = -x[1];\n"
		"  for i in 2:n loop\n"
		"    der(x[i]) = x[i-1] - x[i];\n"
		"  end for;\n"
		"end M;\n";
	Simulation a(job(cascade, 2.0));
	CHECK(a.run());
	double t = a.result().time.back();
	CHECK_NEAR(lastValue(a.result(), "x[1]"), std::exp(-t), 1e-4);
	CHECK_NEAR(lastValue(a.result(), "x[2]"), t * std::exp(-t), 1e-4);
	CHECK_NEAR(lastValue(a.result(), "x[3]"), t * t / 2.0 * std::exp(-t), 1e-4);
	CHECK_NEAR(lastValue(a.result(), "x[4]"), t * t * t / 6.0 * std::exp(-t), 1e-4);

	// Two dimensions, nested iterators, a parameter array, sum(), size() and `each`.
	const char *matrix =
		"model M\n"
		"  parameter Real k[2] = {1, 2};\n"
		"  Real A[2,3](each start = 1);\n"
		"  Real s = sum(A) + size(A, 2);\n"
		"equation\n"
		"  for i in 1:2, j in 1:size(A, 2) loop\n"
		"    der(A[i,j]) = -k[i]*j*A[i,j];\n"
		"  end for;\n"
		"end M;\n";
	Simulation m(job(matrix, 1.0));
	CHECK(m.run());
	t = m.result().time.back();
	double total = 0.0;
	for (int i = 1; i <= 2; i++) {
		for (int j = 1; j <= 3; j++) {
			double expected = std::exp(-i * j * t); // k[i] = i
			char name[16];
			std::snprintf(name, sizeof name, "A[%d,%d]", i, j);
			CHECK_NEAR(lastValue(m.result(), name), expected, 1e-4);
			total += expected;
		}
	}
	CHECK_NEAR(lastValue(m.result(), "s"), total + 3.0, 1e-3);

	// Reverse and stepped ranges; the dimensions in front of the name apply to every declarator.
	const char *ranges =
		"model M\n"
		"  Real[2] x(each start = 1), y(each start = 2);\n"
		"equation\n"
		"  for i in 1:2 loop der(x[i]) = -x[i]; end for;\n"
		"  for i in 2:-1:1 loop der(y[i]) = -2*y[i]; end for;\n"
		"end M;\n";
	Simulation r(job(ranges, 1.0));
	CHECK(r.run());
	t = r.result().time.back();
	CHECK_NEAR(lastValue(r.result(), "x[2]"), std::exp(-t), 1e-4);
	CHECK_NEAR(lastValue(r.result(), "y[1]"), 2.0 * std::exp(-2.0 * t), 1e-4);

	// Initial equations in a loop; when-equations in a loop.
	const char *events =
		"model M\n"
		"  discrete Real a[2](each start = 0);\n"
		"  Real x[2];\n"
		"equation\n"
		"  for i in 1:2 loop\n"
		"    der(x[i]) = -a[i] - x[i];\n"
		"    when time >= i then\n"
		"      a[i] = 2*i;\n"
		"    end when;\n"
		"  end for;\n"
		"initial equation\n"
		"  for i in 1:2 loop x[i] = i; end for;\n"
		"end M;\n";
	Simulation e(job(events, 3.0));
	CHECK(e.run());
	CHECK_EQ(e.result().events.size(), (size_t)2);
	if (e.result().events.size() == 2) {
		CHECK_NEAR(e.result().events[0].time, 1.0, 1e-9);
		CHECK_NEAR(e.result().events[1].time, 2.0, 1e-9);
	}
	t = e.result().time.back();
	// x1: -x - a with a = 2 after t = 1: x1(t) = -2 + (x1(1) + 2) exp(-(t-1)), x1(1) = exp(-1)
	CHECK_NEAR(lastValue(e.result(), "x[1]"), -2.0 + (std::exp(-1.0) + 2.0) * std::exp(-(t - 1.0)), 1e-3);
}

int main() {
	testArrays();
	testSourceLayout();
	testImplicitAlgebraicEquations();
	testBindingEquations();
	testBareFunctionCall();
	testNoEvent();
	testElementaryFunctions();
	testFunctionsInOdes();
	testFunctionsInParametersAndEvents();
	return sastest::finish("modelica_language_test");
}
