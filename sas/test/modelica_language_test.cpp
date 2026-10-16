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

// atan2 is rewritten into half-angle atan forms; the switches between them are state events.
static void testAtan2() {
	const double pi = 3.14159265358979323846;
	// z = atan2(sin(u), cos(u)) is u wrapped into (-pi, pi]: a smooth stretch through x = 0 ...
	const char *smooth =
		"model M\n  Real u(start = -1);\n  Real z;\nequation\n  der(u) = 1;\n  z = atan2(sin(u), cos(u));\nend M;\n";
	Simulation a(job(smooth, 3.5));
	CHECK(a.run());
	CHECK_NEAR(lastValue(a.result(), "z"), -1.0 + a.result().time.back(), 1e-4);
	// ... and across the branch cut at pi, where the value jumps by 2 pi.
	const char *wrap =
		"model M\n  Real u(start = 2);\n  Real z;\nequation\n  der(u) = 1;\n  z = atan2(sin(u), cos(u));\nend M;\n";
	Simulation b(job(wrap, 2.5));
	CHECK(b.run());
	CHECK_NEAR(lastValue(b.result(), "z"), 2.0 + b.result().time.back() - 2.0 * pi, 1e-4);
	// A parameter and a condition use the numeric function.
	const char *num =
		"model M\n  parameter Real p = atan2(1, -1);\n  Real x(start = 1);\nequation\n  der(x) = -x*p;\nend M;\n";
	Simulation c(job(num, 1.0));
	CHECK(c.run());
	CHECK_NEAR(lastValue(c.result(), "x"), std::exp(-0.75 * pi * c.result().time.back()), 1e-4);
}

// when initial() runs once at the start time.
static void testWhenInitial() {
	const char *src =
		"model M\n"
		"  discrete Real a(start = 0);\n"
		"  Real x(start = 1);\n"
		"equation\n"
		"  der(x) = -a*x - x;\n"
		"  when initial() then\n"
		"    reinit(x, 2);\n"
		"    a = 3;\n"
		"  elsewhen time >= 1 then\n"
		"    a = 0;\n"
		"  end when;\n"
		"end M;\n";
	Simulation sim(job(src, 2.0));
	CHECK(sim.run());
	CHECK_EQ(sim.result().events.size(), (size_t)2);
	if (sim.result().events.size() == 2) {
		const SasEventRecord &first = sim.result().events[0];
		CHECK_NEAR(first.time, 0.0, 1e-12);
		CHECK_EQ(first.reinits.size(), (size_t)1);
		CHECK_EQ(first.assignments.size(), (size_t)1);
		CHECK_NEAR(sim.result().events[1].time, 1.0, 1e-9);
	}
	// x starts at 2, decays with rate 4 until t = 1, then with rate 1.
	double t = sim.result().time.back();
	CHECK_NEAR(lastValue(sim.result(), "x"), 2.0 * std::exp(-4.0) * std::exp(-(t - 1.0)), 1e-4);

	// `or` in conditions (it was missing from the keyword table).
	const char *either =
		"model M\n  Real x(start = 1);\n  discrete Real a(start = 0);\nequation\n  der(x) = -x;\n  when x < 0.5 or x > 3 then\n    a = 1;\n  end when;\nend M;\n";
	Simulation o(job(either, 2.0));
	CHECK(o.run());
	CHECK_EQ(o.result().events.size(), (size_t)1);
	if (!o.result().events.empty()) {
		CHECK_NEAR(o.result().events[0].time, std::log(2.0), 1e-6);
	}

	// initial() is only a when-condition by itself.
	Simulation bad(job("model M\n  Real x(start = 1);\n  discrete Real a(start = 0);\nequation\n  der(x) = -x;\n  when initial() or x < 0.5 then\n    a = 1;\n  end when;\nend M;\n", 1.0));
	CHECK(!bad.run());
	CHECK(bad.error().find("initial()") != std::string::npos);
}

// terminate() ends the run successfully at the event; assert() is checked along the solution.
static void testTerminateAndAssert() {
	const char *stop =
		"model M\n  Real x(start = 1);\nequation\n  der(x) = -x;\n"
		"  when x < 0.3 then\n    terminate(\"done\");\n  end when;\nend M;\n";
	Simulation a(job(stop, 5.0));
	CHECK(a.run());
	CHECK(a.result().terminated);
	CHECK(a.result().terminateMessage == "done");
	CHECK_NEAR(a.result().endTime, std::log(1.0 / 0.3), 1e-6);
	CHECK_NEAR(a.result().time.back(), a.result().endTime, 1e-9);
	CHECK_NEAR(lastValue(a.result(), "x"), 0.3, 1e-5);

	// A continuous assertion that becomes false stops the run with an error naming the time.
	const char *bad =
		"model M\n  Real x(start = 1);\nequation\n  der(x) = -x;\n  assert(x > 0.2, \"x too small\");\nend M;\n";
	Simulation b(job(bad, 3.0));
	CHECK(!b.run());
	CHECK(b.error().find("assertion violated") != std::string::npos);
	CHECK(b.error().find("x too small") != std::string::npos);
	CHECK(b.error().find("t=1.609") != std::string::npos); // ln 5

	// At warning level the run continues and the violation is a warning (once).
	const char *warn =
		"model M\n  Real x(start = 1);\nequation\n  der(x) = -x;\n  assert(x > 0.2, \"x low\", AssertionLevel.warning);\nend M;\n";
	Simulation w(job(warn, 3.0));
	CHECK(w.run());
	CHECK(!w.result().terminated);
	int found = 0;
	for (const std::string &m : w.warnings()) {
		found += m.find("x low") != std::string::npos ? 1 : 0;
	}
	CHECK_EQ(found, 1);

	// A true assertion is silent; one that is false at the start fails immediately; one in a when-body is checked at the event.
	Simulation ok(job("model M\n  Real x(start = 1);\nequation\n  der(x) = -x;\n  assert(x > -1, \"never\");\nend M;\n", 1.0));
	CHECK(ok.run());
	Simulation start(job("model M\n  Real x(start = 0.1);\nequation\n  der(x) = -x;\n  assert(x > 0.5, \"at start\");\nend M;\n", 1.0));
	CHECK(!start.run());
	CHECK(start.error().find("at start") != std::string::npos);
	Simulation inWhen(job("model M\n  Real x(start = 1);\nequation\n  der(x) = -x;\n  when time >= 1 then\n    assert(x > 0.5, \"in when\");\n  end when;\nend M;\n", 3.0));
	CHECK(!inWhen.run());
	CHECK(inWhen.error().find("in when") != std::string::npos);

	// Misuse is reported when the model is checked.
	Simulation top(job("model M\n  Real x(start = 1);\nequation\n  der(x) = -x;\n  terminate(\"now\");\nend M;\n", 1.0));
	CHECK(!top.run());
	CHECK(top.error().find("terminate() can only be used in a when-equation") != std::string::npos);
	Simulation noMessage(job("model M\n  Real x(start = 1);\nequation\n  der(x) = -x;\n  assert(x > 0);\nend M;\n", 1.0));
	CHECK(!noMessage.run());
	Simulation inIf(job("model M\n  parameter Boolean b = true;\n  Real x(start = 1);\nequation\n  der(x) = -x;\n  if b then\n    assert(x > 0, \"m\");\n  end if;\nend M;\n", 1.0));
	CHECK(!inIf.run());
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

// The series of exp is exact to high order: with a loose error tolerance the result is still accurate,
// which it is not when the higher-order coefficients are wrong (the step is then limited by the tolerance).
static void testExpSeriesAccuracy() {
	const char *src = "model M\n  Real u(start = 0);\n  Real z;\nequation\n  der(u) = 1;\n  z = exp(u);\nend M;\n";
	ModelicaJob j = job(src, 2.0);
	j.options.errorTol = 1e-3;
	Simulation sim(j);
	CHECK(sim.run());
	CHECK_NEAR(lastValue(sim.result(), "z"), std::exp(2.0), 1e-6);
}

// x^y: small integer powers become products, other powers go through exp and log.
static void testPowerOperator() {
	// der(x) = -x^2 has the solution x = 1/(1 + t).
	Simulation a(job("model M\n  Real x(start = 1);\nequation\n  der(x) = -x^2;\nend M;\n", 2.0));
	CHECK(a.run());
	CHECK_NEAR(lastValue(a.result(), "x"), 1.0 / 3.0, 1e-4);
	const char *src =
		"model M\n  parameter Real q = 2;\n  Real u(start = 0.5);\n  Real c, r, n, h, b, v, e, p;\nequation\n"
		"  der(u) = 0.7;\n  c = u^3;\n  r = u^q;\n  n = u^(-2);\n  h = u^0.5;\n  b = 2^u;\n  v = u^u;\n  e = (u + 1)^2;\n  p = u^2.5;\nend M;\n";
	Simulation sim(job(src, 3.0));
	CHECK(sim.run());
	const double u = 0.5 + 0.7 * sim.result().time.back();
	CHECK_NEAR(lastValue(sim.result(), "c"), std::pow(u, 3), 1e-4);
	CHECK_NEAR(lastValue(sim.result(), "r"), std::pow(u, 2), 1e-4);
	CHECK_NEAR(lastValue(sim.result(), "n"), std::pow(u, -2), 1e-4);
	CHECK_NEAR(lastValue(sim.result(), "h"), std::sqrt(u), 1e-4);
	CHECK_NEAR(lastValue(sim.result(), "b"), std::pow(2.0, u), 1e-4);
	CHECK_NEAR(lastValue(sim.result(), "v"), std::pow(u, u), 1e-4);
	CHECK_NEAR(lastValue(sim.result(), "e"), std::pow(u + 1, 2), 1e-4);
	CHECK_NEAR(lastValue(sim.result(), "p"), std::pow(u, 2.5), 1e-4);
}

// A function of constants in an equation is a constant.
static void testConstantCalls() {
	const char *src =
		"model M\n  parameter Real q = 2;\n  Real u(start = 0.5);\n  Real a, b, c;\nequation\n"
		"  der(u) = 0.7;\n  a = exp(1)*u;\n  b = sin(0.5) + u + exp(q);\n  c = atan2(1, 2) + sqrt(2) + 2^2.5 + 1/(q - 2 + 1);\nend M;\n";
	Simulation sim(job(src, 3.0));
	CHECK(sim.run());
	const double u = 0.5 + 0.7 * sim.result().time.back();
	CHECK_NEAR(lastValue(sim.result(), "a"), std::exp(1.0) * u, 1e-4);
	CHECK_NEAR(lastValue(sim.result(), "b"), std::sin(0.5) + u + std::exp(2.0), 1e-4);
	CHECK_NEAR(lastValue(sim.result(), "c"), std::atan2(1.0, 2.0) + std::sqrt(2.0) + std::pow(2.0, 2.5) + 1.0, 1e-4);
}

// Algorithm sections are turned into equations.
static void testAlgorithmSections() {
	// Sequential assignments, a reassigned variable and for-loops.
	const char *seq =
		"model M\n  parameter Integer n = 4;\n  Real x(start = 1);\n  Real y, s, w;\n  Real v[4];\nequation\n  der(x) = -s;\nalgorithm\n"
		"  y := 2*x;\n  s := y + 1;\n  s := s*0.5;\n  for i in 1:n loop\n    v[i] := i*x;\n  end for;\n  w := 0;\n  for i in 1:n loop\n    w := w + v[i];\n  end for;\nend M;\n";
	Simulation a(job(seq, 2.0));
	CHECK(a.run());
	const double xa = 1.5 * std::exp(-2.0) - 0.5; // x' = -(2x + 1)/2
	CHECK_NEAR(lastValue(a.result(), "x"), xa, 1e-4);
	CHECK_NEAR(lastValue(a.result(), "w"), 10.0 * xa, 1e-3);
	CHECK_NEAR(lastValue(a.result(), "v[3]"), 3.0 * xa, 1e-3);

	// if / elseif / else; the branch that is not taken is not evaluated (log of a positive number only
	// where x < 1), the switches are state events, and a when-statement assigns discrete variables.
	const char *br =
		"model M\n  Real x(start = 0.5);\n  Real y, z;\n  discrete Real cnt(start = 0);\n  discrete Real hit(start = 0);\nequation\n  der(x) = 1;\nalgorithm\n"
		"  if x < 1 then\n    y := log(x);\n    z := 1;\n  elseif x < 2 then\n    y := 0;\n    z := 2;\n  else\n    y := x;\n    z := 3;\n  end if;\n"
		"  when x > 1.5 then\n    cnt := cnt + 1;\n    hit := 5;\n  end when;\n  assert(x < 10, \"too large\");\nend M;\n";
	Simulation b(job(br, 3.0));
	CHECK(b.run());
	CHECK_NEAR(lastValue(b.result(), "y"), 3.5, 1e-4);
	CHECK_NEAR(lastValue(b.result(), "z"), 3.0, 1e-9);
	CHECK_EQ(b.result().events.size(), (size_t)1); // the when at x = 1.5
	bool seen = false;
	for (size_t k = 0; k < b.result().discreteNames.size(); k++) {
		if (b.result().discreteNames[k] == "cnt") {
			CHECK_NEAR(b.result().discrete[k].back(), 1.0, 1e-12); // cnt := cnt + 1 reads pre(cnt)
			seen = true;
		}
		if (b.result().discreteNames[k] == "hit") {
			CHECK_NEAR(b.result().discrete[k].back(), 5.0, 1e-12);
		}
	}
	CHECK(seen);

	// A violated assert stops the run.
	Simulation c(job(br, 12.0));
	CHECK(!c.run());
	CHECK(c.error().find("too large") != string::npos);
}

// User-defined functions are inlined at their calls.
static void testUserFunctions() {
	const char *src =
		"function sq\n  input Real a;\n  output Real y;\nalgorithm\n  y := a*a;\nend sq;\n"
		"model M\n"
		"  function poly\n    input Real a;\n    input Real b = 2;\n    output Real y;\n  protected\n    Real t;\n    Real s = 1;\n  algorithm\n"
		"    t := sq(a) + s;\n    if t > 3 then\n      y := t - b;\n    else\n      y := t + b;\n    end if;\n  end poly;\n"
		"  parameter Real k = 3;\n"
		"  function scale\n    input Real a;\n    output Real y;\n  algorithm\n    y := k*a;\n  end scale;\n"
		"  parameter Real p = sq(2) + 1;\n"
		"  Real x(start = 0.5);\n  Real y, z, w;\nequation\n  der(x) = 1;\n  y = poly(x);\n  z = scale(sq(x)) + poly(1, 5);\n  w = p*x;\nend M;\n";
	Simulation sim(job(src, 2.0));
	CHECK(sim.run());
	const double x = 0.5 + sim.result().time.back();
	const double t = x * x + 1;
	CHECK_NEAR(lastValue(sim.result(), "y"), t > 3 ? t - 2 : t + 2, 1e-4);
	CHECK_NEAR(lastValue(sim.result(), "z"), 3 * x * x + 7, 1e-3);
	CHECK_NEAR(lastValue(sim.result(), "w"), 5 * x, 1e-4);
	// A branch that is not taken is not evaluated: log is only called where its argument is positive.
	const char *guarded =
		"model M\n  function safeLog\n    input Real a;\n    output Real y;\n  algorithm\n    if a > 0 then\n      y := log(a);\n    else\n      y := -10;\n    end if;\n  end safeLog;\n"
		"  Real x(start = 2);\n  Real z;\nequation\n  der(x) = -1;\n  z = safeLog(x);\nend M;\n";
	Simulation g(job(guarded, 3.0));
	CHECK(g.run());
	CHECK_NEAR(lastValue(g.result(), "z"), -10.0, 1e-9);
}

// Packages with constants and functions, several models in one source, and which one is simulated.
// Components: instances of models, blocks and records, flattened to variables named "a.x".
static void testComponents() {
	const char *src =
		"model Decay\n  parameter Real k = 1;\n  Real x(start = 1);\n  Real y;\nequation\n  der(x) = -k*x;\n  y = 2*x;\nend Decay;\n"
		"model Pair\n  Decay a(k = 2), b;\n  Real s;\nequation\n  s = a.y + b.y;\nend Pair;\n"
		"model Top\n  parameter Real q = 3;\n  Pair p(a(k = q), b(x(start = 4)));\n  Real total;\nequation\n  total = p.s;\nend Top;\n";
	Simulation sim(job(src, 1.0));
	CHECK(sim.run());
	CHECK_NEAR(lastValue(sim.result(), "p.a.x"), std::exp(-3.0), 1e-5);
	CHECK_NEAR(lastValue(sim.result(), "p.b.x"), 4.0 * std::exp(-1.0), 1e-5);
	CHECK_NEAR(lastValue(sim.result(), "total"), 2 * std::exp(-3.0) + 8 * std::exp(-1.0), 1e-4);
	CHECK_EQ(sim.modelInfo().name, string("Top"));
	// the job can change the parameters of a component, and wins over the modification in the source
	ModelicaJob changed = job(src, 1.0);
	changed.parameters["p.b.k"] = 0.5;
	changed.parameters["p.a.k"] = 1.0;
	Simulation other(changed);
	CHECK(other.run());
	CHECK_NEAR(lastValue(other.result(), "p.b.x"), 4.0 * std::exp(-0.5), 1e-5);
	CHECK_NEAR(lastValue(other.result(), "p.a.x"), std::exp(-1.0), 1e-5);
	// a record as a parameter, a component used in the declaration of another value
	Simulation rec(job("record Pt\n  Real x;\n  Real y;\nend Pt;\n"
					   "model M\n  parameter Pt c(x = 2, y = 3);\n  parameter Real r = c.x*c.y;\n  Real v(start = 0);\nequation\n  der(v) = r;\nend M;\n",
					   1.0));
	CHECK(rec.run());
	CHECK_NEAR(lastValue(rec.result(), "v"), 6.0, 1e-6);
	// a block whose input is driven by an equation of the model around it; arrays and loops inside a component
	Simulation blk(job("block Gain\n  parameter Real k = 1;\n  input Real u;\n  output Real y;\nequation\n  y = k*u;\nend Gain;\n"
					   "model Vec\n  parameter Integer n = 3;\n  Real w[n](each start = 1);\n  Real tot;\nequation\n  for i in 1:n loop\n    der(w[i]) = -i*w[i];\n  end for;\n  tot = sum(w);\nend Vec;\n"
					   "model M\n  Gain g(k = 3);\n  Vec v;\n  Real t(start = 0), z;\nequation\n  g.u = 2 + t;\n  der(t) = g.y*0.1;\n  z = v.w[2] + g.k;\nend M;\n",
					   1.0));
	CHECK(blk.run());
	CHECK_NEAR(lastValue(blk.result(), "v.w[2]"), std::exp(-2.0), 1e-6);
	CHECK_NEAR(lastValue(blk.result(), "v.tot"), std::exp(-1.0) + std::exp(-2.0) + std::exp(-3.0), 1e-5);
	CHECK_NEAR(lastValue(blk.result(), "g.y"), 3 * (2 + lastValue(blk.result(), "t")), 1e-9);
	CHECK_NEAR(lastValue(blk.result(), "z"), std::exp(-2.0) + 3.0, 1e-5);
	// packages and functions inside a component keep their meaning
	Simulation glob(job("package P\n  constant Real c = 2;\n  function f\n    input Real a;\n    output Real r;\n  algorithm\n    r := a*c;\n  end f;\nend P;\n"
						"model Sub\n  Real x(start = 0);\nequation\n  der(x) = P.f(1) + P.c;\nend Sub;\n"
						"model M\n  Sub s;\nend M;\n",
						1.0));
	CHECK(glob.run());
	CHECK_NEAR(lastValue(glob.result(), "s.x"), 4.0, 1e-6);
	// when and discrete variables inside a component
	Simulation ev(job("model Cnt\n  discrete Integer n(start = 0);\nequation\n  when time > 0.5 then\n    n = pre(n) + 1;\n  end when;\nend Cnt;\n"
					  "model M\n  Cnt c;\nend M;\n",
					  1.0));
	CHECK(ev.run());
	bool sawN = false;
	for (size_t k = 0; k < ev.result().discreteNames.size(); k++) {
		if (ev.result().discreteNames[k] == "c.n") {
			CHECK_NEAR(ev.result().discrete[k].back(), 1.0, 1e-12);
			sawN = true;
		}
	}
	CHECK(sawN);
	// errors
	auto fails = [&](const char *text, const char *needle) {
		Simulation bad(job(text, 1.0));
		CHECK(!bad.check());
		CHECK(bad.error().find(needle) != string::npos);
	};
	const char *decay = "model Decay\n  parameter Real k = 1;\n  Real x(start = 1);\nequation\n  der(x) = -k*x;\nend Decay;\n";
	fails((string(decay) + "model M\n  Decay a(zz = 2);\nend M;\n").c_str(), "no variable or parameter 'zz'");
	fails((string(decay) + "model M\n  Decay a, a;\nend M;\n").c_str(), "Redeclaration of 'a'");
	fails((string(decay) + "model M\n  Decay a[2];\nend M;\n").c_str(), "arrays of components");
	fails("model A\n  A a;\nend A;\nmodel M\n  A a;\nend M;\n", "contains or extends itself");
	fails((string(decay) + "model M\n  Decay a(x = 3);\nend M;\n").c_str(), "not well defined"); // x = 3 is an equation, and x already has one
	fails((string(decay) + "model M\n  Decay a;\nequation\n  a.nope = 1;\nend M;\n").c_str(), "'a.nope'");
	fails("model C\n  parameter Real k;\n  Real x(start = 1);\nequation\n  der(x) = -k*x;\nend C;\nmodel M\n  C c;\nend M;\n", "Parameter 'c.k' has no value");
	fails("model C\n  Real x(start = 1);\n  function f\n    input Real a;\n    output Real b;\n  algorithm\n    b := a;\n  end f;\nequation\n  der(x) = -x;\nend C;\nmodel M\n  C c;\nend M;\n", "function inside a class used as a component");
}

// extends and short class definitions: the base class is read where `extends` stands.
static void testExtends() {
	const char *base = "partial model Base\n  parameter Real k = 1;\n  Real x(start = 1);\nequation\n  der(x) = -k*x;\nend Base;\n";
	Simulation sim(job((string(base) + "model Derived\n  extends Base(k = 3);\n  Real y;\nequation\n  y = 2*x;\nend Derived;\n").c_str(), 1.0));
	CHECK(sim.run());
	CHECK_NEAR(lastValue(sim.result(), "x"), std::exp(-3.0), 1e-5);
	CHECK_NEAR(lastValue(sim.result(), "y"), 2 * std::exp(-3.0), 1e-5);
	CHECK_EQ(sim.modelInfo().name, string("Derived"));
	// the job changes the inherited parameter by its own name
	ModelicaJob changed = job((string(base) + "model Derived\n  extends Base(k = 3);\nend Derived;\n").c_str(), 1.0);
	changed.parameters["k"] = 2.0;
	Simulation other(changed);
	CHECK(other.run());
	CHECK_NEAR(lastValue(other.result(), "x"), std::exp(-2.0), 1e-5);
	// a base that extends another one, used as a component; start values and nested modifications
	const char *chain = "model A\n  parameter Real k = 1;\n  Real x(start = 1);\nequation\n  der(x) = -k*x;\nend A;\n"
						"model B\n  extends A(k = 2, x(start = 5));\n  Real z;\nequation\n  z = x + 1;\nend B;\n"
						"model C\n  extends B;\n  parameter Real extra = k + 1;\nend C;\n";
	Simulation comp(job((string(chain) + "model M\n  C c1, c2(k = 4);\nend M;\n").c_str(), 1.0));
	CHECK(comp.run());
	CHECK_NEAR(lastValue(comp.result(), "c1.x"), 5.0 * std::exp(-2.0), 1e-5);
	CHECK_NEAR(lastValue(comp.result(), "c2.x"), 5.0 * std::exp(-4.0), 1e-5);
	CHECK_NEAR(lastValue(comp.result(), "c2.z"), 5.0 * std::exp(-4.0) + 1.0, 1e-5);
	bool sawExtra = false;
	for (const ModelParameter &p : comp.modelInfo().parameters) {
		if (p.name == "c2.extra") {
			CHECK_NEAR(p.value, 5.0, 0); // computed from the modified k
			sawExtra = true;
		}
	}
	CHECK(sawExtra);
	// short class definitions, as components and as the model to simulate
	const char *shortSrc = "model Decay\n  parameter Real k = 1;\n  Real x(start = 1);\nequation\n  der(x) = -k*x;\nend Decay;\n"
						   "model Fast = Decay(k = 10) \"fast decay\";\n"
						   "model M\n  Fast f;\nend M;\n";
	Simulation sh(job(shortSrc, 0.5));
	CHECK(sh.run());
	CHECK_NEAR(lastValue(sh.result(), "f.x"), std::exp(-5.0), 1e-5);
	ModelicaJob top = job(shortSrc, 0.5);
	top.className = "Fast";
	Simulation st(top);
	CHECK(st.run());
	CHECK_NEAR(lastValue(st.result(), "x"), std::exp(-5.0), 1e-5);
	CHECK_EQ(st.modelInfo().models.size(), (size_t)3);
	// a record extended by a record
	Simulation rec(job("record P\n  Real a;\nend P;\nrecord Q\n  extends P;\n  Real b;\nend Q;\n"
					   "model M\n  parameter Q q(a = 2, b = 3);\n  Real v(start = 0);\nequation\n  der(v) = q.a*q.b;\nend M;\n",
					   1.0));
	CHECK(rec.run());
	CHECK_NEAR(lastValue(rec.result(), "v"), 6.0, 1e-6);
	// errors
	auto fails = [&](const string &text, const char *needle) {
		Simulation bad(job(text.c_str(), 1.0));
		CHECK(!bad.check());
		CHECK(bad.error().find(needle) != string::npos);
	};
	fails("model M\n  extends Nope;\nend M;\n", "unknown base class 'Nope'");
	fails("model M\n  extends M;\nend M;\n", "contains or extends itself");
	fails("model A\n  extends B;\nend A;\nmodel B\n  extends A;\nend B;\nmodel M\n  extends A;\nend M;\n", "contains or extends itself");
	fails(string(base) + "model M\n  extends Base(zz = 1);\nend M;\n", "no variable or parameter 'zz'");
	fails(string(base) + "model M\n  Real y;\nequation\n  y = 1;\n  extends Base;\nend M;\n", "`extends` must come before");
	fails(string(base) + "model M\n  Base b;\nend M;\n", "partial class 'Base'");
	fails("model A\n  function f\n    input Real a;\n    output Real b;\n  algorithm\n    b := a;\n  end f;\n  Real x(start = 1);\nequation\n  der(x) = -x;\nend A;\nmodel M\n  extends A;\nend M;\n", "function inside a class used as a component");
}

// Types made with `type T = Real(...)`: unit and start values come with the type, the variable can change them.
static void testTypeAliases() {
	const char *src =
		"type Volt = Real(unit = \"V\", start = 5);\n"
		"type Level = Volt(start = 2) \"a level\";\n"
		"type Count = Integer;\n"
		"package Phys\n  type Rate = Real(unit = \"1/s\");\n  type Slow = Rate(start = 3);\nend Phys;\n"
		"model M\n"
		"  parameter Phys.Rate k = 2;\n"
		"  Volt a;\n  Level b;\n  Level c(start = 7);\n  Phys.Slow d, e;\n  Volt arr[2](each start = 1);\n  discrete Count n(start = 4);\n"
		"equation\n  der(a) = -k*a;\n  der(b) = -b;\n  der(c) = -c;\n  der(d) = -d;\n  der(e) = -e;\n  der(arr) = -arr;\n  when time > 0.5 then\n    n = pre(n) + 1;\n  end when;\nend M;\n";
	Simulation sim(job(src, 1.0));
	CHECK(sim.run());
	const ModelInfo &info = sim.modelInfo();
	auto var = [&](const string &name) {
		for (const ModelVariable &v : info.variables) {
			if (v.name == name) {
				return v;
			}
		}
		CHECK(false);
		return ModelVariable();
	};
	CHECK_EQ(var("a").unit, string("V"));
	CHECK_NEAR(var("a").start, 5.0, 0);
	CHECK_NEAR(var("b").start, 2.0, 0); // the outer type changes the start of the inner one
	CHECK_EQ(var("b").unit, string("V"));
	CHECK_NEAR(var("c").start, 7.0, 0); // the variable has the last word
	CHECK_EQ(var("d").unit, string("1/s"));
	CHECK_NEAR(var("d").start, 3.0, 0);
	CHECK_NEAR(var("e").start, 3.0, 0);
	CHECK_NEAR(var("arr[2]").start, 1.0, 0);
	CHECK_EQ(info.parameters[0].unit, string("1/s"));
	CHECK_NEAR(lastValue(sim.result(), "a"), 5.0 * std::exp(-2.0), 1e-5);
	bool sawCount = false;
	for (size_t k = 0; k < sim.result().discreteNames.size(); k++) {
		if (sim.result().discreteNames[k] == "n") {
			CHECK_NEAR(sim.result().discrete[k].back(), 5.0, 1e-12);
			sawCount = true;
		}
	}
	CHECK(sawCount);
	// errors
	Simulation unknown(job("model M\n  Nope x;\nequation\n  der(x) = 1;\nend M;\n", 1.0));
	CHECK(!unknown.check());
	CHECK(unknown.error().find("unknown type 'Nope'") != string::npos);
	Simulation twice(job("type T = Real;\ntype T = Real;\nmodel M\n  Real x(start = 1);\nequation\n  der(x) = -x;\nend M;\n", 1.0));
	CHECK(!twice.check());
	CHECK(twice.error().find("defined twice") != string::npos);
	Simulation later(job("type A = B;\ntype B = Real;\nmodel M\n  Real x(start = 1);\nequation\n  der(x) = -x;\nend M;\n", 1.0));
	CHECK(!later.check());
	CHECK(later.error().find("unknown type 'B'") != string::npos);
	Simulation arrays(job("type V3 = Real[3];\nmodel M\n  Real x(start = 1);\nequation\n  der(x) = -x;\nend M;\n", 1.0));
	CHECK(!arrays.check());
	CHECK(arrays.error().find("array types") != string::npos);
	Simulation record(job("type R = Foo.Bar(x = 1);\nmodel M\n  Real x(start = 1);\nequation\n  der(x) = -x;\nend M;\n", 1.0));
	CHECK(!record.check());
	// a type after the model
	Simulation after(job("model M\n  U x(start = 1);\nequation\n  der(x) = -x;\nend M;\ntype U = Real(unit = \"m\");\n", 1.0));
	CHECK(after.run());
	CHECK_EQ(after.modelInfo().variables[0].unit, string("m"));
}

static void testClassesAndPackages() {
	const char *src =
		"package Phys\n"
		"  constant Real g = 9.81;\n  constant Real two_g = 2*g;\n  constant Integer n = 3;\n"
		"  function drag\n    input Real v;\n    output Real f;\n  algorithm\n    f := two_g*v/n;\n  end drag;\n"
		"  function twice\n    input Real v;\n    output Real f;\n  algorithm\n    f := drag(v) + drag(v);\n  end twice;\n"
		"  package Inner\n    constant Real c = 0.5*g;\n    function half\n      input Real v;\n      output Real y;\n    algorithm\n      y := c*v;\n    end half;\n  end Inner;\n"
		"end Phys;\n"
		"model Unused\n  Real q(start = 1);\nequation\n  der(q) = 1;\nend Unused;\n"
		"model Fall\n  parameter Real m = Phys.g / 9.81;\n  Real v(start = 0);\n  Real h(start = 10);\n  Real d, e, f;\nequation\n"
		"  der(v) = -Phys.g*m;\n  der(h) = v;\n  d = Phys.twice(2);\n  e = Phys.Inner.half(4);\n  f = Phys.Inner.c + Phys.n;\nend Fall;\n";
	Simulation sim(job(src, 1.0));
	CHECK(sim.run());
	CHECK_NEAR(lastValue(sim.result(), "v"), -9.81, 1e-4);
	CHECK_NEAR(lastValue(sim.result(), "d"), 2 * (2 * 9.81 * 2 / 3), 1e-6);
	CHECK_NEAR(lastValue(sim.result(), "e"), 0.5 * 9.81 * 4, 1e-6);
	CHECK_NEAR(lastValue(sim.result(), "f"), 0.5 * 9.81 + 3, 1e-6);
	CHECK_EQ(sim.modelInfo().name, string("Fall")); // the last model unless one is named
	CHECK_EQ(sim.modelInfo().models.size(), (size_t)2);
	// another model of the source can be chosen
	ModelicaJob chosen = job(src, 1.0);
	chosen.className = "Unused";
	Simulation other(chosen);
	CHECK(other.run());
	CHECK_NEAR(lastValue(other.result(), "q"), 2.0, 1e-6);
	CHECK_EQ(other.modelInfo().name, string("Unused"));
	// an unknown name is an error, and the list of models is still offered
	chosen.className = "Nope";
	Simulation missing(chosen);
	CHECK(!missing.check());
	CHECK(missing.error().find("no model named 'Nope'") != string::npos);
	CHECK(missing.hasModelInfo() && !missing.modelInfo().complete);
	CHECK_EQ(missing.modelInfo().models.size(), (size_t)2);
	// functions and packages may stand after the model that uses them
	const char *late =
		"model M\n  Real x(start = 1);\nequation\n  der(x) = -P.k*sq(x);\nend M;\n"
		"function sq\n  input Real a;\n  output Real y;\nalgorithm\n  y := a*a;\nend sq;\n"
		"package P\n  constant Real k = 2;\nend P;\n";
	Simulation after(job(late, 0.5));
	CHECK(after.run());
	CHECK_NEAR(lastValue(after.result(), "x"), 1.0 / (1.0 + 2.0 * 0.5), 1e-4);
	// an error in a package is reported with its line
	Simulation bad(job("package P\n  constant Real a = 1;\n  constant Real a = 2;\nend P;\nmodel M\n  Real x(start = 1);\nequation\n  der(x) = -x;\nend M;\n", 1.0));
	CHECK(!bad.check());
	CHECK(bad.error().find("defined twice") != string::npos);
	Simulation unknown(job("model M\n  Real x(start = 1);\nequation\n  der(x) = -Q.k*x;\nend M;\n", 1.0));
	CHECK(!unknown.check());
	CHECK(unknown.error().find("'Q.k'") != string::npos);
	Simulation nonConstant(job("package P\n  constant Real a = b;\nend P;\nmodel M\n  Real x(start = 1);\nequation\n  der(x) = -x;\nend M;\n", 1.0));
	CHECK(!nonConstant.check());
}

// Whole-array equations are expanded elementwise.
static void testArrayExpressions() {
	const char *src =
		"model M\n  Real x[3](start = {1, 2, 3});\n  Real y[3], z[3], w[3], c[2,2];\n  Real s;\nequation\n"
		"  der(x) = -x;\n  y = 2*x + {1, 1, 1};\n  z = sin(x) - x/2;\n  w = if s > 4 then x else -x;\n  s = sum(x);\n  c = {{1, 2}, {3, 4}}*s;\nend M;\n";
	Simulation sim(job(src, 1.0));
	CHECK(sim.run());
	const double e = std::exp(-1.0);
	CHECK_NEAR(lastValue(sim.result(), "x[2]"), 2 * e, 1e-4);
	CHECK_NEAR(lastValue(sim.result(), "y[3]"), 6 * e + 1, 1e-4);
	CHECK_NEAR(lastValue(sim.result(), "z[1]"), std::sin(e) - e / 2, 1e-4);
	const double s = 6 * e;
	CHECK_NEAR(lastValue(sim.result(), "w[1]"), s > 4 ? e : -e, 1e-4);
	CHECK_NEAR(lastValue(sim.result(), "c[2,1]"), 3 * s, 1e-4);
	// Assignments in an algorithm section and reinit() of an array work the same way.
	const char *alg =
		"model M\n  Real x[2](each start = 1);\n  Real y[2];\n  discrete Real k(start = 0);\nequation\n  der(x) = -x;\n"
		"  when time >= 0.5 then\n    reinit(x, {5, 6});\n  end when;\nalgorithm\n  y := 2*x;\nend M;\n";
	Simulation b(job(alg, 1.0));
	CHECK(b.run());
	CHECK_NEAR(lastValue(b.result(), "y[2]"), 12 * std::exp(-0.5), 1e-3);
}

// A run ends exactly at the requested end time, also when the last steps are short (exp grows fast).
static void testEndTimeIsReached() {
	const char *src = "model M\n  Real u(start = -1);\n  Real z;\nequation\n  der(u) = 1;\n  z = exp(u);\nend M;\n";
	const double ends[] = {3.0, 2.5, 0.7, 1.2345};
	for (double e : ends) {
		Simulation sim(job(src, e));
		CHECK(sim.run());
		if (!sim.result().time.empty()) {
			CHECK_NEAR(sim.result().time.back(), e, 1e-12);
			CHECK_NEAR(lastValue(sim.result(), "z"), std::exp(-1.0 + e), 1e-4 * std::exp(-1.0 + e));
		}
	}
	// The same holds for a plain ODE and for a hybrid model with events.
	Simulation ode(job("model M\n  Real x(start = 1);\nequation\n  der(x) = -x;\nend M;\n", 2.0));
	CHECK(ode.run());
	CHECK_NEAR(ode.result().time.back(), 2.0, 1e-12);
	const char *hyb = "model M\n  discrete Real a(start = 1);\n  Real x(start = 1);\nequation\n  der(x) = -a*x;\n  when time >= 0.5 then\n    a = 2;\n  end when;\nend M;\n";
	Simulation h(job(hyb, 1.0));
	CHECK(h.run());
	CHECK_NEAR(h.result().time.back(), 1.0, 1e-12);
}

// A function evaluated outside its domain stops the run with a message instead of returning garbage.
static void checkDomainError(const char *name, const char *body, double u0, double rate, const char *expectInMessage) {
	string src = string("model M\n  Real u(start = ") + to_string(u0) + ");\n  Real z;\nequation\n  der(u) = " + to_string(rate) + ";\n  z = " + body + ";\nend M;\n";
	Simulation sim(job(src.c_str(), 3.0));
	if (sim.run()) {
		cerr << "FAILED " << name << ": the run was accepted" << endl;
		CHECK(false);
		return;
	}
	if (sim.error().find(expectInMessage) == string::npos || sim.error().find("cannot advance") == string::npos) {
		cerr << "FAILED " << name << ": message \"" << sim.error() << "\"" << endl;
		CHECK(false);
	}
}

static void testDomainChecks() {
	checkDomainError("log", "log(u)", 1.0, -1.0, "log");
	checkDomainError("sqrt", "sqrt(u)", 1.0, -1.0, "sqrt");
	checkDomainError("asin", "asin(u)", 0.0, 1.0, "asin");
	checkDomainError("acos", "acos(u)", 0.0, -1.0, "acos");
	checkDomainError("division", "1/u", 1.0, -1.0, "denominator");
	// The same inside a hybrid model.
	const char *hyb = "model M\n  discrete Real a(start = 1);\n  Real u(start = 1);\n  Real z;\nequation\n  der(u) = -a;\n  z = log(u);\n  when time >= 0.2 then\n    a = 2;\n  end when;\nend M;\n";
	Simulation h(job(hyb, 3.0));
	CHECK(!h.run());
	CHECK(h.error().find("log") != string::npos);
	// A guard that keeps the argument inside the domain works: the branch that is not taken is not evaluated.
	const char *guard = "model M\n  Real x(start = 2);\n  Real z, r;\nequation\n  der(x) = -1;\n  z = if x > 0 then log(x) else -10;\n  r = if x > 0 then 1/x else 0;\nend M;\n";
	Simulation gd(job(guard, 3.0));
	CHECK(gd.run());
	CHECK_NEAR(lastValue(gd.result(), "z"), -10.0, 1e-9);
	CHECK_NEAR(lastValue(gd.result(), "r"), 0.0, 1e-9);
	// A function that touches the edge of its domain without leaving it still works.
	Simulation ok(job("model M\n  Real u(start = 0);\n  Real z;\nequation\n  der(u) = 1;\n  z = sqrt(u + 1);\nend M;\n", 3.0));
	CHECK(ok.run());
	CHECK_NEAR(lastValue(ok.result(), "z"), 2.0, 1e-4);
}

int main() {
	testAlgorithmSections();
	testUserFunctions();
	testClassesAndPackages();
	testTypeAliases();
	testComponents();
	testExtends();
	testArrayExpressions();
	testExpSeriesAccuracy();
	testPowerOperator();
	testConstantCalls();
	testEndTimeIsReached();
	testDomainChecks();
	testAtan2();
	testTerminateAndAssert();
	testWhenInitial();
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
