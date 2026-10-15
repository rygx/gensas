// End-to-end tests of the time-event simulator: parse -> check -> simulate, compared with closed-form
// solutions of piecewise-linear models.
#include "sas/HybridSimulator.h"
#include "sas/ModelChecker.h"
#include "sas/SasExpr.h"
#include "sas/SasLexico.h"
#include "sas/test/test_util.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

using namespace che::core;
using std::shared_ptr;
using std::string;

static SolverFeatures hybridFeatures() {
	SolverFeatures f;
	f.discreteVariables = f.ifExpressions = f.ifEquations = f.whenEquations = true;
	f.piecewiseFunctions = f.eventFunctions = f.simulationTime = true;
	f.timeEventsOnly = true;
	f.stateEvents = true;
	return f;
}

struct Run {
	HybridSimulator *sim = NULL;
	SasSolutionSet *sol = NULL;
	~Run() {
		delete sol;
		delete sim;
	}
};

static bool simulate(Run &run, const string &src, double endTime) {
	SasModelParser parser;
	parser.parseString(src.c_str());
	shared_ptr<SasModel> model = GlobalPool::getInstance().createNewTopModel();
	parser.currentModel = model;
	parseSingleSasModel(&parser);
	CHECK_EQ(parser.errCount, 0);
	std::ostringstream os;
	int nErr = checkModel(*model, hybridFeatures(), "test.mo", os);
	if (nErr != 0) {
		std::cerr << os.str();
	}
	CHECK_EQ(nErr, 0);
	if (parser.errCount != 0 || nErr != 0) {
		return false;
	}
	run.sim = new HybridSimulator(model);
	SasComputationOptions options;
	options.endTime = endTime;
	options.segment = 1.0;
	options.alphaTol = 1e-4;
	options.errorTol = 1e-6;
	options.nLvl = 15;
	run.sol = run.sim->run(options);
	if (run.sol == NULL) {
		std::cerr << "simulation failed: " << run.sim->error() << std::endl;
	}
	CHECK(run.sol != NULL);
	return run.sol != NULL;
}

static double finalValue(Run &run, const char *name) {
	std::map<string, double>::const_iterator it = run.sim->finalValues().find(name);
	CHECK(it != run.sim->finalValues().end());
	return it != run.sim->finalValues().end() ? it->second : NAN;
}

// A discrete gain changes at t = 1: x(2) = exp(-1) * exp(-3).
static void testDiscreteGainSwitch() {
	Run run;
	if (!simulate(run,
				  "model M\n"
				  "  discrete Real a(start = 1);\n"
				  "  Real x(start = 1);\n"
				  "equation\n"
				  "  der(x) = -a*x;\n"
				  "  when time >= 1 then\n"
				  "    a = 3;\n"
				  "  end when;\n"
				  "end M;\n",
				  2.0)) {
		return;
	}
	CHECK_EQ(run.sim->events().size(), (size_t)1);
	if (run.sim->events().size() == 1) {
		CHECK_NEAR(run.sim->events()[0].time, 1.0, 1e-12);
	}
	CHECK(run.sim->segmentCount() >= 2);
	CHECK_NEAR(finalValue(run, "a"), 3.0, 0.0);
	CHECK_NEAR(finalValue(run, "x"), std::exp(-4.0), 1e-4);
	CHECK_NEAR(run.sol->solutionLink.back()->absEnd, 2.0, 1e-9);
	// Consecutive segments cover the whole horizon without gaps.
	double prevEnd = 0.0;
	for (const shared_ptr<SasSolution> &s : run.sol->solutionLink) {
		CHECK_NEAR(s->absStart, prevEnd, 1e-9);
		prevEnd = s->absEnd;
	}
}

// reinit at sample instants: x doubles at t = 1 and t = 2; x(2.5) = 4 exp(-2.5).
static void testReinitAtSamples() {
	Run run;
	if (!simulate(run,
				  "model M\n"
				  "  Real x(start = 1);\n"
				  "equation\n"
				  "  der(x) = -x;\n"
				  "  when sample(1, 1) then\n"
				  "    reinit(x, 2*x);\n"
				  "  end when;\n"
				  "end M;\n",
				  2.5)) {
		return;
	}
	CHECK_EQ(run.sim->events().size(), (size_t)2);
	CHECK_NEAR(finalValue(run, "x"), 4.0 * std::exp(-2.5), 1e-3);
}

// if-equation with a discrete condition: der(x) = -x until t = 1, then der(x) = -2x.
static void testIfEquationSwitch() {
	Run run;
	if (!simulate(run,
				  "model M\n"
				  "  Boolean slow(start = true);\n"
				  "  Real x(start = 1);\n"
				  "equation\n"
				  "  if slow then\n"
				  "    der(x) = -x;\n"
				  "  else\n"
				  "    der(x) = -2*x;\n"
				  "  end if;\n"
				  "  when time >= 1 then\n"
				  "    slow = false;\n"
				  "  end when;\n"
				  "end M;\n",
				  2.0)) {
		return;
	}
	CHECK_NEAR(finalValue(run, "x"), std::exp(-3.0), 1e-4);
	CHECK_NEAR(finalValue(run, "slow"), 0.0, 0.0);
}

// if-expression selected by an Integer counter that advances at every sample.
static void testIfExpressionWithCounter() {
	Run run;
	if (!simulate(run,
				  "model M\n"
				  "  Integer n(start = 0);\n"
				  "  Real x(start = 1);\n"
				  "equation\n"
				  "  der(x) = if n < 1 then -x else 0;\n"
				  "  when sample(1, 1) then\n"
				  "    n = pre(n) + 1;\n"
				  "  end when;\n"
				  "end M;\n",
				  3.0)) {
		return;
	}
	CHECK_NEAR(finalValue(run, "n"), 2.0, 0.0);
	CHECK_NEAR(finalValue(run, "x"), std::exp(-1.0), 1e-4); // frozen after the first sample
}

// A DAE: algebraic y = a*x with a discrete gain.
static void testAlgebraicVariable() {
	Run run;
	if (!simulate(run,
				  "model M\n"
				  "  discrete Real a(start = 1);\n"
				  "  Real x(start = 1);\n"
				  "  Real y(start = 1);\n"
				  "equation\n"
				  "  y = a*x;\n"
				  "  der(x) = -y;\n"
				  "  when time >= 1 then\n"
				  "    a = 3;\n"
				  "  end when;\n"
				  "end M;\n",
				  2.0)) {
		return;
	}
	CHECK_NEAR(finalValue(run, "x"), std::exp(-4.0), 1e-3);
	CHECK_NEAR(finalValue(run, "y"), 3.0 * std::exp(-4.0), 1e-3);
}

// Without any event the simulator behaves like a plain solve.
static void testNoEvents() {
	Run run;
	if (!simulate(run,
				  "model M\n"
				  "  Real x(start = 1);\n"
				  "equation\n"
				  "  der(x) = -x;\n"
				  "  when time >= 10 then\n"
				  "    reinit(x, 0);\n"
				  "  end when;\n"
				  "end M;\n",
				  2.0)) {
		return;
	}
	CHECK(run.sim->events().empty()); // the event lies beyond the horizon
	CHECK_NEAR(finalValue(run, "x"), std::exp(-2.0), 1e-4);
}

// Bouncing ball with restitution e = 0.5 (reinit at the state event h < 0). The bounces are at
// t1 = sqrt(2h0/g) and t1 + 2 e t1.
static void testBouncingBall() {
	Run run;
	if (!simulate(run,
				  "model Ball\n"
				  "  parameter Real g = 9.81;\n"
				  "  parameter Real e = 0.5;\n"
				  "  Real h(start = 1);\n"
				  "  Real v(start = 0);\n"
				  "equation\n"
				  "  der(h) = v;\n"
				  "  der(v) = -g;\n"
				  "  when h < 0 then\n"
				  "    reinit(v, -e*pre(v));\n"
				  "  end when;\n"
				  "end Ball;\n",
				  1.0)) {
		return;
	}
	const double g = 9.81;
	const double e = 0.5;
	double t1 = std::sqrt(2.0 / g);
	double t2 = t1 * (1.0 + 2.0 * e);
	CHECK_EQ(run.sim->events().size(), (size_t)2);
	if (run.sim->events().size() == 2) {
		CHECK_NEAR(run.sim->events()[0].time, t1, 1e-6);
		CHECK_NEAR(run.sim->events()[1].time, t2, 1e-6);
		CHECK_NEAR(run.sim->events()[0].reinits[0].second, e * g * t1, 1e-5);
	}
	double v2 = e * e * g * t1;
	double tau = 1.0 - t2;
	CHECK_NEAR(finalValue(run, "h"), v2 * tau - 0.5 * g * tau * tau, 1e-4);
	CHECK_NEAR(finalValue(run, "v"), v2 - g * tau, 1e-4);
}

// Thermostat: two state events (T < Tlow, T > Thigh) switching a Boolean that selects the heating mode.
static void testThermostat() {
	Run run;
	if (!simulate(run,
				  "model Thermostat\n"
				  "  parameter Real Tlow = 18;\n"
				  "  parameter Real Thigh = 22;\n"
				  "  Real T(start = 20);\n"
				  "  Boolean on(start = true);\n"
				  "equation\n"
				  "  der(T) = if on then 6 - 0.1*T else 1 - 0.1*T;\n"
				  "  when T < Tlow then\n"
				  "    on = true;\n"
				  "  elsewhen T > Thigh then\n"
				  "    on = false;\n"
				  "  end when;\n"
				  "end Thermostat;\n",
				  6.0)) {
		return;
	}
	// on: T -> 60, off: T -> 10.
	double t1 = -10.0 * std::log(38.0 / 40.0);		// 20 -> 22, switch off
	double t2 = t1 + -10.0 * std::log(8.0 / 12.0);	// 22 -> 18, switch on
	double t3 = t2 + -10.0 * std::log(38.0 / 42.0); // 18 -> 22, switch off
	CHECK_EQ(run.sim->events().size(), (size_t)3);
	if (run.sim->events().size() == 3) {
		CHECK_NEAR(run.sim->events()[0].time, t1, 1e-5);
		CHECK_NEAR(run.sim->events()[1].time, t2, 1e-5);
		CHECK_NEAR(run.sim->events()[2].time, t3, 1e-5);
	}
	CHECK_NEAR(finalValue(run, "on"), 0.0, 0.0);
	CHECK_NEAR(finalValue(run, "T"), 10.0 + 12.0 * std::exp(-0.1 * (6.0 - t3)), 1e-3);
}

// An if-expression on a state condition switches the dynamics without any action: x ramps to 1 and stays.
static void testStateDependentIfExpression() {
	Run run;
	if (!simulate(run,
				  "model M\n"
				  "  Real x(start = 0);\n"
				  "equation\n"
				  "  der(x) = if x < 1 then 1 else 0;\n"
				  "end M;\n",
				  2.0)) {
		return;
	}
	CHECK_NEAR(finalValue(run, "x"), 1.0, 1e-6);
	CHECK(run.sim->segmentCount() >= 2);
}

// when-conditions re-arm: x = sin(t) crosses 0.5 upward at pi/6 and 2 pi + pi/6.
static void testRearm() {
	Run run;
	if (!simulate(run,
				  "model M\n"
				  "  Real x(start = 0);\n"
				  "  Real v(start = 1);\n"
				  "  Integer n(start = 0);\n"
				  "equation\n"
				  "  der(x) = v;\n"
				  "  der(v) = -x;\n"
				  "  when x > 0.5 then\n"
				  "    n = pre(n) + 1;\n"
				  "  end when;\n"
				  "end M;\n",
				  8.0)) {
		return;
	}
	const double pi = 3.14159265358979323846;
	CHECK_EQ(run.sim->events().size(), (size_t)2);
	if (run.sim->events().size() == 2) {
		CHECK_NEAR(run.sim->events()[0].time, pi / 6.0, 1e-5);
		CHECK_NEAR(run.sim->events()[1].time, 2.0 * pi + pi / 6.0, 1e-5);
	}
	CHECK_NEAR(finalValue(run, "n"), 2.0, 0.0);
	CHECK_NEAR(finalValue(run, "x"), std::sin(8.0), 1e-4);
}

// Event iteration: the first event sets b, which makes the second when-condition true at the same time.
static void testEventIteration() {
	Run run;
	if (!simulate(run,
				  "model M\n"
				  "  Real x(start = 0);\n"
				  "  Boolean b(start = false);\n"
				  "  discrete Real c(start = 0);\n"
				  "equation\n"
				  "  der(x) = 1;\n"
				  "  when x > 1 then\n"
				  "    b = true;\n"
				  "  end when;\n"
				  "  when b then\n"
				  "    c = 5;\n"
				  "  end when;\n"
				  "end M;\n",
				  2.0)) {
		return;
	}
	CHECK_EQ(run.sim->events().size(), (size_t)2);
	if (run.sim->events().size() == 2) {
		CHECK_NEAR(run.sim->events()[0].time, 1.0, 1e-6);
		CHECK_NEAR(run.sim->events()[1].time, run.sim->events()[0].time, 1e-12);
	}
	CHECK_NEAR(finalValue(run, "c"), 5.0, 0.0);
}

// Without damping loss the ball bounces infinitely often in finite time (Zeno behavior). The simulator
// must stop with an error instead of falling through the floor.
static void testZenoIsReported() {
	SasModelParser parser;
	parser.parseString(
		"model Ball\n"
		"  Real h(start = 1);\n"
		"  Real v(start = 0);\n"
		"equation\n"
		"  der(h) = v;\n"
		"  der(v) = -9.81;\n"
		"  when h < 0 then\n"
		"    reinit(v, -0.5*pre(v));\n"
		"  end when;\n"
		"end Ball;\n");
	shared_ptr<SasModel> model = GlobalPool::getInstance().createNewTopModel();
	parser.currentModel = model;
	parseSingleSasModel(&parser);
	CHECK_EQ(parser.errCount, 0);
	HybridSimulator sim(model);
	SasComputationOptions options;
	options.endTime = 3.0; // the accumulation point is at about 1.35
	SasSolutionSet *sol = sim.run(options);
	CHECK(sol == NULL);
	CHECK(sim.error().find("Zeno") != string::npos);
	delete sol;
}

int main() {
	testZenoIsReported();
	testBouncingBall();
	testThermostat();
	testStateDependentIfExpression();
	testRearm();
	testEventIteration();
	testDiscreteGainSwitch();
	testReinitAtSamples();
	testIfEquationSwitch();
	testIfExpressionWithCounter();
	testAlgebraicVariable();
	testNoEvents();
	return sastest::finish("hybrid_sim_test");
}
