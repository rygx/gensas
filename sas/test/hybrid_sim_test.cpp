// End-to-end tests of the time-event simulator: parse -> check -> simulate, compared with closed-form
// solutions of piecewise-linear models.
#include "sas/HybridSimulator.h"
#include "sas/ModelChecker.h"
#include "sas/SasExpr.h"
#include "sas/SasLexico.h"
#include "sas/SasResult.h"
#include "sas/test/test_util.h"

#include "json/json.h"

#include <cmath>
#include <cstdio>
#include <fstream>
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

static bool simulate(Run &run, const string &src, double endTime, PiecewiseMode mode = PW_EVENTS, double eps = 1e-6) {
	SasModelParser parser;
	parser.parseString(src.c_str());
	shared_ptr<SasModel> model = GlobalPool::getInstance().createNewTopModel();
	parser.currentModel = model;
	parseSingleSasModel(&parser);
	CHECK_EQ(parser.errCount, 0);
	rewritePiecewiseFunctions(*model, mode, eps);
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
// A variable that has to stay at a switching point chatters: v' = -sign(v) reaches 0 and the events pile up.
static void testSignChatterIsReported() {
	SasModelParser parser;
	parser.parseString(
		"model M\n"
		"  Real v(start = 0.5);\n"
		"equation\n"
		"  der(v) = -sign(v);\n"
		"end M;\n");
	shared_ptr<SasModel> model = GlobalPool::getInstance().createNewTopModel();
	parser.currentModel = model;
	parseSingleSasModel(&parser);
	CHECK_EQ(parser.errCount, 0);
	rewritePiecewiseFunctions(*model);
	HybridSimulator sim(model);
	SasComputationOptions options;
	options.endTime = 2.0;
	SasSolutionSet *sol = sim.run(options);
	CHECK(sol == NULL);
	CHECK(!sim.error().empty());
	delete sol;
}

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

// Writes the solution sampled every `interval` as JSON.
static void writeJson(che::core::SasSolutionSet &set, const char *fileName, double interval) {
	che::core::SasResult result = che::core::buildResult(set, interval);
	std::string err;
	CHECK(che::core::makeResultWriter("json")->write(result, fileName, err));
}

// The JSON result carries the discrete variables (sampled per solution segment, so the two samples at an
// event time show the values before and after) and the list of events.
static void testJsonResult() {
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
	const char *fileName = "hybrid_sim_result.json";
	writeJson(*run.sol, fileName, 0.25);
	std::ifstream in(fileName);
	Json::Value root;
	Json::CharReaderBuilder builder;
	std::string errs;
	CHECK(Json::parseFromStream(builder, in, &root, &errs));
	CHECK(root.isMember("time") && root.isMember("x") && root.isMember("discrete") && root.isMember("events"));
	if (!(root.isMember("time") && root.isMember("x") && root.isMember("discrete") && root.isMember("events"))) {
		return;
	}
	const Json::Value &time = root["time"];
	const Json::Value &a = root["discrete"]["a"];
	CHECK_EQ(a.size(), time.size());
	CHECK(time.size() >= 8);
	int nBefore = 0;
	int nAfter = 0;
	for (Json::ArrayIndex i = 0; i < time.size() && i < a.size(); i++) {
		double ti = time[i].asDouble();
		double ai = a[i].asDouble();
		if (ti < 1.0 - 1e-9) {
			CHECK_NEAR(ai, 1.0, 0.0);
			nBefore++;
		} else if (ti > 1.0 + 1e-9) {
			CHECK_NEAR(ai, 3.0, 0.0);
			nAfter++;
		}
	}
	CHECK(nBefore >= 3 && nAfter >= 3);
	// At t = 1 two samples exist: before (a = 1) and after (a = 3) the event.
	int atEvent = 0;
	for (Json::ArrayIndex i = 0; i < time.size() && i < a.size(); i++) {
		if (std::fabs(time[i].asDouble() - 1.0) < 1e-9) {
			CHECK_NEAR(a[i].asDouble(), atEvent == 0 ? 1.0 : 3.0, 0.0);
			atEvent++;
		}
	}
	CHECK_EQ(atEvent, 2);
	CHECK_EQ(root["events"].size(), (Json::ArrayIndex)1);
	CHECK_NEAR(root["events"][0]["time"].asDouble(), 1.0, 1e-12);
	CHECK_NEAR(root["events"][0]["assignments"]["a"].asDouble(), 3.0, 0.0);
}

// An if-equation that swaps which variable is the state: mode A has der(x) = -x and y = 2x, mode B has
// der(y) = -y and x = 0.5 y. The variable set is the same, the order and the roles differ. x = exp(-t)
// before t = 1; after it y = 2 exp(-1) exp(-(t-1)) = 2 exp(-t), so x is exp(-t) in both modes.
static void testStateAlgebraicSwap() {
	Run run;
	if (!simulate(run,
				  "model M\n"
				  "  Boolean a(start = true);\n"
				  "  Real x(start = 1);\n"
				  "  Real y(start = 2);\n"
				  "equation\n"
				  "  if a then\n"
				  "    der(x) = -x;\n"
				  "    y = 2*x;\n"
				  "  else\n"
				  "    der(y) = -y;\n"
				  "    x = 0.5*y;\n"
				  "  end if;\n"
				  "  when time >= 1 then\n"
				  "    a = false;\n"
				  "  end when;\n"
				  "end M;\n",
				  2.0)) {
		return;
	}
	CHECK_NEAR(finalValue(run, "x"), std::exp(-2.0), 1e-3);
	CHECK_NEAR(finalValue(run, "y"), 2.0 * std::exp(-2.0), 1e-3);
	// Results are matched by name across the two computation models.
	const char *fileName = "hybrid_sim_swap.json";
	writeJson(*run.sol, fileName, 0.5);
	std::ifstream in(fileName);
	Json::Value root;
	Json::CharReaderBuilder builder;
	std::string errs;
	CHECK(Json::parseFromStream(builder, in, &root, &errs));
	if (!root.isMember("x") || !root.isMember("y") || !root.isMember("time")) {
		CHECK(false);
		return;
	}
	for (Json::ArrayIndex i = 0; i < root["time"].size(); i++) {
		double ti = root["time"][i].asDouble();
		CHECK_NEAR(root["x"][i].asDouble(), std::exp(-ti), 2e-3);
		CHECK_NEAR(root["y"][i].asDouble(), 2.0 * std::exp(-ti), 4e-3);
	}
}

// A pulse that starts and ends between two guard samples (the default spacing is 1/48 s on a 1 s
// segment). The parabola x(t) = h - 400 (t - c)^2 is above zero for |t - c| < sqrt(h / 400) only.
static string pulseModel(const char *sign, const char *guard, const char *body) {
	const double c = 49.0 / 96.0; // halfway between the samples 24/48 and 25/48
	const double h = 0.005;
	char buf[1024];
	std::snprintf(buf, sizeof(buf),
				  "model M\n"
				  "  Real x(start = %.15g);\n"
				  "  Real v(start = %.15g);\n"
				  "  Real y(start = 0);\n"
				  "  discrete Integer n(start = 0);\n"
				  "equation\n"
				  "  der(x) = v;\n"
				  "  der(v) = %s800;\n"
				  "%s"
				  "%s"
				  "end M;\n",
				  *sign == '-' ? h - 400.0 * c * c : -h + 400.0 * c * c,
				  *sign == '-' ? 800.0 * c : -800.0 * c, sign, guard, body);
	return buf;
}

static void testPulseBetweenSamples() {
	Run run;
	string src = pulseModel("-",
							"  der(y) = 0;\n"
							"  when x > 0 then\n"
							"    n = pre(n) + 1;\n"
							"  end when;\n",
							"");
	if (!simulate(run, src, 1.0)) {
		return;
	}
	CHECK_EQ(run.sim->events().size(), (size_t)1);
	if (run.sim->events().size() == 1) {
		CHECK_NEAR(run.sim->events()[0].time, 49.0 / 96.0 - std::sqrt(0.005 / 400.0), 1e-9);
	}
	CHECK_NEAR(finalValue(run, "n"), 1.0, 0.0);
}

// The margin dips below zero and comes back between two samples: a when-condition re-arms and
// fires on the way up.
static void testDipBetweenSamples() {
	Run run;
	string src = pulseModel("",
							"  der(y) = 0;\n"
							"  when x > 0 then\n"
							"    n = pre(n) + 1;\n"
							"  end when;\n",
							"");
	if (!simulate(run, src, 1.0)) {
		return;
	}
	CHECK_EQ(run.sim->events().size(), (size_t)1);
	if (run.sim->events().size() == 1) {
		CHECK_NEAR(run.sim->events()[0].time, 49.0 / 96.0 + std::sqrt(0.005 / 400.0), 1e-9);
	}
	CHECK_NEAR(finalValue(run, "n"), 1.0, 0.0);
}

// An if-equation whose guard is true only during the pulse: y integrates the pulse length.
static void testPulseInIfEquation() {
	Run run;
	string src = pulseModel("-",
							"  if x > 0 then\n"
							"    der(y) = 1;\n"
							"  else\n"
							"    der(y) = 0;\n"
							"  end if;\n",
							"");
	if (!simulate(run, src, 1.0)) {
		return;
	}
	CHECK_NEAR(finalValue(run, "y"), 2.0 * std::sqrt(0.005 / 400.0), 1e-6);
}

// edge() and change() in when-conditions: b toggles at t = 0.5, 1.5, 2.5; it rises at 0.5 and 2.5.
static void testEdgeAndChange() {
	Run run;
	if (!simulate(run,
				  "model M\n"
				  "  Boolean b(start = false);\n"
				  "  Integer n(start = 0);\n"
				  "  Integer m(start = 0);\n"
				  "  Real x(start = 1);\n"
				  "equation\n"
				  "  der(x) = -x;\n"
				  "  when sample(0.5, 1) then\n"
				  "    b = not pre(b);\n"
				  "  end when;\n"
				  "  when edge(b) then\n"
				  "    n = pre(n) + 1;\n"
				  "  end when;\n"
				  "  when change(b) then\n"
				  "    m = pre(m) + 1;\n"
				  "  end when;\n"
				  "end M;\n",
				  3.0)) {
		return;
	}
	CHECK_NEAR(finalValue(run, "n"), 2.0, 0.0);
	CHECK_NEAR(finalValue(run, "m"), 3.0, 0.0);
	CHECK_NEAR(finalValue(run, "b"), 1.0, 0.0);
	CHECK_NEAR(finalValue(run, "x"), std::exp(-3.0), 1e-4);
	// The counters are updated at the toggle instants (with the toggle in the same event record).
	std::vector<double> nTimes;
	for (const HybridEvent &ev : run.sim->events()) {
		for (const auto &a : ev.assignments) {
			if (a.first == "n") {
				nTimes.push_back(ev.time);
			}
		}
	}
	CHECK_EQ(nTimes.size(), (size_t)2);
	if (nTimes.size() == 2) {
		CHECK_NEAR(nTimes[0], 0.5, 1e-9);
		CHECK_NEAR(nTimes[1], 2.5, 1e-9);
	}
}

// change() of a Real discrete variable assigned by a time event, together with a state event.
static void testChangeOfDiscreteReal() {
	Run run;
	if (!simulate(run,
				  "model M\n"
				  "  discrete Real a(start = 1);\n"
				  "  Integer m(start = 0);\n"
				  "  Real x(start = 1);\n"
				  "equation\n"
				  "  der(x) = -x;\n"
				  "  when time >= 1 then\n"
				  "    a = 3;\n"
				  "  end when;\n"
				  "  when x < 0.2 then\n"
				  "    a = 3;\n"
				  "  end when;\n"
				  "  when change(a) then\n"
				  "    m = pre(m) + 1;\n"
				  "  end when;\n"
				  "end M;\n",
				  3.0)) {
		return;
	}
	// x < 0.2 at t = ln 5 = 1.609 assigns a = 3 again (no change); only the first assignment counts.
	CHECK_NEAR(finalValue(run, "m"), 1.0, 0.0);
}

// abs() in a continuous equation: y = 1 - t, x' = |y|, so x(2) = 1; the mode switches at t = 1.
static void testAbsInContinuousEquation() {
	Run run;
	if (!simulate(run,
				  "model M\n"
				  "  Real x(start = 0);\n"
				  "  Real y(start = 1);\n"
				  "equation\n"
				  "  der(x) = abs(y);\n"
				  "  der(y) = -1;\n"
				  "end M;\n",
				  2.0)) {
		return;
	}
	CHECK_EQ(run.sim->events().size(), (size_t)0); // state events change the mode, they are not when-events
	CHECK_NEAR(finalValue(run, "y"), -1.0, 1e-6);
	CHECK_NEAR(finalValue(run, "x"), 1.0, 1e-6);
	CHECK(run.sim->segmentCount() >= 2); // the mode switched at t = 1
}

// min and max: y' = min(t, 1) + max(-t, -2); y(3) = (0.5 + 2) - (2 + 2).
static void testMinMax() {
	Run run;
	if (!simulate(run,
				  "model M\n"
				  "  Real x(start = 0);\n"
				  "  Real y(start = 0);\n"
				  "equation\n"
				  "  der(x) = 1;\n"
				  "  der(y) = min(x, 1) + max(-x, -2);\n"
				  "end M;\n",
				  3.0)) {
		return;
	}
	// min(t,1): integral 0..3 = 0.5 + 2 = 2.5; max(-t,-2): -(integral of min(t,2)) = -(2 + 2) = -4.
	CHECK_NEAR(finalValue(run, "x"), 3.0, 1e-6);
	CHECK_NEAR(finalValue(run, "y"), 2.5 - 4.0, 1e-5);
}

// sign() is piecewise constant: x' = sign(y), y = 1 - t, so x rises until t = 1 and falls back to 0 at t = 2.
static void testSign() {
	Run run;
	if (!simulate(run,
				  "model M\n"
				  "  Real x(start = 0);\n"
				  "  Real y(start = 1);\n"
				  "equation\n"
				  "  der(x) = sign(y);\n"
				  "  der(y) = -1;\n"
				  "end M;\n",
				  2.0)) {
		return;
	}
	CHECK_NEAR(finalValue(run, "x"), 0.0, 1e-6);
}

// Nested calls and an expression argument: x' = abs(abs(y) - 0.5), y = 1 - t; the inner margin changes
// sign at t = 0.5 and t = 1.5.
static void testNestedPiecewise() {
	Run run;
	if (!simulate(run,
				  "model M\n"
				  "  Real x(start = 0);\n"
				  "  Real y(start = 1);\n"
				  "equation\n"
				  "  der(x) = abs(abs(y) - 0.5);\n"
				  "  der(y) = -1;\n"
				  "end M;\n",
				  2.0)) {
		return;
	}
	// The integral of | |1-t| - 0.5 | over [0, 2] is 0.25 on [0, 1] plus 0.25 on [1, 2].
	CHECK_NEAR(finalValue(run, "x"), 0.5, 1e-6);
}

// Pieces inside an if-equation branch and a when-body keep working.
static void testPiecewiseWithOtherConstructs() {
	Run run;
	if (!simulate(run,
				  "model M\n"
				  "  discrete Real a(start = 1);\n"
				  "  Real x(start = 1);\n"
				  "equation\n"
				  "  if a > 0 then\n"
				  "    der(x) = -abs(x);\n"
				  "  else\n"
				  "    der(x) = abs(x);\n"
				  "  end if;\n"
				  "  when time >= 1 then\n"
				  "    a = -abs(a);\n"
				  "  end when;\n"
				  "end M;\n",
				  2.0)) {
		return;
	}
	CHECK_NEAR(finalValue(run, "a"), -1.0, 0.0);
	CHECK_NEAR(finalValue(run, "x"), 1.0, 1e-4); // exp(-1) * exp(1)
}

// ---- Smooth approximations of abs/min/max/sign (no events) ----

static string absModel() {
	return "model M\n"
		   "  Real x(start = 0);\n"
		   "  Real y(start = 1);\n"
		   "equation\n"
		   "  der(x) = abs(y);\n"
		   "  der(y) = -1;\n"
		   "end M;\n";
}

// x' = abs(y), y = 1 - t: the exact x(2) is 1; the smooth version is off by about eps^2 ln(1/eps).
static double smoothAbsError(double eps) {
	Run run;
	if (!simulate(run, absModel(), 2.0, PW_SMOOTH, eps)) {
		return NAN;
	}
	CHECK_EQ(run.sim->events().size(), (size_t)0);
	return std::fabs(finalValue(run, "x") - 1.0);
}

static void testSmoothAbsConverges() {
	double coarse = smoothAbsError(1e-2);
	double fine = smoothAbsError(1e-3);
	CHECK(coarse < 2e-3);
	CHECK(fine < 2e-4);
	CHECK(fine < coarse); // the error shrinks with eps
}

static void testSmoothMinMaxSign() {
	{
		Run run;
		if (simulate(run,
					 "model M\n"
					 "  Real x(start = 0);\n"
					 "  Real y(start = 0);\n"
					 "equation\n"
					 "  der(x) = 1;\n"
					 "  der(y) = min(x, 1) + max(-x, -2);\n"
					 "end M;\n",
					 3.0, PW_SMOOTH, 1e-3)) {
			CHECK_NEAR(finalValue(run, "y"), 2.5 - 4.0, 2e-3);
		}
	}
	{
		Run run;
		if (simulate(run,
					 "model M\n"
					 "  Real x(start = 0);\n"
					 "  Real y(start = 1);\n"
					 "equation\n"
					 "  der(x) = sign(y);\n"
					 "  der(y) = -1;\n"
					 "end M;\n",
					 2.0, PW_SMOOTH, 1e-3)) {
			CHECK_NEAR(finalValue(run, "x"), 0.0, 5e-3);
		}
	}
}

// The model that chatters with events (sliding friction) runs through with smoothing and settles
// where sign(v) = 0.5, i.e. v = eps / sqrt(3).
static void testSmoothSlidingFriction() {
	Run run;
	const double eps = 1e-2;
	if (!simulate(run,
				  "model M\n"
				  "  Real v(start = 0.5);\n"
				  "equation\n"
				  "  der(v) = 0.5 - sign(v);\n"
				  "end M;\n",
				  2.0, PW_SMOOTH, eps)) {
		return;
	}
	CHECK_NEAR(finalValue(run, "v"), eps / std::sqrt(3.0), 1e-3);
}

// Equations listed in a different order per mode.
static void testReorderedEquations() {
	Run run;
	if (!simulate(run,
				  "model M\n"
				  "  Boolean m(start = true);\n"
				  "  Real x(start = 1);\n"
				  "  Real y(start = 2);\n"
				  "equation\n"
				  "  if m then\n"
				  "    der(x) = -x;\n"
				  "    der(y) = -y;\n"
				  "  else\n"
				  "    der(y) = -y;\n"
				  "    der(x) = -2*x;\n"
				  "  end if;\n"
				  "  when time >= 1 then\n"
				  "    m = false;\n"
				  "  end when;\n"
				  "end M;\n",
				  2.0)) {
		return;
	}
	CHECK_NEAR(finalValue(run, "x"), std::exp(-3.0), 1e-4);
	CHECK_NEAR(finalValue(run, "y"), 2.0 * std::exp(-2.0), 1e-4);
}

int main() {
	testSmoothAbsConverges();
	testSmoothMinMaxSign();
	testSmoothSlidingFriction();
	testSignChatterIsReported();
	testAbsInContinuousEquation();
	testMinMax();
	testSign();
	testNestedPiecewise();
	testPiecewiseWithOtherConstructs();
	testPulseBetweenSamples();
	testDipBetweenSamples();
	testPulseInIfEquation();
	testEdgeAndChange();
	testChangeOfDiscreteReal();
	testStateAlgebraicSwap();
	testReorderedEquations();
	testJsonResult();
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
