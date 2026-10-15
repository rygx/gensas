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

int main() {
	testDiscreteGainSwitch();
	testReinitAtSamples();
	testIfEquationSwitch();
	testIfExpressionWithCounter();
	testAlgebraicVariable();
	testNoEvents();
	return sastest::finish("hybrid_sim_test");
}
