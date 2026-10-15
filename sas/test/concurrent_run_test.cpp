// Independent simulations in one process must not share mutable state: several threads each parse,
// check and solve their own model at the same time, and every result must equal the closed-form
// solution (and the result of the same model run alone).
#include "sas/HybridSimulator.h"
#include "sas/ModelChecker.h"
#include "sas/SasComputation.h"
#include "sas/SasExpr.h"
#include "sas/SasLexico.h"
#include "sas/test/test_util.h"
#include "util/SafeArmadillo.h"

#include <cmath>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace che::core;
using std::string;

namespace {

	const int kThreads = 8;
	const int kRounds = 4;

	SolverFeatures hybridFeatures() {
		SolverFeatures f;
		f.discreteVariables = f.ifExpressions = f.ifEquations = f.whenEquations = true;
		f.piecewiseFunctions = f.eventFunctions = f.simulationTime = true;
		f.timeEventsOnly = true;
		f.stateEvents = true;
		return f;
	}

	SasComputationOptions makeOptions(double endTime) {
		SasComputationOptions o;
		o.endTime = endTime;
		o.segment = 1.0;
		o.alphaTol = 1e-4;
		o.errorTol = 1e-6;
		o.nLvl = 15;
		return o;
	}

	// Plain ODE x' = -k x, x(0) = 1 through the non-hybrid path. Returns x(2), or NAN on failure.
	double runOde(double k) {
		std::ostringstream src;
		src << "model Decay\n  parameter Real k = " << k << ";\n  Real x(start = 1);\nequation\n  der(x) = -k*x;\nend Decay;\n";
		SasModelParser parser;
		parser.parseString(src.str().c_str());
		parser.currentModel = parser.pool.createNewTopModel();
		parseSingleSasModel(&parser);
		if (parser.errCount != 0) {
			return NAN;
		}
		SasComputationModel comp;
		comp.generateCompModel(parser.pool);
		comp.preprocessModels();
		comp.generateDAEs();
		SasComputationOptions options = makeOptions(2.0);
		SasSolutionSet *sol = comp.solve(options);
		double v = NAN;
		if (sol != NULL) {
			auto found = sol->findSasSolution(2.0);
			if (!found.empty()) {
				std::shared_ptr<SasSolution> seg = found.back();
				v = seg->solution->getSolValue(2.0 - seg->absStart)(0);
			}
		}
		delete sol;
		return v;
	}

	// Hybrid model with a discrete gain switch at t = 1: x(2) = exp(-1) * exp(-(1+d)) with d = gain.
	double runHybrid(double gain) {
		std::ostringstream src;
		src << "model M\n  discrete Real a(start = 1);\n  Real x(start = 1);\nequation\n  der(x) = -a*x;\n"
			<< "  when time >= 1 then\n    a = " << gain << ";\n  end when;\nend M;\n";
		SasModelParser parser;
		parser.parseString(src.str().c_str());
		std::shared_ptr<SasModel> model = parser.pool.createNewTopModel();
		parser.currentModel = model;
		parseSingleSasModel(&parser);
		if (parser.errCount != 0) {
			return NAN;
		}
		rewritePiecewiseFunctions(*model);
		std::ostringstream os;
		if (checkModel(*model, hybridFeatures(), "test.mo", os) != 0) {
			return NAN;
		}
		HybridSimulator sim(model);
		SasSolutionSet *sol = sim.run(makeOptions(2.0));
		double v = NAN;
		if (sol != NULL) {
			std::map<string, double>::const_iterator it = sim.finalValues().find("x");
			if (it != sim.finalValues().end()) {
				v = it->second;
			}
		}
		delete sol;
		return v;
	}

} // namespace

int main() {
	std::vector<double> odeOut(kThreads * kRounds, NAN), hybOut(kThreads * kRounds, NAN);
	std::vector<std::thread> threads;
	for (int t = 0; t < kThreads; t++) {
		threads.emplace_back([&, t]() {
			for (int r = 0; r < kRounds; r++) {
				odeOut[t * kRounds + r] = runOde(0.5 + 0.25 * t);
				hybOut[t * kRounds + r] = runHybrid(2.0 + t);
			}
		});
	}
	for (std::thread &th : threads) {
		th.join();
	}
	for (int t = 0; t < kThreads; t++) {
		for (int r = 0; r < kRounds; r++) {
			CHECK_NEAR(odeOut[t * kRounds + r], std::exp(-2.0 * (0.5 + 0.25 * t)), 1e-4);
			CHECK_NEAR(hybOut[t * kRounds + r], std::exp(-1.0) * std::exp(-(2.0 + t)), 1e-4);
		}
	}
	return sastest::finish("concurrent_run_test");
}
