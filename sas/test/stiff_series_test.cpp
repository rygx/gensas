// Stiff problems: the series of a fast decaying variable is computed backwards in the order
// (SasComputationOptions::stableSeries), so that the step is not limited to about 1/|d|.
//   - a linear fast mode y' = -1000 (y - c) coupled to an oscillator: exact solution by expm();
//   - Robertson's problem: the same result with and without the option, in far fewer steps.
#include "util/SafeArmadillo.h"
#include "sas/SasComputation.h"
#include "sas/SasExpr.h"
#include "sas/SasLexico.h"
#include "sas/test/test_util.h"

#include <cmath>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

using namespace che::core;

static const char *kLinear =
	"model L\n"
	"  Real y(start = 1);\n"
	"  Real c(start = 1);\n"
	"  Real s(start = 0);\n"
	"equation\n"
	"  der(c) = -s;\n"
	"  der(s) = c;\n"
	"  der(y) = -1000 * (y - c);\n"
	"end L;\n";

// The same, but the fast variable starts far from the manifold y = c (a transient): the series of the manifold must
// not be used.
static const char *kLinearOff =
	"model LO\n"
	"  Real y(start = 5);\n"
	"  Real c(start = 1);\n"
	"  Real s(start = 0);\n"
	"equation\n"
	"  der(c) = -s;\n"
	"  der(s) = c;\n"
	"  der(y) = -1000 * (y - c);\n"
	"end LO;\n";

static const char *kRobertson =
	"model R\n"
	"  Real y1(start = 1);\n"
	"  Real y2(start = 0);\n"
	"  Real y3(start = 0);\n"
	"equation\n"
	"  der(y1) = -0.04 * y1 + 1e4 * y2 * y3;\n"
	"  der(y2) = 0.04 * y1 - 1e4 * y2 * y3 - 3e7 * y2 * y2;\n"
	"  der(y3) = 3e7 * y2 * y2;\n"
	"end R;\n";

// Robertson's problem with the third species as an algebraic variable (a DAE, so the effective self-derivative of a
// fast variable is formed with the algebraic variables eliminated).
static const char *kRobertsonDae =
	"model RD\n"
	"  Real y1(start = 1);\n"
	"  Real y2(start = 0);\n"
	"  Real y3(start = 0);\n"
	"equation\n"
	"  der(y1) = -0.04 * y1 + 1e4 * y2 * y3;\n"
	"  der(y2) = 0.04 * y1 - 1e4 * y2 * y3 - 3e7 * y2 * y2;\n"
	"  0 = y1 + y2 + y3 - 1;\n"
	"end RD;\n";

struct Result {
	size_t steps = 0;
	int manifoldSteps = 0; // steps whose series was computed backwards for a fast variable
	std::vector<double> x; // final state in the order of `names`
};

static Result run(const char *model, const std::vector<std::string> &names, double endTime, double errorTol, bool stable, bool adaptive = true) {
	Result res;
	SasModelParser parser;
	parser.parseString(model);
	parser.currentModel = parser.pool.createNewTopModel();
	parseSingleSasModel(&parser);
	CHECK_EQ(parser.errCount, 0);
	SasComputationModel compModel;
	compModel.generateCompModel(parser.pool);
	compModel.preprocessModels();
	compModel.generateDAEs();
	SasComputationOptions options;
	options.endTime = endTime;
	options.segment = endTime;
	options.alphaTol = 1e-4;
	options.errorTol = errorTol;
	options.nLvl = 30;
	options.stableSeries = stable;
	options.stableSeriesAdaptive = adaptive;
	std::unique_ptr<SasSolutionSet> sol(compModel.solve(options));
	CHECK(sol != nullptr);
	if (sol == nullptr || sol->solutionLink.empty()) {
		return res;
	}
	res.steps = sol->solutionLink.size();
	res.manifoldSteps = compModel.manifoldSeriesSteps();
	SasComputationModel *cm = sol->solutionLink.back()->linkCompModel;
	auto found = sol->findSasSolution(endTime);
	CHECK(!found.empty());
	if (found.empty()) {
		return res;
	}
	std::shared_ptr<SasSolution> st = found.back();
	arma::vec x = st->solution->getSolValue(endTime - st->absStart);
	for (const std::string &n : names) {
		int idx = -1;
		for (size_t i = 0; i < cm->idTable.size(); i++) {
			if (n == (const char *)cm->idTable[i]->str) {
				idx = (int)i;
			}
		}
		CHECK(idx >= 0 && idx < (int)x.n_rows);
		res.x.push_back(idx >= 0 && idx < (int)x.n_rows ? x(idx) : NAN);
	}
	return res;
}

int main() {
	// a start off the manifold (y(0) = 5, c = cos t, s = sin t): the transient of y must be computed, inside the boundary
	// layer and after it, with and without the stable series; y = yp(t) + (5 - yp(0)) exp(-1000 t)
	for (double T : {0.002, 10.0}) {
		auto yExact = [](double t) {
			auto yp = [](double u) { return (1e6 * std::cos(u) + 1000.0 * std::sin(u)) / (1e6 + 1.0); };
			return yp(t) + (5.0 - yp(0.0)) * std::exp(-1000.0 * t);
		};
		Result on = run(kLinearOff, {"y", "c", "s"}, T, 1e-6, true);
		Result off = run(kLinearOff, {"y", "c", "s"}, T, 1e-6, false);
		CHECK_EQ((int)on.x.size(), 3);
		CHECK_EQ((int)off.x.size(), 3);
		if (on.x.size() == 3 && off.x.size() == 3) {
			CHECK_NEAR(on.x[0], yExact(T), 1e-5);
			CHECK_NEAR(off.x[0], yExact(T), 1e-5);
			CHECK_NEAR(on.x[1], std::cos(T), 1e-5);
			CHECK_NEAR(on.x[2], std::sin(T), 1e-5);
			CHECK_NEAR(off.x[1], std::cos(T), 1e-5);
			CHECK_NEAR(off.x[2], std::sin(T), 1e-5);
		}
	}

	// linear: exact solution of x' = A x
	{
		arma::mat A = {{-1000.0, 1000.0, 0.0}, {0.0, 0.0, -1.0}, {0.0, 1.0, 0.0}}; // y, c, s
		arma::vec x0 = {1.0, 1.0, 0.0};
		const double T = 10.0;
		arma::vec ref = arma::expmat(A * T) * x0;
		Result on = run(kLinear, {"y", "c", "s"}, T, 1e-6, true);
		Result off = run(kLinear, {"y", "c", "s"}, T, 1e-6, false);
		CHECK_EQ((int)on.x.size(), 3);
		CHECK_EQ((int)off.x.size(), 3);
		for (int k = 0; k < 3 && (int)on.x.size() == 3 && (int)off.x.size() == 3; k++) {
			CHECK_NEAR(on.x[k], ref(k), 1e-4);
			CHECK_NEAR(off.x[k], ref(k), 1e-4);
		}
		// the step with the option is not limited by the fast mode
		CHECK(on.steps > 0 && on.steps <= 10);
		CHECK(on.steps * 5 < off.steps);
	}
	// Robertson: same result, far fewer steps
	{
		const double T = 100.0;
		Result on = run(kRobertson, {"y1", "y2", "y3"}, T, 1e-6, true);
		Result off = run(kRobertson, {"y1", "y2", "y3"}, T, 1e-6, false);
		CHECK_EQ((int)on.x.size(), 3);
		CHECK_EQ((int)off.x.size(), 3);
		for (int k = 0; k < 3 && (int)on.x.size() == 3 && (int)off.x.size() == 3; k++) {
			CHECK_NEAR(on.x[k], off.x[k], 1e-4);
		}
		CHECK(on.steps > 0 && on.steps * 10 < off.steps);
		// the invariant y1 + y2 + y3 = 1
		if (on.x.size() == 3) {
			CHECK_NEAR(on.x[0] + on.x[1] + on.x[2], 1.0, 1e-6);
		}
	}
	// Robertson as a DAE: the same result in far fewer steps; y1(100) = 0.61723488 (reference solution)
	{
		const double T = 100.0;
		Result on = run(kRobertsonDae, {"y1", "y2", "y3"}, T, 1e-6, true);
		Result off = run(kRobertsonDae, {"y1", "y2", "y3"}, T, 1e-6, false);
		CHECK_EQ((int)on.x.size(), 3);
		CHECK_EQ((int)off.x.size(), 3);
		if (on.x.size() == 3 && off.x.size() == 3) {
			CHECK_NEAR(on.x[0], 0.617234882, 1e-5);
			CHECK_NEAR(on.x[0] + on.x[1] + on.x[2], 1.0, 1e-9);
			CHECK_NEAR(on.x[0], off.x[0], 1e-4);
		}
		CHECK(on.steps > 0 && on.steps * 10 < off.steps);
	}
	// The manifold series is tried only where it pays. A fast mode (d = -2000) under a forcing at 300 rad/s: the
	// step is limited by the forcing, the manifold series saves 1.75 times the steps and costs several passes of
	// the series per step, so the adaptive rule keeps using the forward series almost always, and the result is
	// the same as with the series at every step.
	{
		const char *forced =
			"model F\n"
			"  Real y(start = 1);\n"
			"  Real c(start = 1);\n"
			"  Real s(start = 0);\n"
			"equation\n"
			"  der(c) = -300 * s;\n"
			"  der(s) = 300 * c;\n"
			"  der(y) = -2000 * (y - c);\n"
			"end F;\n";
		const double T = 1.0;
		Result adaptive = run(forced, {"y", "c", "s"}, T, 1e-3, true, true);
		Result always = run(forced, {"y", "c", "s"}, T, 1e-3, true, false);
		Result forward = run(forced, {"y", "c", "s"}, T, 1e-3, false, false);
		CHECK_EQ((int)adaptive.x.size(), 3);
		CHECK_EQ((int)always.x.size(), 3);
		CHECK_EQ((int)forward.x.size(), 3);
		for (int k = 0; k < 3 && adaptive.x.size() == 3 && always.x.size() == 3 && forward.x.size() == 3; k++) {
			CHECK_NEAR(adaptive.x[k], always.x[k], 2e-3);
			CHECK_NEAR(adaptive.x[k], forward.x[k], 2e-3);
		}
		CHECK(always.manifoldSteps > 20);
		CHECK(adaptive.manifoldSteps * 3 < always.manifoldSteps);
		CHECK(adaptive.steps <= forward.steps);
	}
	// A linear fast mode on the manifold: the adaptive rule finds that the manifold series pays and keeps it (a few
	// steps, as with the series at every step, plus at most one trial)
	{
		const double T = 10.0;
		Result adaptive = run(kLinear, {"y", "c", "s"}, T, 1e-6, true, true);
		Result always = run(kLinear, {"y", "c", "s"}, T, 1e-6, true, false);
		CHECK(adaptive.steps <= always.steps + 2);
		CHECK(adaptive.manifoldSteps >= 1);
	}
	return sastest::finish("stiff_series_test");
}
