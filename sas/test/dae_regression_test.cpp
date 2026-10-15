// Regression test for the solver with an algebraic equation: der(x) = -y, y = k*x has the solution
// x = x0 * exp(-k t). Guards against evaluating the Pade solution with stale cached coefficients
// during the initial AE calibration.
#include "sas/SasComputation.h"
#include "sas/SasExpr.h"
#include "sas/SasLexico.h"
#include "sas/test/test_util.h"

#include <cmath>
#include <memory>

using namespace che::core;

static const char *kModel =
	"model Dae\n"
	"  parameter Real k = 2;\n"
	"  Real x(start = 1);\n"
	"  Real y(start = 1);\n"
	"equation\n"
	"  y = k*x;\n"
	"  der(x) = -y;\n"
	"end Dae;\n";

// y(0) is inconsistent with y = k*x, so the solver has to calibrate it first.
static const char *kInconsistent =
	"model Dae\n"
	"  parameter Real k = 2;\n"
	"  Real x(start = 1);\n"
	"  Real y(start = 5);\n"
	"equation\n"
	"  y = k*x;\n"
	"  der(x) = -y;\n"
	"end Dae;\n";

static void solveAndCheck(const char *src) {
	SasModelParser parser;
	parser.parseString(src);
	std::shared_ptr<SasModel> model = GlobalPool::getInstance().createNewTopModel();
	parser.currentModel = model;
	parseSingleSasModel(&parser);
	CHECK_EQ(parser.errCount, 0);

	SasComputationModel compModel;
	compModel.generateCompModel(model);
	compModel.preprocessModels();
	compModel.generateDAEs();
	CHECK_EQ(compModel.nX, 1);
	CHECK_EQ(compModel.nY, 1);
	CHECK_EQ(compModel.nDE, 1);
	CHECK_EQ(compModel.nAE, 1);

	SasComputationOptions options;
	options.endTime = 2.0;
	SasSolutionSet *sol = compModel.solve(options);
	CHECK(sol != nullptr);
	if (sol == nullptr || sol->solutionLink.empty()) {
		delete sol;
		return;
	}
	CHECK_NEAR(sol->solutionLink.back()->absEnd, options.endTime, 1e-3);
	std::shared_ptr<SasSolution> last = sol->solutionLink.back();
	arma::vec end = last->solution->getSolValue(last->absEnd - last->absStart);
	CHECK_NEAR(end(0), std::exp(-4.0), 1e-4);
	CHECK_NEAR(end(1), 2.0 * std::exp(-4.0), 1e-4);
	delete sol;
}

int main() {
	solveAndCheck(kModel);
	solveAndCheck(kInconsistent);
	return sastest::finish("dae_regression_test");
}
