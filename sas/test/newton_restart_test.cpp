// The algebraic variables at the start of a run (or after an event) are found by Newton's method from the given
// values. When these are far from the solution the plain iteration diverges; the damped iteration finds the
// solution. Model: x' = -x + y with 0 = atan(y) - 0.1 and the start value y = 3 (the plain Newton step from
// there goes to y = -8.5 and away).
#include "util/SafeArmadillo.h"
#include "sas/SasComputation.h"
#include "sas/SasExpr.h"
#include "sas/SasLexico.h"
#include "sas/test/test_util.h"

#include <cmath>
#include <cstring>
#include <memory>

using namespace che::core;

static const char *kModel =
	"model NR\n"
	"  Real x(start = 1);\n"
	"  Real y(start = 3);\n"
	"equation\n"
	"  der(x) = -x + y;\n"
	"  0 = atan(y) - 0.1;\n"
	"end NR;\n";

int main() {
	SasModelParser parser;
	parser.parseString(kModel);
	parser.currentModel = parser.pool.createNewTopModel();
	parseSingleSasModel(&parser);
	CHECK_EQ(parser.errCount, 0);
	SasComputationModel compModel;
	compModel.generateCompModel(parser.pool);
	compModel.preprocessModels();
	compModel.generateDAEs();
	SasComputationOptions options;
	options.endTime = 2.0;
	options.segment = 2.0;
	options.alphaTol = 1e-4;
	options.errorTol = 1e-8;
	options.nLvl = 20;
	std::unique_ptr<SasSolutionSet> sol(compModel.solve(options));
	CHECK(sol != nullptr);
	if (sol == nullptr || sol->solutionLink.empty()) {
		return sastest::finish("newton_restart_test");
	}
	SasComputationModel *cm = sol->solutionLink.back()->linkCompModel;
	int ix = -1, iy = -1;
	for (size_t i = 0; i < cm->idTable.size(); i++) {
		if (std::strcmp((const char *)cm->idTable[i]->str, "x") == 0)
			ix = (int)i;
		if (std::strcmp((const char *)cm->idTable[i]->str, "y") == 0)
			iy = (int)i;
	}
	CHECK(ix >= 0 && iy >= 0);
	auto found = sol->findSasSolution(2.0);
	CHECK(!found.empty());
	if (found.empty() || ix < 0 || iy < 0) {
		return sastest::finish("newton_restart_test");
	}
	std::shared_ptr<SasSolution> st = found.back();
	arma::vec v = st->solution->getSolValue(2.0 - st->absStart);
	CHECK(ix < (int)v.n_rows && iy < (int)v.n_rows);
	if (ix >= (int)v.n_rows || iy >= (int)v.n_rows) {
		return sastest::finish("newton_restart_test");
	}
	const double ys = std::tan(0.1);
	CHECK_NEAR(v(iy), ys, 1e-8);
	// x' = -x + ys, x(0) = 1
	CHECK_NEAR(v(ix), ys + (1.0 - ys) * std::exp(-2.0), 1e-6);
	return sastest::finish("newton_restart_test");
}
