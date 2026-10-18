// Regression test: an undamped oscillator (c' = -s, s' = c) solved over a long span with high orders. The
// residual used to choose the step must be the defect of the Pade approximant that is evaluated; a residual
// built from a separately constructed derivative approximant is blind to a spurious pole and let steps through
// whose result was far off (s = 549 instead of -0.54 at order 30).
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
	"model osc\n"
	"  Real c;\n"
	"  Real s;\n"
	"initial equation\n"
	"  c=1.0;\n"
	"  s=0.0;\n"
	"equation\n"
	"  der(c)=-s;\n"
	"  der(s)=c;\n"
	"end osc;\n";

static void run(int nLvl, double errorTol) {
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
	options.endTime = 10.0;
	options.segment = 10.0;
	options.alphaTol = 1e-4;
	options.errorTol = errorTol;
	options.nLvl = nLvl;
	std::unique_ptr<SasSolutionSet> sol(compModel.solve(options));
	CHECK(sol != nullptr);
	if (sol == nullptr || sol->solutionLink.empty()) {
		return;
	}
	SasComputationModel *cm = sol->solutionLink.back()->linkCompModel;
	int ic = -1, is = -1;
	for (size_t i = 0; i < cm->idTable.size(); i++) {
		if (std::strcmp((const char *)cm->idTable[i]->str, "c") == 0)
			ic = (int)i;
		if (std::strcmp((const char *)cm->idTable[i]->str, "s") == 0)
			is = (int)i;
	}
	CHECK(ic >= 0 && is >= 0);
	auto found = sol->findSasSolution(10.0);
	CHECK(!found.empty());
	if (found.empty() || ic < 0 || is < 0) {
		return;
	}
	std::shared_ptr<SasSolution> st = found.back();
	arma::vec x = st->solution->getSolValue(10.0 - st->absStart);
	CHECK(ic < (int)x.n_rows && is < (int)x.n_rows);
	if (ic >= (int)x.n_rows || is >= (int)x.n_rows) {
		return;
	}
	CHECK_NEAR(x(ic), std::cos(10.0), 1e-3);
	CHECK_NEAR(x(is), std::sin(10.0), 1e-3);
}

int main() {
	for (int nLvl : {15, 30, 40}) {
		for (double tol : {1e-5, 1e-8}) {
			run(nLvl, tol);
		}
	}
	return sastest::finish("oscillator_regression_test");
}
