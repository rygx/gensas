// End-to-end regression test of the Modelica path (parse -> DAE -> SAS solve)
// on resources/mofile/test_solve_ode.mo, a 3-state linear ODE x' = A x.
// The exact solution x(t) = expm(A t) x0 is computed with Armadillo and the SAS
// solution is compared against it. One model per process (GlobalPool is a singleton).
#include "util/SafeArmadillo.h"
#include "sas/SasComputation.h"
#include "sas/SasExpr.h"
#include "sas/SasLexico.h"
#include "sas/test/test_util.h"

#include <cstring>
#include <iostream>
#include <memory>
#include <string>

using namespace che::core;
using namespace arma;
using std::shared_ptr;

static const char *kModel =
	"model exODE01\n"
	"  Real x1 \"x1\";\n"
	"  Real x2 \"x2\";\n"
	"  Real x3 \"x3\";\n"
	"  parameter Real k1=-1.0 \"k1\";\n"
	"  parameter Real k2=-2.0 \"k2\";\n"
	"  parameter Real k3=-3.0 \"k3\";\n"
	"  parameter Real yk0=0.3 \"yk0\";\n"
	"  parameter Real yk1=0.4 \"yk1\";\n"
	"  parameter Real yk2=0.2 \"yk2\";\n"
	"  parameter Real yk3=0.1 \"yk3\";\n"
	"initial equation\n"
	"  x1=3.0;\n"
	"  x2=3.0;\n"
	"  x3=3.0;\n"
	"equation\n"
	"  der(x1)=k1*x1+yk0*x2+yk1*x3;\n"
	"  der(x2)=k2*x2+yk0*x3+yk2*x1;\n"
	"  der(x3)=k3*x3+yk0*x1+yk3*x2;\n"
	"end exODE01;\n";

int main() {
	SasModelParser parser;
	parser.parseString(kModel);
	parser.currentModel = GlobalPool::getInstance().createNewTopModel();
	parseSingleSasModel(&parser);
	CHECK_EQ(parser.errCount, 0);

	SasComputationModel compModel;
	compModel.generateCompModel(GlobalPool::getInstance());
	compModel.preprocessModels();
	compModel.generateDAEs();
	CHECK_EQ(compModel.nDE, 3);
	CHECK_EQ(compModel.nAE, 0);
	CHECK_EQ(compModel.nX, 3);
	CHECK_EQ(compModel.nY, 0);

	SasComputationOptions options;
	options.endTime = 5.0;
	options.segment = 1.0;
	options.alphaTol = 1e-4;
	options.errorTol = 1e-6;
	options.nLvl = 15;
	SasSolutionSet *sol = compModel.solve(options);
	CHECK(sol != nullptr);
	if (sol == nullptr || sol->solutionLink.empty()) {
		return sastest::finish("ode_regression_test");
	}
	CHECK_NEAR(sol->solutionLink.back()->absEnd, options.endTime, 1e-3);

	// State order in the solution vector: the order of the der() equations.
	SasComputationModel *cm = sol->solutionLink.back()->linkCompModel;
	int idx[3] = {-1, -1, -1};
	const char *names[3] = {"x1", "x2", "x3"};
	for (size_t i = 0; i < cm->idTable.size(); i++) {
		for (int k = 0; k < 3; k++) {
			if (std::strcmp((const char *)cm->idTable[i]->str, names[k]) == 0) {
				idx[k] = (int)i;
			}
		}
	}
	for (int k = 0; k < 3; k++) {
		CHECK(idx[k] >= 0);
	}

	// x' = A x with A[i][j] read from the model above.
	mat A = {{-1.0, 0.3, 0.4},
			 {0.2, -2.0, 0.3},
			 {0.3, 0.1, -3.0}};
	vec x0 = {3.0, 3.0, 3.0};
	for (double t : {0.0, 0.25, 0.5, 1.0, 2.0, 3.5, 5.0}) {
		auto found = sol->findSasSolution(t);
		CHECK(!found.empty());
		if (found.empty()) {
			continue;
		}
		std::shared_ptr<SasSolution> s = found.back();
		vec x = s->solution->getSolValue(t - s->absStart);
		vec ref = expmat(A * t) * x0;
		for (int k = 0; k < 3; k++) {
			if (idx[k] >= 0) {
				CHECK_NEAR(x(idx[k]), ref(k), 1e-4);
			}
		}
	}

	delete sol;
	return sastest::finish("ode_regression_test");
}
