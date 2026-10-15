// MAT v5 result writer, read back with matio.
#include "sas/HybridSimulator.h"
#include "sas/ModelChecker.h"
#include "sas/SasExpr.h"
#include "sas/SasLexico.h"
#include "sas/SasResult.h"
#include "sas/SasResultMat.h"
#include "sas/test/test_util.h"

#include "matio.h"

#include <cmath>
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

static const char *kModel =
	"model M\n"
	"  discrete Real a(start = 1);\n"
	"  Real x(start = 1);\n"
	"equation\n"
	"  der(x) = -a*x;\n"
	"  when time >= 1 then\n"
	"    a = 3;\n"
	"  end when;\n"
	"end M;\n";

static std::vector<double> readDoubles(mat_t *fp, const char *name, size_t &rows, size_t &cols) {
	std::vector<double> out;
	rows = cols = 0;
	matvar_t *v = Mat_VarRead(fp, name);
	CHECK(v != NULL);
	if (v == NULL) {
		return out;
	}
	CHECK(v->class_type == MAT_C_DOUBLE);
	CHECK_EQ(v->rank, 2);
	rows = v->dims[0];
	cols = v->dims[1];
	const double *data = (const double *)v->data;
	out.assign(data, data + rows * cols);
	Mat_VarFree(v);
	return out;
}

// Row r of a char matrix, trailing blanks removed.
static std::string charRow(mat_t *fp, const char *name, size_t r) {
	matvar_t *v = Mat_VarRead(fp, name);
	CHECK(v != NULL);
	if (v == NULL) {
		return "";
	}
	CHECK(v->class_type == MAT_C_CHAR);
	size_t rows = v->dims[0];
	size_t cols = v->dims[1];
	std::string s;
	for (size_t c = 0; c < cols && r < rows; c++) {
		s += ((const char *)v->data)[c * rows + r];
	}
	Mat_VarFree(v);
	while (!s.empty() && s.back() == ' ') {
		s.pop_back();
	}
	return s;
}

static void testMatFile() {
	Run run;
	if (!simulate(run, kModel, 2.0)) {
		return;
	}
	SasResult r = buildResult(*run.sol, 0.5);
	std::string err;
	CHECK(makeMatResultWriter()->write(r, "result_test.mat", err));
	mat_t *fp = Mat_Open("result_test.mat", MAT_ACC_RDONLY);
	CHECK(fp != NULL);
	if (fp == NULL) {
		return;
	}
	CHECK_EQ(Mat_GetVersion(fp), (int)MAT_FT_MAT5);
	size_t rows, cols;
	std::vector<double> v = readDoubles(fp, "gensas_version", rows, cols);
	CHECK_EQ(v.size(), (size_t)1);
	if (v.size() == 1) {
		CHECK_NEAR(v[0], kResultFormatVersion, 0.0);
	}
	std::vector<double> t = readDoubles(fp, "t", rows, cols);
	CHECK_EQ(rows, (size_t)1);
	CHECK_EQ(t.size(), r.time.size());
	std::vector<double> s = readDoubles(fp, "s", rows, cols);
	CHECK_EQ(rows, (size_t)1);
	CHECK_EQ(cols, r.time.size());
	for (size_t j = 0; j < s.size() && j < r.time.size(); j++) {
		CHECK_NEAR(s[j], r.values[0][j], 0.0);
		CHECK_NEAR(t[j], r.time[j], 0.0);
	}
	std::vector<double> d = readDoubles(fp, "d", rows, cols);
	CHECK_EQ(rows, (size_t)1);
	CHECK_EQ(cols, r.time.size());
	for (size_t j = 0; j < d.size() && j < r.time.size(); j++) {
		CHECK_NEAR(d[j], r.discrete[0][j], 0.0);
	}
	std::vector<double> et = readDoubles(fp, "event_t", rows, cols);
	CHECK_EQ(et.size(), (size_t)1);
	if (et.size() == 1) {
		CHECK_NEAR(et[0], 1.0, 1e-12);
	}
	CHECK(charRow(fp, "names", 0) == "x");
	CHECK(charRow(fp, "dnames", 0) == "a");
	Mat_Close(fp);

	// Empty results are an error, not an empty file.
	SasResult empty;
	CHECK(!makeMatResultWriter()->write(empty, "result_test_empty.mat", err));
	CHECK(!err.empty());
}

int main() {
	testMatFile();
	return sastest::finish("result_mat_test");
}
