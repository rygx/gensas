// Result sampling and result writers (JSON, CSV).
#include "sas/HybridSimulator.h"
#include "sas/ModelChecker.h"
#include "sas/SasExpr.h"
#include "sas/SasLexico.h"
#include "sas/SasResult.h"
#include "sas/test/test_util.h"

#include "json/json.h"

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
	shared_ptr<SasModel> model = parser.pool.createNewTopModel();
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

static std::vector<std::string> readLines(const char *fileName) {
	std::vector<std::string> lines;
	std::ifstream in(fileName);
	std::string line;
	while (std::getline(in, line)) {
		lines.push_back(line);
	}
	return lines;
}

static void testBuildResult() {
	Run run;
	if (!simulate(run, kModel, 2.0)) {
		return;
	}
	SasResult r = buildResult(*run.sol, 0.5);
	CHECK(!r.empty());
	CHECK_EQ(r.names.size(), (size_t)1);
	CHECK_EQ(r.discreteNames.size(), (size_t)1);
	CHECK_EQ(r.values.size(), (size_t)1);
	CHECK_EQ(r.discrete.size(), (size_t)1);
	CHECK_EQ(r.values[0].size(), r.time.size());
	CHECK_EQ(r.discrete[0].size(), r.time.size());
	CHECK_EQ(r.events.size(), (size_t)1);
	CHECK_NEAR(r.startTime, 0.0, 0.0);
	CHECK_NEAR(r.endTime, 2.0, 1e-9);
	if (r.names.size() == 1) {
		CHECK(r.names[0] == "x");
	}
	// Non-decreasing times; x follows the closed-form solution.
	for (size_t j = 0; j < r.time.size(); j++) {
		if (j > 0) {
			CHECK(r.time[j] >= r.time[j - 1]);
		}
		double t = r.time[j];
		double expected = t < 1.0 ? std::exp(-t) : std::exp(-1.0) * std::exp(-3.0 * (t - 1.0));
		CHECK_NEAR(r.values[0][j], expected, 1e-4);
	}
	// A non-positive interval gives an empty result.
	CHECK(buildResult(*run.sol, 0.0).empty());
}

static void testJsonVersion() {
	Run run;
	if (!simulate(run, kModel, 2.0)) {
		return;
	}
	SasResult r = buildResult(*run.sol, 0.5);
	std::string err;
	CHECK(makeResultWriter("json")->write(r, "result_test.json", err));
	std::ifstream in("result_test.json");
	Json::Value root;
	Json::CharReaderBuilder builder;
	std::string errs;
	CHECK(Json::parseFromStream(builder, in, &root, &errs));
	CHECK(root.isMember("gensas"));
	CHECK_EQ(root["gensas"]["version"].asInt(), kResultFormatVersion);
	CHECK(root["gensas"]["format"].asString() == "result");
	CHECK_EQ(root["time"].size(), (Json::ArrayIndex)r.time.size());
	CHECK_EQ(root["x"].size(), (Json::ArrayIndex)r.time.size());
	CHECK_EQ(root["events"].size(), (Json::ArrayIndex)1);
	// Unwritable path reports an error instead of failing silently.
	CHECK(!makeResultWriter("json")->write(r, "/nonexistent-dir/x.json", err));
	CHECK(!err.empty());
}

static void testCsv() {
	Run run;
	if (!simulate(run, kModel, 2.0)) {
		return;
	}
	SasResult r = buildResult(*run.sol, 0.5);
	std::string err;
	CHECK(makeResultWriter("csv")->write(r, "result_test.csv", err));
	std::vector<std::string> lines = readLines("result_test.csv");
	CHECK_EQ(lines.size(), r.time.size() + 1);
	if (lines.empty()) {
		return;
	}
	CHECK(lines[0] == "time,x,a");
	// The second sample at the event time shows the new value of a.
	bool sawBefore = false;
	bool sawAfter = false;
	for (size_t j = 0; j < r.time.size() && j + 1 < lines.size(); j++) {
		std::istringstream is(lines[j + 1]);
		std::string cell;
		std::vector<double> cells;
		while (std::getline(is, cell, ',')) {
			cells.push_back(std::stod(cell));
		}
		CHECK_EQ(cells.size(), (size_t)3);
		if (cells.size() == 3) {
			CHECK_NEAR(cells[0], r.time[j], 1e-12);
			CHECK_NEAR(cells[1], r.values[0][j], 1e-12 * (1 + std::fabs(r.values[0][j])));
			if (std::fabs(cells[0] - 1.0) < 1e-9) {
				sawBefore = sawBefore || cells[2] == 1.0;
				sawAfter = sawAfter || cells[2] == 3.0;
			}
		}
	}
	CHECK(sawBefore && sawAfter);
	std::vector<std::string> ev = readLines("result_test.events.csv");
	CHECK_EQ(ev.size(), (size_t)2);
	if (ev.size() == 2) {
		CHECK(ev[0] == "time,kind,name,value");
		CHECK(ev[1].find(",assign,a,3") != std::string::npos);
	}
	CHECK(makeResultWriter("hdf5") == nullptr);
}

int main() {
	testBuildResult();
	testJsonVersion();
	testCsv();
	return sastest::finish("result_test");
}
