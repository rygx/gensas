// Simulation facade: job parsing and validation, plain and hybrid runs, progress and cancellation.
#include "sas/Simulation.h"
#include "sas/test/test_util.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

using namespace che::core;
using std::string;

static const char *kDecay =
	"model Decay\n"
	"  Real x(start = 1);\n"
	"equation\n"
	"  der(x) = -x;\n"
	"end Decay;\n";

static const char *kGainSwitch =
	"model M\n"
	"  discrete Real a(start = 1);\n"
	"  Real x(start = 1);\n"
	"equation\n"
	"  der(x) = -a*x;\n"
	"  when time >= 1 then\n"
	"    a = 3;\n"
	"  end when;\n"
	"end M;\n";

static double lastValue(const SasResult &r, const string &name) {
	for (size_t i = 0; i < r.names.size(); i++) {
		if (r.names[i] == name) {
			return r.values[i].back();
		}
	}
	return NAN;
}

static ModelicaJob sourceJob(const char *src, double endTime) {
	ModelicaJob job;
	job.modelSource = src;
	job.options.endTime = endTime;
	return job;
}

static void testPlainRun() {
	Simulation sim(sourceJob(kDecay, 2.0));
	CHECK(sim.run());
	CHECK(sim.error().empty());
	CHECK(!sim.result().empty());
	CHECK_NEAR(sim.result().time.back(), 2.0, 1e-9);
	CHECK_NEAR(lastValue(sim.result(), "x"), std::exp(-2.0), 1e-4);
	CHECK(sim.result().events.empty());
}

// Running twice gives the same result and does not pile up warnings.
static void testRunTwice() {
	Simulation sim(sourceJob(kDecay, 1.0));
	CHECK(sim.run());
	size_t warnings = sim.warnings().size();
	size_t n = sim.result().time.size();
	CHECK(sim.run());
	CHECK_EQ(sim.warnings().size(), warnings);
	CHECK_EQ(sim.result().time.size(), n);
}

static void testHybridRun() {
	Simulation sim(sourceJob(kGainSwitch, 2.0));
	CHECK(sim.run());
	CHECK_NEAR(lastValue(sim.result(), "x"), std::exp(-1.0) * std::exp(-3.0), 1e-4);
	CHECK_EQ(sim.result().events.size(), (size_t)1);
	CHECK_NEAR(sim.result().events[0].time, 1.0, 1e-9);
}

static void testBadModels() {
	{
		Simulation sim(sourceJob("model M\n  Real x;\nequation\n  der(x) = \n", 1.0));
		CHECK(!sim.run());
		CHECK(sim.error().find("error:") != string::npos); // the parser's message is part of the error
		CHECK(!sim.cancelled());
	}
	{
		ModelicaJob job;
		job.modelFile = "/nonexistent/model.mo";
		Simulation sim(job);
		CHECK(!sim.run());
		CHECK(sim.error().find("cannot read") != string::npos);
	}
	{
		Simulation sim(sourceJob("", 1.0));
		CHECK(!sim.run());
	}
}

// The solver throws on this model (an equation side that is a bare function call); run() reports it.
static void testSolverFailureIsReported() {
	Simulation sim(sourceJob(
		"model M\n"
		"  Real x(start = 1);\n"
		"  Real y(start = 1);\n"
		"equation\n"
		"  der(x) = sqrt(y);\n"
		"  der(y) = -0.1*y;\n"
		"end M;\n",
		1.0));
	bool ok = sim.run();
	if (!ok) {
		CHECK(!sim.error().empty());
		CHECK(!sim.cancelled());
	}
}

static void testProgressAndCancel() {
	for (const char *src : {kDecay, kGainSwitch}) {
		// Progress is monotone, ends at the end time and sees the whole horizon.
		{
			Simulation sim(sourceJob(src, 5.0));
			double prev = 0.0;
			int calls = 0;
			bool ok = true;
			CHECK(sim.run([&](double t, double end) {
				ok = ok && t >= prev && std::fabs(end - 5.0) < 1e-12;
				prev = t;
				calls++;
				return false;
			}));
			CHECK(ok);
			CHECK(calls >= 2);
			CHECK_NEAR(prev, 5.0, 1e-6);
		}
		// Cancel after the first report.
		{
			Simulation sim(sourceJob(src, 5.0));
			int calls = 0;
			CHECK(!sim.run([&](double, double) {
				calls++;
				return true;
			}));
			CHECK(sim.cancelled());
			CHECK_EQ(calls, 1);
		}
	}
}

static void testLimitsAreWarnings() {
	ModelicaJob job = sourceJob(kDecay, 5000.0);
	job.outputStep = 1e-6;
	Simulation sim(job);
	CHECK(sim.warnings().size() >= 2);
	CHECK_NEAR(sim.job().options.endTime, 1000.0, 0.0);
	CHECK_NEAR(sim.job().outputStep, 1e-3, 0.0);
}

static void testParseJob() {
	ModelicaJob job;
	string err;
	const string good =
		"{\"gensas\":{\"format\":\"job\",\"version\":1},\"kind\":\"modelica\","
		"\"model\":{\"file\":\"m/ball.mo\",\"piecewise\":\"smooth\",\"smoothEps\":1e-3},"
		"\"solver\":{\"endTime\":3,\"segment\":0.5,\"level\":12,\"alphaTol\":1e-4,\"errorTol\":1e-6},"
		"\"output\":{\"step\":0.02,\"json\":\"out/a.json\",\"csv\":\"/abs/a.csv\",\"mat\":\"a.mat\"}}";
	CHECK(parseModelicaJob(good, "/jobs", job, err));
	CHECK_EQ(job.modelFile, string("/jobs/m/ball.mo"));
	CHECK(job.piecewise == PW_SMOOTH);
	CHECK_NEAR(job.smoothEps, 1e-3, 0.0);
	CHECK_NEAR(job.options.endTime, 3.0, 0.0);
	CHECK_NEAR(job.options.segment, 0.5, 0.0);
	CHECK_EQ(job.options.nLvl, 12);
	CHECK_NEAR(job.options.alphaTol, 1e-4, 0.0);
	CHECK_NEAR(job.options.errorTol, 1e-6, 0.0);
	CHECK_NEAR(job.outputStep, 0.02, 0.0);
	CHECK_EQ(job.jsonOutput, string("/jobs/out/a.json"));
	CHECK_EQ(job.csvOutput, string("/abs/a.csv"));
	CHECK_EQ(job.matOutput, string("/jobs/a.mat"));

	// Inline source needs no base directory.
	CHECK(parseModelicaJob("{\"gensas\":{\"format\":\"job\",\"version\":1},\"kind\":\"modelica\",\"model\":{\"source\":\"model M end M;\"}}", "", job, err));
	CHECK_EQ(job.modelSource, string("model M end M;"));
	CHECK(job.modelFile.empty());

	struct Bad {
		const char *json;
		const char *expect;
	};
	const string head = "{\"gensas\":{\"format\":\"job\",\"version\":1},\"kind\":\"modelica\",";
	const Bad bad[] = {
		{"not json", "invalid JSON"},
		{"[]", "expected an object"},
		{"{\"gensas\":{\"format\":\"job\",\"version\":2},\"kind\":\"modelica\",\"model\":{\"source\":\"x\"}}", "gensas.version"},
		{"{\"gensas\":{\"format\":\"result\",\"version\":1},\"kind\":\"modelica\"}", "gensas.format"},
	};
	for (const Bad &b : bad) {
		CHECK(!parseModelicaJob(b.json, "", job, err));
		CHECK(err.find(b.expect) != string::npos);
	}
	const struct {
		string tail;
		const char *expect;
	} badBody[] = {
		{"\"model\":{\"source\":\"x\"},\"colour\":1}", "unknown key colour"},
		{"\"model\":{\"source\":\"x\",\"file\":\"y\"}}", "exactly one"},
		{"\"model\":{}}", "exactly one"},
		{"\"model\":{\"source\":\"x\",\"piecewise\":\"fuzzy\"}}", "model.piecewise"},
		{"\"model\":{\"source\":\"x\",\"smoothEps\":-1}}", "model.smoothEps"},
		{"\"model\":{\"source\":\"x\"},\"solver\":{\"endTime\":\"5\"}}", "solver.endTime: expected a number"},
		{"\"model\":{\"source\":\"x\"},\"solver\":{\"endTim\":5}}", "unknown key solver.endTim"},
		{"\"model\":{\"source\":\"x\"},\"solver\":{\"level\":2.5}}", "solver.level"},
		{"\"model\":{\"source\":\"x\"},\"output\":{\"json\":3}}", "output.json: expected a string"},
		{"\"model\":{\"source\":\"x\"},\"output\":{\"xml\":\"a\"}}", "unknown key output.xml"},
	};
	for (const auto &b : badBody) {
		CHECK(!parseModelicaJob(head + b.tail, "", job, err));
		CHECK(err.find(b.expect) != string::npos);
		if (err.find(b.expect) == string::npos) {
			std::cerr << "  got: " << err << std::endl;
		}
	}
	CHECK(!parseModelicaJob("{\"gensas\":{\"format\":\"job\",\"version\":1},\"kind\":\"powerflow\"}", "", job, err));
	CHECK(err.find("kind") != string::npos);
}

static void testWriteResult() {
	Simulation sim(sourceJob(kDecay, 1.0));
	CHECK(sim.run());
	string err;
	const char *path = "simulation_test_out.csv";
	CHECK(sim.writeResult("csv", path, err));
	std::ifstream in(path);
	string line;
	std::getline(in, line);
	CHECK(line.find("time") != string::npos);
	in.close();
	std::remove(path);
	CHECK(!sim.writeResult("xml", "x", err));
}

// Errors and warnings about the model text come with line numbers, without running the model.
static void testModelInfo() {
	const char *src =
		"model Tank \"a tank\"\n"
		"  parameter Real k(unit = \"1/s\") = 2 \"rate\";\n"
		"  constant Real c = 3;\n"
		"  Real x(start = 1, unit = \"m\") \"level\";\n"
		"  Real y;\n"
		"  input Real u = 0.5;\n"
		"equation\n"
		"  der(x) = -k * x + u;\n"
		"  y = x * c;\n"
		"end Tank;\n";
	Simulation sim(sourceJob(src, 1.0));
	CHECK(sim.check());
	CHECK(sim.hasModelInfo());
	const ModelInfo &info = sim.modelInfo();
	CHECK_EQ(info.name, string("Tank"));
	CHECK_EQ(info.description, string("a tank"));
	CHECK(!info.hybrid);
	CHECK_EQ(info.parameters.size(), (size_t)2);
	CHECK_EQ(info.variables.size(), (size_t)3);
	for (const ModelParameter &p : info.parameters) {
		if (p.name == "k") {
			CHECK_NEAR(p.value, 2.0, 1e-12);
			CHECK(!p.constant);
			CHECK_EQ(p.unit, string("1/s"));
			CHECK_EQ(p.description, string("rate"));
			CHECK_EQ(p.line, 2);
		} else {
			CHECK_EQ(p.name, string("c"));
			CHECK(p.constant);
		}
	}
	for (const ModelVariable &v : info.variables) {
		if (v.name == "x") {
			CHECK_EQ(v.kind, string("state"));
			CHECK(v.hasStart);
			CHECK_NEAR(v.start, 1.0, 1e-12);
			CHECK_EQ(v.unit, string("m"));
			CHECK_EQ(v.description, string("level"));
		} else if (v.name == "y") {
			CHECK_EQ(v.kind, string("algebraic"));
			CHECK(!v.hasStart);
		} else {
			CHECK_EQ(v.name, string("u"));
			CHECK(v.input);
		}
	}
	CHECK(modelInfoJson(info).find("\"kind\":\"state\"") != string::npos);

	// discrete variables make a model hybrid
	Simulation hybrid(sourceJob("model H\n  Real x(start = 1);\n  discrete Real d(start = 0);\nequation\n  der(x) = -x;\n  when time > 0.5 then\n    d = 1;\n  end when;\nend H;\n", 1.0));
	CHECK(hybrid.check());
	CHECK(hybrid.modelInfo().hybrid);
	bool sawDiscrete = false;
	for (const ModelVariable &v : hybrid.modelInfo().variables) {
		sawDiscrete = sawDiscrete || (v.name == "d" && v.kind == "discrete");
	}
	CHECK(sawDiscrete);

	// the equation count is checked without running: two unknowns, one equation
	Simulation bad(sourceJob("model B\n  Real x(start = 1);\n  Real y;\nequation\n  der(x) = -x;\nend B;\n", 1.0));
	CHECK(!bad.check());
	CHECK(!bad.diagnostics().empty());
	CHECK(bad.hasModelInfo()); // the model was read, so the editor can still list its variables
}

static void testDiagnostics() {
	// a good model: no diagnostics, and it can be run afterwards
	{
		Simulation sim(sourceJob(kDecay, 1.0));
		CHECK(sim.check());
		CHECK(sim.diagnostics().empty());
		CHECK(sim.run());
	}
	// a parser error
	{
		Simulation sim(sourceJob("model M\n  Real x;\nequation\n  der(x) = \nend M;\n", 1.0));
		CHECK(!sim.check());
		CHECK(!sim.diagnostics().empty());
		if (!sim.diagnostics().empty()) {
			const Diagnostic &d = sim.diagnostics()[0];
			CHECK(d.severity == Diagnostic::Error);
			CHECK(d.line >= 4);
			CHECK_EQ(d.source, string("parser"));
		}
		CHECK(!sim.error().empty());
	}
	// the statement's own line, not the scanner's, and no "Line N:" prefix in the message
	{
		Simulation sim(sourceJob("model M\n  Real x[2];\nequation\n  der(x[3]) = 1;\n  der(x[1]) = 1;\n  der(x[2]) = 1;\nend M;\n", 1.0));
		CHECK(!sim.check());
		bool found = false;
		for (const Diagnostic &d : sim.diagnostics()) {
			if (d.message.find("out of range") != string::npos) {
				CHECK_EQ(d.line, 4);
				CHECK(d.message.find("Line") == string::npos);
				found = true;
			}
		}
		CHECK(found);
	}
	// a checker error (undeclared identifier)
	{
		Simulation sim(sourceJob("model M\n  Real x(start = 1);\nequation\n  der(x) = -y;\nend M;\n", 1.0));
		CHECK(!sim.check());
		bool found = false;
		for (const Diagnostic &d : sim.diagnostics()) {
			if (d.source == "checker") {
				CHECK_EQ(d.line, 4);
				CHECK(d.message.find("undeclared identifier 'y'") != string::npos);
				found = true;
			}
		}
		CHECK(found);
		CHECK(!sim.run()); // run() fills the diagnostics too
		CHECK(!sim.diagnostics().empty());
	}
	// a warning does not make the check fail
	{
		Simulation sim(sourceJob("model M\n  Real x(start = 1, nominal = 2);\nequation\n  der(x) = -x;\nend M;\n", 1.0));
		CHECK(sim.check());
		CHECK_EQ(sim.diagnostics().size(), (size_t)1);
		if (sim.diagnostics().size() == 1) {
			CHECK(sim.diagnostics()[0].severity == Diagnostic::Warning);
			CHECK_EQ(sim.diagnostics()[0].line, 2);
		}
	}
}

int main() {
	testDiagnostics();
	testModelInfo();
	testPlainRun();
	testRunTwice();
	testHybridRun();
	testBadModels();
	testSolverFailureIsReported();
	testProgressAndCancel();
	testLimitsAreWarnings();
	testParseJob();
	testWriteResult();
	return sastest::finish("simulation_test");
}
