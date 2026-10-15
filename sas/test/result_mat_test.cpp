// MAT v5 writer, read back with a separate minimal reader written from the file format description.
#include "sas/HybridSimulator.h"
#include "sas/ModelChecker.h"
#include "sas/SasExpr.h"
#include "sas/SasLexico.h"
#include "sas/SasResult.h"
#include "io/Mat5Writer.h"
#include "sas/SasResultMat.h"
#include "sas/test/test_util.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

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

// ---- Minimal level 5 reader (little-endian, uncompressed), independent of the writer ----

struct MatVar {
	unsigned arrayClass = 0;
	std::vector<int> dims;
	std::vector<unsigned char> data;
};

static bool readMat5(const char *path, std::map<std::string, MatVar> &out, std::string &err) {
	std::ifstream f(path, std::ios::binary);
	std::vector<unsigned char> b((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
	if (b.size() < 128 || std::string((const char *)b.data(), 19) != "MATLAB 5.0 MAT-file") {
		err = "bad header";
		return false;
	}
	if (b[124] != 0x00 || b[125] != 0x01 || b[126] != 'I' || b[127] != 'M') {
		err = "bad version or byte order mark";
		return false;
	}
	auto u32 = [&](size_t pos) {
		return (unsigned)b[pos] | ((unsigned)b[pos + 1] << 8) | ((unsigned)b[pos + 2] << 16) | ((unsigned)b[pos + 3] << 24);
	};
	size_t pos = 128;
	while (pos < b.size()) {
		if (pos + 8 > b.size() || u32(pos) != 14) {
			err = "expected a matrix element";
			return false;
		}
		size_t end = pos + 8 + u32(pos + 4);
		if (end > b.size() || end % 8 != 0) {
			err = "bad matrix element size";
			return false;
		}
		size_t p = pos + 8;
		MatVar v;
		std::string name;
		for (int part = 0; part < 4; part++) {
			if (p + 8 > end) {
				err = "truncated sub-element";
				return false;
			}
			unsigned type = u32(p);
			size_t n = u32(p + 4);
			size_t padded = (n + 7) & ~(size_t)7;
			if (p + 8 + padded > end) {
				err = "sub-element exceeds its matrix";
				return false;
			}
			const unsigned char *d = b.data() + p + 8;
			if (part == 0) { // array flags
				if (type != 6 || n != 8) {
					err = "bad array flags";
					return false;
				}
				v.arrayClass = u32(p + 8) & 0xff;
			} else if (part == 1) { // dimensions
				if (type != 5) {
					err = "bad dimensions type";
					return false;
				}
				for (size_t i = 0; i < n / 4; i++) {
					v.dims.push_back((int)u32(p + 8 + 4 * i));
				}
			} else if (part == 2) { // name
				if (type != 1) {
					err = "bad name type";
					return false;
				}
				name.assign((const char *)d, n);
			} else { // real part
				v.data.assign(d, d + n);
			}
			p += 8 + padded;
		}
		if (p != end) {
			err = "unexpected data after the real part";
			return false;
		}
		out[name] = v;
		pos = end;
	}
	return true;
}

static std::vector<double> asDoubles(const MatVar &v) {
	std::vector<double> d(v.data.size() / 8);
	if (!d.empty()) {
		std::memcpy(d.data(), v.data.data(), d.size() * 8);
	}
	return d;
}

// Row r of a char matrix, trailing blanks removed.
static std::string charRow(const MatVar &v, size_t r) {
	size_t rows = (size_t)v.dims[0];
	size_t cols = (size_t)v.dims[1];
	std::string s;
	for (size_t c = 0; c < cols && r < rows; c++) {
		s += (char)v.data[c * rows + r];
	}
	while (!s.empty() && s.back() == ' ') {
		s.pop_back();
	}
	return s;
}

static const MatVar *get(const std::map<std::string, MatVar> &m, const char *name) {
	std::map<std::string, MatVar>::const_iterator it = m.find(name);
	CHECK(it != m.end());
	return it == m.end() ? NULL : &it->second;
}

static void testMatFile() {
	Run run;
	if (!simulate(run, kModel, 2.0)) {
		return;
	}
	SasResult r = buildResult(*run.sol, 0.5);
	std::string err;
	CHECK(makeMatResultWriter()->write(r, "result_test.mat", err));
	std::map<std::string, MatVar> m;
	CHECK(readMat5("result_test.mat", m, err));
	if (!err.empty()) {
		std::cerr << err << std::endl;
		return;
	}
	const MatVar *ver = get(m, "gensas_version");
	const MatVar *t = get(m, "t");
	const MatVar *s = get(m, "s");
	const MatVar *d = get(m, "d");
	const MatVar *et = get(m, "event_t");
	const MatVar *names = get(m, "names");
	const MatVar *dnames = get(m, "dnames");
	if (!ver || !t || !s || !d || !et || !names || !dnames) {
		return;
	}
	CHECK_EQ(ver->arrayClass, 6u);
	CHECK_NEAR(asDoubles(*ver).at(0), kResultFormatVersion, 0.0);
	CHECK_EQ(t->dims[0], 1);
	CHECK_EQ((size_t)t->dims[1], r.time.size());
	CHECK_EQ(s->dims[0], 1);
	CHECK_EQ((size_t)s->dims[1], r.time.size());
	CHECK_EQ(d->dims[0], 1);
	std::vector<double> tv = asDoubles(*t), sv = asDoubles(*s), dv = asDoubles(*d), ev = asDoubles(*et);
	for (size_t j = 0; j < r.time.size() && j < tv.size() && j < sv.size() && j < dv.size(); j++) {
		CHECK_NEAR(tv[j], r.time[j], 0.0);
		CHECK_NEAR(sv[j], r.values[0][j], 0.0);
		CHECK_NEAR(dv[j], r.discrete[0][j], 0.0);
	}
	CHECK_EQ(ev.size(), (size_t)1);
	if (ev.size() == 1) {
		CHECK_NEAR(ev[0], 1.0, 1e-12);
	}
	CHECK_EQ(names->arrayClass, 4u);
	CHECK(charRow(*names, 0) == "x");
	CHECK(charRow(*dnames, 0) == "a");

	// Empty results are an error, not an empty file.
	SasResult empty;
	CHECK(!makeMatResultWriter()->write(empty, "result_test_empty.mat", err));
	CHECK(!err.empty());
}

// The writer itself: odd name lengths (padding), a larger matrix, char rows of different length.
static void testWriterRoundTrip() {
	che::io::Mat5Writer w;
	CHECK(w.open("mat5_roundtrip.mat"));
	std::vector<double> big(37 * 53);
	for (size_t i = 0; i < big.size(); i++) {
		big[i] = std::sin(0.001 * (double)i) * 1e5 - (double)i;
	}
	double scalar = -2.5;
	CHECK(w.writeMatrix("a", 1, 1, &scalar));
	CHECK(w.writeMatrix("odd_name_x", 37, 53, big.data()));
	CHECK(w.writeMatrix("seven77", 0, 0, NULL));
	CHECK(w.writeChars("txt", std::vector<std::string>{"x", "longer.name", "mid"}));
	CHECK(w.close());
	CHECK(w.error().empty());

	std::map<std::string, MatVar> m;
	std::string err;
	CHECK(readMat5("mat5_roundtrip.mat", m, err));
	if (!err.empty()) {
		std::cerr << err << std::endl;
		return;
	}
	const MatVar *a = get(m, "a");
	const MatVar *b = get(m, "odd_name_x");
	const MatVar *e = get(m, "seven77");
	const MatVar *txt = get(m, "txt");
	if (!a || !b || !e || !txt) {
		return;
	}
	CHECK_NEAR(asDoubles(*a).at(0), -2.5, 0.0);
	CHECK_EQ(b->dims[0], 37);
	CHECK_EQ(b->dims[1], 53);
	std::vector<double> bv = asDoubles(*b);
	CHECK_EQ(bv.size(), big.size());
	for (size_t i = 0; i < bv.size() && i < big.size(); i++) {
		CHECK_NEAR(bv[i], big[i], 0.0);
	}
	CHECK_EQ(e->data.size(), (size_t)0);
	CHECK(charRow(*txt, 0) == "x");
	CHECK(charRow(*txt, 1) == "longer.name");
	CHECK(charRow(*txt, 2) == "mid");

	// Errors are reported, not ignored.
	che::io::Mat5Writer bad;
	CHECK(!bad.open("/nonexistent-dir/x.mat"));
	CHECK(!bad.error().empty());
	che::io::Mat5Writer longName;
	CHECK(longName.open("mat5_longname.mat"));
	CHECK(!longName.writeMatrix(std::string(64, 'n'), 1, 1, &scalar));
	CHECK(!longName.close());
}

int main() {
	testMatFile();
	testWriterRoundTrip();
	return sastest::finish("result_mat_test");
}
