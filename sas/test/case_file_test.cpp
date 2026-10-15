// GenSAS case files (binary and JSON encodings) and the PSAT reader on top of them.
#include "io/CaseFile.h"
#include "io/MatPsatDataRW.h"
#include "sas/test/test_util.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace che::io;

static Matrix makeMatrix(size_t rows, size_t cols, double seed) {
	Matrix m;
	m.dims[0] = rows;
	m.dims[1] = cols;
	for (size_t c = 0; c < cols; c++) {
		for (size_t r = 0; r < rows; r++) {
			m.values.push_back(seed + 100.0 * (double)r + (double)c / 3.0);
		}
	}
	return m;
}

static MatrixSet sampleSet() {
	MatrixSet s;
	s.set("bus", makeMatrix(3, 6, 1.0));
	s.set("one", makeMatrix(1, 1, -2.5e-7));
	s.set("empty", makeMatrix(0, 0, 0.0));
	s.set("wide", makeMatrix(2, 40, 1e9));
	return s;
}

static void sameSet(const MatrixSet &a, const MatrixSet &b) {
	CHECK_EQ(a.all().size(), b.all().size());
	for (const auto &kv : a.all()) {
		const Matrix *o = b.find(kv.first);
		CHECK(o != NULL);
		if (o == NULL) {
			continue;
		}
		CHECK_EQ(o->dims[0], kv.second.dims[0]);
		CHECK_EQ(o->dims[1], kv.second.dims[1]);
		CHECK_EQ(o->values.size(), kv.second.values.size());
		for (size_t i = 0; i < o->values.size() && i < kv.second.values.size(); i++) {
			CHECK_NEAR(o->values[i], kv.second.values[i], 0.0); // exact round trip
		}
	}
}

static void testRoundTrip() {
	MatrixSet in = sampleSet();
	std::string err;
	MatrixSet bin, js;
	CHECK(writeCaseFile("case_test.case", in, err));
	CHECK(readCaseFile("case_test.case", bin, err));
	sameSet(in, bin);
	CHECK(writeCaseJson("case_test.case.json", in, err));
	CHECK(readCaseFile("case_test.case.json", js, err));
	sameSet(in, js);
}

static std::string slurp(const char *path) {
	std::ifstream f(path, std::ios::binary);
	return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

static void spit(const char *path, const std::string &s) {
	std::ofstream f(path, std::ios::binary);
	f << s;
}

static void testBadFiles() {
	std::string err;
	MatrixSet out;
	CHECK(!readCaseFile("does-not-exist.case", out, err));
	CHECK(!err.empty());

	std::string good = slurp("case_test.case");
	CHECK(good.size() > 100);
	spit("case_trunc.case", good.substr(0, good.size() - 16)); // data cut off
	CHECK(!readCaseFile("case_trunc.case", out, err));
	CHECK(err.find("truncated") != std::string::npos);
	spit("case_hdr.case", good.substr(0, 20)); // header cut off
	CHECK(!readCaseFile("case_hdr.case", out, err));
	std::string badMagic = good;
	badMagic[3] = 'X';
	spit("case_magic.case", badMagic);
	CHECK(!readCaseFile("case_magic.case", out, err));

	spit("case_v2.json", "{\"gensas\":{\"format\":\"case\",\"version\":2},\"matrices\":{}}");
	CHECK(!readCaseFile("case_v2.json", out, err));
	CHECK(err.find("version") != std::string::npos);
	spit("case_other.json", "{\"hello\":1}");
	CHECK(!readCaseFile("case_other.json", out, err));
	spit("case_ragged.json", "{\"gensas\":{\"format\":\"case\",\"version\":1},\"matrices\":{\"bus\":[[1,2],[3]]}}");
	CHECK(!readCaseFile("case_ragged.json", out, err));
	CHECK(err.find("different lengths") != std::string::npos);
	spit("case_text.json", "{\"gensas\":{\"format\":\"case\",\"version\":1},\"matrices\":{\"bus\":[[1,\"a\"]]}}");
	CHECK(!readCaseFile("case_text.json", out, err));
}

// The reader maps columns to fields; optional matrices may be missing, bus is mandatory.
static void testPsatReader() {
	MatrixSet s;
	Matrix bus = makeMatrix(2, 6, 0.0);
	// rows: [number, baseV, vMag, vAng, area, region]
	bus.values = {1, 2, 138, 230, 1.02, 0.99, 0.0, -0.1, 1, 1, 1, 1};
	s.set("bus", bus);
	Matrix pq = makeMatrix(1, 9, 0.0);
	pq.values = {2, 100, 230, 0.5, 0.2, 1.1, 0.9, 1, 1};
	s.set("pq", pq);
	std::string err;
	CHECK(writeCaseFile("case_psat.case", s, err));

	che::io::chedata::MatPsatReader reader;
	che::io::chedata::PsatDataSet data;
	CHECK_EQ(reader.parse("case_psat.case", &data), CHE_IO_SUCCESS);
	CHECK_EQ(data.nBus, 2);
	CHECK_EQ(data.nPq, 1);
	CHECK(data.nSw <= 0); // absent
	if (data.nBus == 2 && data.nPq == 1) {
		CHECK_EQ(data.buses[1].busNumber, 2);
		CHECK_NEAR(data.buses[1].baseV, 230.0, 0.0);
		CHECK_NEAR(data.buses[0].vMag, 1.02, 0.0);
		CHECK_NEAR(data.pqs[0].P, 0.5, 0.0);
		CHECK_NEAR(data.pqs[0].vMin, 0.9, 0.0);
		CHECK_EQ(data.pqs[0].busNumber, 2);
	}

	// Wrong column count and a missing bus matrix are rejected.
	MatrixSet bad;
	bad.set("bus", makeMatrix(2, 5, 0.0));
	che::io::chedata::PsatDataSet data2;
	CHECK(reader.parse(bad, &data2) != CHE_IO_SUCCESS);
	MatrixSet nobus;
	nobus.set("pq", pq);
	che::io::chedata::PsatDataSet data3;
	CHECK(reader.parse(nobus, &data3) != CHE_IO_SUCCESS);
}

int main() {
	testRoundTrip();
	testBadFiles();
	testPsatReader();
	return sastest::finish("case_file_test");
}
