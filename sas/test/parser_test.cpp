// Regression tests for the Modelica model parser (sas/ModelParser.cpp).
// They pin the behaviour of the parser for the model class supported today,
// so that parser changes can be made safely.
#include "sas/SasExpr.h"
#include "sas/SasLexico.h"
#include "sas/SasModel.h"
#include "sas/test/test_util.h"

#include <cstring>
#include <memory>
#include <string>

using namespace che::core;
using std::shared_ptr;

static const char *kOdeModel =
	"model exODE01\n"
	"  Real x1 \"x1\";\n"
	"  Real x2 \"x2\";\n"
	"  parameter Real k1=-1.0 \"k1\";\n"
	"  parameter Real k2=-2.0 \"k2\";\n"
	"initial equation\n"
	"  x1=3.0;\n"
	"  x2=2.0;\n"
	"equation\n"
	"  der(x1)=k1*x1+0.5*x2;\n"
	"  der(x2)=k2*x2;\n"
	"end exODE01;\n";

struct Parsed {
	SasModelParser parser;
	shared_ptr<SasModel> model;
};

static void parseInto(Parsed &p, const std::string &src) {
	p.parser.parseString(src.c_str());
	p.model = GlobalPool::getInstance().createNewTopModel();
	p.parser.currentModel = p.model;
	parseSingleSasModel(&p.parser);
}

static shared_ptr<IdPearl> find(const shared_ptr<SasModel> &m, const char *name) {
	return m->lookupIdFromModel(name, (int)strlen(name));
}

static void testValidModel() {
	Parsed p;
	parseInto(p, kOdeModel);
	CHECK_EQ(p.parser.errCount, 0);
	CHECK_EQ(p.model->errCount, 0);
	CHECK(p.model->modelInfo != nullptr);
	CHECK(std::string((const char *)p.model->modelInfo->str) == "exODE01");

	shared_ptr<IdPearl> x1 = find(p.model, "x1");
	CHECK(x1 != nullptr);
	if (x1) {
		CHECK_EQ(x1->idTy, ID_VAR);
		CHECK(x1->initVal);
		CHECK_NEAR(x1->value.d, 3.0, 0.0);
	}
	shared_ptr<IdPearl> x2 = find(p.model, "x2");
	CHECK(x2 != nullptr);
	if (x2) {
		CHECK_NEAR(x2->value.d, 2.0, 0.0);
	}
	shared_ptr<IdPearl> k1 = find(p.model, "k1");
	CHECK(k1 != nullptr);
	if (k1) {
		CHECK_EQ(k1->idTy, ID_PARAM);
		CHECK_NEAR(k1->value.d, -1.0, 0.0);
	}
	shared_ptr<IdPearl> k2 = find(p.model, "k2");
	CHECK(k2 != nullptr);
	if (k2) {
		CHECK_NEAR(k2->value.d, -2.0, 0.0);
	}

	CHECK_EQ(p.model->modelTrees.size(), (size_t)2);
	for (const auto &tree : p.model->modelTrees) {
		CHECK(tree->pHead != nullptr);
		if (tree->pHead) {
			CHECK_EQ((int)tree->pHead->op, (int)OP_EQN);
			// left-hand side is der(x)
			CHECK(tree->pHead->subs[0] != nullptr);
			if (tree->pHead->subs[0]) {
				CHECK_EQ((int)tree->pHead->subs[0]->op, (int)OP_CALL);
			}
		}
	}
}

static void testMissingEnd() {
	std::string src(kOdeModel);
	src = src.substr(0, src.find("end exODE01;"));
	Parsed p;
	parseInto(p, src);
	CHECK(p.parser.errCount > 0);
}

static void testBadEquation() {
	std::string src(kOdeModel);
	src.replace(src.find("der(x2)=k2*x2;"), 14, "der(x2)=;");
	Parsed p;
	parseInto(p, src);
	CHECK(p.parser.errCount > 0);
}

static void testMissingParameterValue() {
	std::string src(kOdeModel);
	src.replace(src.find("k1=-1.0"), 7, "k1");
	Parsed p;
	parseInto(p, src);
	CHECK(p.parser.errCount > 0);
}

int main() {
	testValidModel();
	testMissingEnd();
	testBadEquation();
	testMissingParameterValue();
	return sastest::finish("parser_test");
}
