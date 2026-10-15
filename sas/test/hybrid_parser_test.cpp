// Tests for the hybrid-Modelica subset accepted by the parser (Integer/Boolean
// variables, attributes, parameter expressions, if-expressions, if/when equations)
// and for the model checker that reports what the solver cannot run yet.
#include "sas/ModelChecker.h"
#include "sas/SasExpr.h"
#include "sas/SasLexico.h"
#include "sas/SasModel.h"
#include "sas/test/test_util.h"

#include <cstring>
#include <memory>
#include <sstream>
#include <string>

using namespace che::core;
using std::shared_ptr;
using std::string;

struct Parsed {
	SasModelParser parser;
	shared_ptr<SasModel> model;
};

static void parseInto(Parsed &p, const string &src) {
	p.parser.parseString(src.c_str());
	p.model = p.parser.pool.createNewTopModel();
	p.parser.currentModel = p.model;
	parseSingleSasModel(&p.parser);
}

static shared_ptr<IdPearl> find(const Parsed &p, const char *name) {
	return p.model->lookupIdFromModel(name, (int)strlen(name));
}

static SolverFeatures allFeatures() {
	SolverFeatures f;
	f.discreteVariables = f.ifExpressions = f.ifEquations = f.whenEquations = true;
	f.piecewiseFunctions = f.eventFunctions = f.simulationTime = f.terminateAndAssert = true;
	return f;
}

// Runs the checker; returns the number of errors and stores the text in `text`.
static int check(Parsed &p, const SolverFeatures &f, string &text) {
	std::ostringstream os;
	int n = checkModel(*p.model, f, "test.mo", os);
	text = os.str();
	return n;
}

static bool contains(const string &text, const char *what) {
	return text.find(what) != string::npos;
}

static const char *kBall =
	"model Ball\n"
	"  Real h(start=1.0) \"height\";\n"
	"  Real v(start=0.0) \"velocity\";\n"
	"  parameter Real g=9.81;\n"
	"  parameter Real e=0.8;\n"
	"equation\n"
	"  der(h)=v;\n"
	"  der(v)=-g;\n"
	"  when h<0 then\n"
	"    reinit(v, -e*pre(v));\n"
	"  end when;\n"
	"end Ball;\n";

static const char *kThermostat =
	"model Thermostat\n"
	"  Real T(start=20.0);\n"
	"  Boolean on(start=true);\n"
	"  parameter Real k=0.1;\n"
	"  parameter Real Ta=10.0;\n"
	"  parameter Real P=5.0;\n"
	"  parameter Real Tlow=18.0;\n"
	"  parameter Real Thigh=22.0;\n"
	"equation\n"
	"  der(T) = if on then -k*(T-Ta)+P else -k*(T-Ta);\n"
	"  when T < Tlow then\n"
	"    on = true;\n"
	"  elsewhen T > Thigh then\n"
	"    on = false;\n"
	"  end when;\n"
	"end Thermostat;\n";

static void testBouncingBallParses() {
	Parsed p;
	parseInto(p, kBall);
	CHECK_EQ(p.parser.errCount, 0);
	shared_ptr<IdPearl> h = find(p, "h");
	CHECK(h != nullptr);
	if (h) {
		CHECK_EQ((int)h->idTy, (int)ID_VAR);
		CHECK(h->initVal);
		CHECK_NEAR(h->value.d, 1.0, 0.0);
	}
	CHECK_EQ(p.model->modelTrees.size(), (size_t)3);
	// Third statement is the when-equation, with the reinit() call as its body.
	AstNode *when = p.model->modelTrees.back()->pHead;
	CHECK(when != nullptr);
	if (when) {
		CHECK_EQ((int)when->nodeKind, (int)NK_WhenEquation);
		CHECK_EQ(when->line, 9);
		CHECK(when->subs[0] != nullptr && when->subs[0]->op == OP_LESS);
		AstNode *body = when->subs[1] ? when->subs[1]->subs[0] : nullptr;
		CHECK(body != nullptr);
		if (body) {
			CHECK_EQ((int)body->op, (int)OP_CALL);
			CHECK(body->next == nullptr);
		}
		CHECK(when->subs[1] != nullptr && when->subs[1]->subs[1] == nullptr);
	}
}

static void testBouncingBallChecker() {
	Parsed p;
	parseInto(p, kBall);
	string text;
	// The current solver cannot run when-equations, reinit() or pre() yet.
	int n = check(p, SolverFeatures(), text);
	CHECK_EQ(n, 3);
	CHECK(contains(text, "test.mo:9: error: when-equations are not supported"));
	CHECK(contains(text, "function reinit() is not supported"));
	CHECK(contains(text, "function pre() is not supported"));
	CHECK_EQ(check(p, allFeatures(), text), 0);
}

static void testThermostat() {
	Parsed p;
	parseInto(p, kThermostat);
	CHECK_EQ(p.parser.errCount, 0);
	shared_ptr<IdPearl> on = find(p, "on");
	CHECK(on != nullptr);
	if (on) {
		CHECK_EQ((int)on->idTy, (int)ID_DISC);
		CHECK_EQ((int)on->ty, (int)TY_BOOL);
		CHECK(on->initVal);
		CHECK_EQ(on->value.i, 1);
	}
	CHECK_EQ(p.model->modelTrees.size(), (size_t)2);
	// der(T) = if ... : right-hand side is an if-expression.
	AstNode *eq = p.model->modelTrees.front()->pHead;
	CHECK(eq != nullptr && eq->op == OP_EQN);
	if (eq) {
		CHECK(eq->subs[1] != nullptr && eq->subs[1]->op == OP_IFEXPR);
	}
	// elsewhen is chained under the first when.
	AstNode *when = p.model->modelTrees.back()->pHead;
	CHECK(when != nullptr && when->nodeKind == NK_WhenEquation);
	if (when) {
		CHECK(when->subs[1] != nullptr && when->subs[1]->subs[1] != nullptr);
		CHECK_EQ((int)when->subs[1]->subs[1]->nodeKind, (int)NK_WhenEquation);
		CHECK_EQ(when->subs[1]->subs[1]->line, 13);
	}
	string text;
	CHECK(check(p, SolverFeatures(), text) >= 3);
	CHECK(contains(text, "discrete variable 'on'"));
	CHECK(contains(text, "if-expressions are not supported"));
	CHECK_EQ(check(p, allFeatures(), text), 0);
}

static void testElseIfExpression() {
	Parsed p;
	parseInto(p, "model M\n  Real x(start=0);\n  parameter Real a=1;\nequation\n"
				 "  der(x) = if x < 1 then a elseif x < 2 then 2*a else 3*a;\nend M;\n");
	CHECK_EQ(p.parser.errCount, 0);
	AstNode *rhs = p.model->modelTrees.front()->pHead->subs[1];
	CHECK(rhs != nullptr && rhs->op == OP_IFEXPR);
	if (rhs && rhs->subs[1]) {
		CHECK_EQ((int)rhs->subs[1]->op, (int)OP_BRANCH);
		// else part is a nested if-expression (the elseif)
		CHECK(rhs->subs[1]->subs[1] != nullptr && rhs->subs[1]->subs[1]->op == OP_IFEXPR);
	}
}

static void testIfEquation() {
	Parsed p;
	parseInto(p, "model M\n  Real x(start=0);\n  Real y;\nequation\n"
				 "  der(x) = 1;\n"
				 "  if x < 1 then\n    y = 1;\n  elseif x < 2 then\n    y = 2;\n  else\n    y = 3;\n  end if;\n"
				 "end M;\n");
	CHECK_EQ(p.parser.errCount, 0);
	CHECK_EQ(p.model->modelTrees.size(), (size_t)2);
	AstNode *n = p.model->modelTrees.back()->pHead;
	CHECK(n != nullptr && n->nodeKind == NK_IfEquation);
	if (n) {
		CHECK_EQ(n->line, 6);
		AstNode *thenList = n->subs[1]->subs[0];
		CHECK(thenList != nullptr && thenList->op == OP_EQN);
		// elseif is desugared into a nested if as the sole statement of the else list
		AstNode *elseList = n->subs[1]->subs[1];
		CHECK(elseList != nullptr && elseList->nodeKind == NK_IfEquation);
		if (elseList) {
			CHECK_EQ(elseList->line, 8);
			CHECK(elseList->subs[1]->subs[1] != nullptr && elseList->subs[1]->subs[1]->op == OP_EQN);
		}
	}
	string text;
	CHECK_EQ(check(p, SolverFeatures(), text), 2); // if at line 6 and the nested elseif at line 8
	CHECK_EQ(check(p, allFeatures(), text), 0);
}

static void testAttributesAndExpressions() {
	Parsed p;
	parseInto(p, "model M\n"
				 "  parameter Real a = 2.0;\n"
				 "  parameter Real b = a*3 + sqrt(4.0);\n"
				 "  parameter Integer n = 3;\n"
				 "  Real x(start=b/2, unit=\"m\", fixed=true) \"position\";\n"
				 "initial equation\n"
				 "  der(x) = -a;\n"
				 "equation\n"
				 "  der(x) = -a*x;\n"
				 "end M;\n");
	CHECK_EQ(p.parser.errCount, 0);
	shared_ptr<IdPearl> b = find(p, "b");
	CHECK(b != nullptr);
	if (b) {
		CHECK_NEAR(b->value.d, 8.0, 1e-12);
	}
	shared_ptr<IdPearl> n = find(p, "n");
	CHECK(n != nullptr);
	if (n) {
		CHECK_EQ((int)n->idTy, (int)ID_PARAM);
		CHECK_EQ(n->value.i, 3);
	}
	shared_ptr<IdPearl> x = find(p, "x");
	CHECK(x != nullptr);
	if (x) {
		CHECK_NEAR(x->value.d, 4.0, 1e-12);
		CHECK(x->initVal);
		CHECK(x->initDer);
		CHECK_NEAR(x->derVal.d, -2.0, 1e-12);
		CHECK(x->unit != nullptr && std::string((const char *)x->unit) == "m");
	}
	// Attribute names are not model symbols.
	CHECK(find(p, "start") == nullptr);
	CHECK(find(p, "fixed") == nullptr);
	string text;
	CHECK_EQ(check(p, SolverFeatures(), text), 0);
}

static void testSectionOrder() {
	// No initial equation section at all; several equation sections.
	Parsed p;
	parseInto(p, "model M\n  Real x(start=1);\n  Real y(start=2);\nequation\n  der(x) = -x;\nequation\n  der(y) = -y;\nend M;\n");
	CHECK_EQ(p.parser.errCount, 0);
	CHECK_EQ(p.model->modelTrees.size(), (size_t)2);
}

static void testInitialEquationExpression() {
	Parsed p;
	parseInto(p, "model M\n  parameter Real k=2;\n  Real x;\ninitial equation\n  x = 2*k + 1;\nequation\n  der(x) = -x;\nend M;\n");
	CHECK_EQ(p.parser.errCount, 0);
	shared_ptr<IdPearl> x = find(p, "x");
	CHECK(x != nullptr);
	if (x) {
		CHECK(x->initVal);
		CHECK_NEAR(x->value.d, 5.0, 1e-12);
	}
}

static int parseErrors(const string &src) {
	Parsed p;
	parseInto(p, src);
	return p.parser.errCount;
}

static void testSyntaxErrors() {
	CHECK(parseErrors("model M\n  Real x;\n  Real x;\nequation\n  der(x)=1;\nend M;\n") > 0);													  // redeclaration
	CHECK(parseErrors("model M\n  parameter Real k;\nequation\nend M;\n") > 0);																	  // parameter without value
	CHECK(parseErrors("model M\n  parameter Real k=1;\n  parameter Real j=k*time;\nequation\nend M;\n") > 0);									  // non-constant parameter
	CHECK(parseErrors("model M\n  discrete Real d=1;\nequation\nend M;\n") > 0);																  // binding on a discrete variable
	CHECK(parseErrors("model M\n  Real x=1;\nequation\n  der(x)=1;\nend M;\n") == 0);															  // binding on a continuous variable is an equation
	CHECK(parseErrors("model M\n  Real x;\nequation\n  reinit(x, 1);\nend M;\n") > 0);															  // reinit outside when
	CHECK(parseErrors("model M\n  Real x;\nequation\n  when x>0 then\n    when x>1 then\n      x=1;\n    end when;\n  end when;\nend M;\n") > 0); // nested when
	CHECK(parseErrors("model M\n  Real x;\nequation\n  if x>0 then\n    x=1;\n  else\n    x=2;\nend M;\n") > 0);								  // missing end if
	CHECK(parseErrors("model M\n  Real x;\nequation\n  der(x)=;\nend M;\n") > 0);																  // missing expression
	CHECK(parseErrors("model M\n  Real x;\nequation\n  x + 1;\nend M;\n") > 0);																	  // not an equation
	CHECK(parseErrors("model M\n  Real x;\ninitial equation\n  x = x + 1;\nequation\n  der(x)=1;\nend M;\n") > 0);								  // non-constant initial value
	CHECK(parseErrors("model M\n  Real x;\nequation\n  der(x)=1;\n  Real y;\nend M;\n") > 0);													  // declaration after equations
}

// Comments, annotations, several declarators in one declaration, input/output prefixes.
static void testSourceSyntax() {
	CHECK_EQ(parseErrors("model M \"text\"\n  /* block\n     comment */\n  Real x(start=1); // line\n  /**/\nequation\n  der(x) = -x; /* inline */\nend M;\n"), 0);
	CHECK_EQ(parseErrors("model M\n  Real x(start=1) \"a\", y(start=2), z;\nequation\n  der(x) = -x;\n  der(y) = -y;\n  z = x + y;\nend M;\n"), 0);
	CHECK_EQ(parseErrors("model M\n  parameter Real a=1, b=2;\n  Real x(start=1);\nequation\n  der(x) = -a*b*x;\nend M;\n"), 0);
	CHECK_EQ(parseErrors("model M\n  output Real y;\n  input Real u = 1;\n  Real x(start=1);\nequation\n  der(x) = -x + u;\n  y = 2*x;\nend M;\n"), 0);
	CHECK_EQ(parseErrors("model M\n  Real x(start=1) annotation(Dialog(group=\"g\"));\n  annotation(Icon(graphics={Rectangle(extent={{-1,-1},{1,1}})}));\nequation\n  der(x) = -x;\n  annotation(experiment(StopTime=1));\nend M;\n"), 0);
	CHECK_EQ(parseErrors("model M\n  Real x(start=1) annotation(Dialog(group=\"g\")), y(start=2) \"d\" annotation(a=1);\nequation\n  der(x) = -x;\n  der(y) = -y;\nend M;\n"), 0);
	CHECK(parseErrors("model M\n  discrete Real a(start=0);\nequation\n  when initial then\n    a = 1;\n  end when;\nend M;\n") > 0); // initial without ()
	CHECK(parseErrors("model M\n  input Real u;\nequation\nend M;\n") > 0);															  // input without a value
	CHECK(parseErrors("model M\n  Real x;\n  /* open\nequation\n  der(x) = -x;\nend M;\n") > 0);									  // unterminated comment
	CHECK(parseErrors("model M\n  extends Base;\n  Real x;\nequation\n  der(x) = -x;\nend M;\n") > 0);
	CHECK(parseErrors("model M\n  Real x;\nalgorithm\n  x := 1;\nend M;\n") > 0);
}

// Arrays and for-equations: what is accepted and what is rejected with an error.
static void testArraySyntax() {
	const string head = "model M\n  parameter Integer n = 3;\n  Real x[n](each start = 1);\n  Real y[2,2](each start = 0);\nequation\n";
	CHECK_EQ(parseErrors(head + "  for i in 1:n loop der(x[i]) = -x[i]; end for;\n  for i in 1:2, j in 1:2 loop der(y[i,j]) = -y[i,j]; end for;\nend M;\n"), 0);
	CHECK_EQ(parseErrors(head + "  der(x[1]) = -x[2]; der(x[2]) = 0; der(x[3]) = 0;\n  for i in 1:2 loop der(y[i,1]) = x[i+1]; der(y[i,2]) = 0; end for;\nend M;\n"), 0);
	CHECK(parseErrors(head + "  der(x) = -x;\nend M;\n") > 0);																		   // array without subscripts
	CHECK(parseErrors(head + "  der(x[4]) = -1;\nend M;\n") > 0);																	   // out of range
	CHECK(parseErrors(head + "  der(y[1]) = -1;\nend M;\n") > 0);																	   // too few subscripts
	CHECK(parseErrors(head + "  der(n[1]) = -1;\nend M;\n") > 0);																	   // not an array
	CHECK(parseErrors(head + "  for i in 1:n loop der(x[i+3]) = 1; end for;\nend M;\n") > 0);										   // out of range in a loop
	CHECK(parseErrors(head + "  for i in 1:q loop der(x[1]) = 1; end for;\nend M;\n") > 0);											   // range not constant
	CHECK(parseErrors(head + "  for i in 1:0 loop der(x[1]) = 1; end for;\n  der(x[1])=1; der(x[2])=1; der(x[3])=1;\nend M;\n") == 0); // empty loop
	CHECK(parseErrors("model M\n  Real x[0];\nequation\nend M;\n") > 0);															   // empty array
	CHECK(parseErrors("model M\n  Real x[2](start = {1, 2, 3});\nequation\n  der(x[1]) = 1;\nend M;\n") > 0);						   // wrong number of values
	CHECK(parseErrors("model M\n  Real x[2];\n  Real x[3];\nequation\nend M;\n") > 0);												   // redeclaration
	CHECK(parseErrors("model M\n  Real x[2];\nequation\n  for i in 1:2 loop\n    der(x[i]) = 1;\nend M;\n") > 0);					   // missing end for
}

static void testCheckerSemantics() {
	{
		Parsed p;
		parseInto(p, "model M\n  Real x(start=1);\nequation\n  der(x) = -y;\nend M;\n");
		CHECK_EQ(p.parser.errCount, 0);
		string text;
		CHECK_EQ(check(p, SolverFeatures(), text), 1);
		CHECK(contains(text, "test.mo:4: error: undeclared identifier 'y'"));
	}
	{
		Parsed p;
		parseInto(p, "model M\n  Real x(start=1);\nequation\n  der(x) = foo(x);\nend M;\n");
		string text;
		CHECK_EQ(check(p, SolverFeatures(), text), 1);
		CHECK(contains(text, "unknown function 'foo'"));
	}
	{
		Parsed p;
		parseInto(p, "model M\n  Real x(start=1);\nequation\n  der(x) = -time*x;\nend M;\n");
		string text;
		CHECK_EQ(check(p, SolverFeatures(), text), 1);
		CHECK(contains(text, "'time' is not supported"));
		SolverFeatures f;
		f.simulationTime = true;
		CHECK_EQ(check(p, f, text), 0);
	}
}

int main() {
	testBouncingBallParses();
	testBouncingBallChecker();
	testThermostat();
	testElseIfExpression();
	testIfEquation();
	testAttributesAndExpressions();
	testSectionOrder();
	testInitialEquationExpression();
	testSyntaxErrors();
	testSourceSyntax();
	testArraySyntax();
	testCheckerSemantics();
	return sastest::finish("hybrid_parser_test");
}
