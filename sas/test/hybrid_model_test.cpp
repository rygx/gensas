// Tests for the hybrid model layer: expression evaluation, specialization for a discrete state,
// time-event scheduling and the stricter checker rules of the time-event simulator.
#include "sas/HybridModel.h"
#include "sas/ModelChecker.h"
#include "sas/SasExpr.h"
#include "sas/SasLexico.h"
#include "sas/SasModel.h"
#include "sas/test/test_util.h"

#include <cmath>
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
	p.model = GlobalPool::getInstance().createNewTopModel();
	p.parser.currentModel = p.model;
	parseSingleSasModel(&p.parser);
}

static SolverFeatures hybridFeatures() {
	SolverFeatures f;
	f.discreteVariables = f.ifExpressions = f.ifEquations = f.whenEquations = true;
	f.piecewiseFunctions = f.eventFunctions = f.simulationTime = true;
	f.timeEventsOnly = true;
	return f;
}

static int check(Parsed &p, string &text) {
	std::ostringstream os;
	int n = checkModel(*p.model, hybridFeatures(), "test.mo", os);
	text = os.str();
	return n;
}

static bool contains(const string &text, const char *what) {
	return text.find(what) != string::npos;
}

// Counts nodes matching: identifier `name` (if non-null) or op `op` (if name is null).
static int countNodes(const AstNode *p, const char *name, Op op) {
	int n = 0;
	for (; p != NULL; p = p->next) {
		if (name != NULL) {
			if (p->op == OP_ID && p->ty != TY_FUNCTION && p->value.p != NULL && std::strcmp((const char *)p->value.p, name) == 0) {
				n++;
			}
		} else if (p->op == op) {
			n++;
		}
		n += countNodes(p->subs[0], name, op) + countNodes(p->subs[1], name, op);
	}
	return n;
}

static const char *kSwitching =
	"model Switching\n"
	"  parameter Real a = 2;\n"
	"  discrete Real d(start = 3);\n"
	"  Boolean b(start = true);\n"
	"  Integer n(start = 0);\n"
	"  Real x(start = 1);\n"
	"equation\n"
	"  der(x) = if b then -d*x else a*x;\n"
	"  when sample(0.5, 1) then\n"
	"    n = pre(n) + 1;\n"
	"  end when;\n"
	"  when time >= 1.5 then\n"
	"    d = pre(d) + 1;\n"
	"    b = false;\n"
	"  end when;\n"
	"end Switching;\n";

static void testEnvAndEval() {
	Parsed p;
	parseInto(p, kSwitching);
	CHECK_EQ(p.parser.errCount, 0);
	HybridEnv env = makeHybridEnv(*p.model);
	CHECK_NEAR(env.values["a"], 2.0, 0.0);
	CHECK_NEAR(env.values["d"], 3.0, 0.0);
	CHECK_NEAR(env.values["b"], 1.0, 0.0);
	CHECK_NEAR(env.values["n"], 0.0, 0.0);
	CHECK(env.discrete.count("d") == 1 && env.discrete.count("b") == 1 && env.discrete.count("a") == 0);
	CHECK(isHybridModel(*p.model));

	// The right-hand side of the first equation is an if-expression: b -> -d*x branch needs x, so
	// evaluate the pieces that only use discrete data.
	AstNode *rhs = p.model->modelTrees.front()->pHead->subs[1];
	CHECK(rhs != NULL && rhs->op == OP_IFEXPR);
	if (rhs) {
		double v = 0;
		string err;
		CHECK(env.eval(rhs->subs[0], v, err));
		CHECK_NEAR(v, 1.0, 0.0);
		env.values["b"] = 0.0;
		CHECK(env.eval(rhs->subs[0], v, err));
		CHECK_NEAR(v, 0.0, 0.0);
		// x is not known to the environment.
		CHECK(!env.eval(rhs, v, err));
		CHECK(contains(err, "no value for"));
	}
}

static void testSpecialize() {
	{
		Parsed p;
		parseInto(p, kSwitching);
		HybridEnv env = makeHybridEnv(*p.model);
		string err;
		CHECK(specializeModel(*p.model, env, err));
		CHECK_EQ(p.model->modelTrees.size(), (size_t)1); // when-equations removed
		const AstNode *eq = p.model->modelTrees.front()->pHead;
		CHECK_EQ(countNodes(eq, NULL, OP_IFEXPR), 0);
		CHECK_EQ(countNodes(eq, "d", OP_NONE), 0);
		CHECK_EQ(countNodes(eq, "b", OP_NONE), 0);
		CHECK_EQ(countNodes(eq, "a", OP_NONE), 0); // b is true: the `-d*x` branch
		CHECK_EQ(countNodes(eq, "x", OP_NONE), 2); // der(x) and x
		CHECK_EQ(countNodes(eq, NULL, OP_CONST), 1);
	}
	{
		Parsed p;
		parseInto(p, kSwitching);
		HybridEnv env = makeHybridEnv(*p.model);
		env.values["b"] = 0.0;
		string err;
		CHECK(specializeModel(*p.model, env, err));
		const AstNode *eq = p.model->modelTrees.front()->pHead;
		CHECK_EQ(countNodes(eq, "a", OP_NONE), 1);
		CHECK_EQ(countNodes(eq, "d", OP_NONE), 0);
		CHECK_EQ(countNodes(eq, NULL, OP_IFEXPR), 0);
	}
}

static void testIfEquationSpecialization() {
	const char *src =
		"model M\n"
		"  Boolean m(start = true);\n"
		"  Real x(start = 1);\n"
		"  Real y(start = 0);\n"
		"equation\n"
		"  if m then\n"
		"    der(x) = -x;\n"
		"    der(y) = x;\n"
		"  else\n"
		"    der(x) = -2*x;\n"
		"    der(y) = 0;\n"
		"  end if;\n"
		"  when time >= 1 then\n"
		"    m = false;\n"
		"  end when;\n"
		"end M;\n";
	for (int mode = 1; mode >= 0; mode--) {
		Parsed p;
		parseInto(p, src);
		CHECK_EQ(p.parser.errCount, 0);
		HybridEnv env = makeHybridEnv(*p.model);
		env.values["m"] = mode;
		string err;
		CHECK(specializeModel(*p.model, env, err));
		CHECK_EQ(p.model->modelTrees.size(), (size_t)2);
		int muls = 0;
		for (const shared_ptr<AstTree> &t : p.model->modelTrees) {
			CHECK(t->pHead != NULL && t->pHead->op == OP_EQN && t->pHead->nodeKind != NK_IfEquation);
			muls += countNodes(t->pHead, NULL, OP_MUL);
		}
		CHECK_EQ(muls, mode == 1 ? 0 : 1);
	}
}

static void testTimeEvents() {
	Parsed p;
	parseInto(p, kSwitching);
	HybridEnv env = makeHybridEnv(*p.model);
	std::vector<WhenClause> clauses = collectWhenClauses(*p.model);
	CHECK_EQ(clauses.size(), (size_t)2);
	string err;
	double te = 0;

	CHECK(nextEventTime(clauses, env, 0.0, te, err));
	CHECK_NEAR(te, 0.5, 1e-12);

	HybridEvent ev;
	CHECK(fireEvents(clauses, env, 0.5, false, ev, err));
	CHECK_EQ(ev.assignments.size(), (size_t)1);
	CHECK_NEAR(env.values["n"], 1.0, 0.0);

	CHECK(nextEventTime(clauses, env, 0.5, te, err));
	CHECK_NEAR(te, 1.5, 1e-12);
	CHECK(fireEvents(clauses, env, 1.5, false, ev, err)); // both clauses fire
	CHECK_EQ(ev.assignments.size(), (size_t)3);
	CHECK_NEAR(env.values["n"], 2.0, 0.0);
	CHECK_NEAR(env.values["d"], 4.0, 0.0);
	CHECK_NEAR(env.values["b"], 0.0, 0.0);

	CHECK(nextEventTime(clauses, env, 1.5, te, err));
	CHECK_NEAR(te, 2.5, 1e-12); // the time >= 1.5 clause is done
	CHECK(fireEvents(clauses, env, 2.5, false, ev, err));
	CHECK_NEAR(env.values["n"], 3.0, 0.0);
	CHECK_NEAR(env.values["d"], 4.0, 0.0);

	// Nothing triggers between events.
	CHECK(fireEvents(clauses, env, 2.7, false, ev, err));
	CHECK(ev.assignments.empty() && ev.reinits.empty());
}

static void testInitialSampleAndReinit() {
	const char *src =
		"model M\n"
		"  Integer k(start = 0);\n"
		"  Real x(start = 1);\n"
		"equation\n"
		"  der(x) = -x;\n"
		"  when sample(0, 2) then\n"
		"    k = pre(k) + 1;\n"
		"    reinit(x, 2*x + k);\n"
		"  end when;\n"
		"  when time >= 1 then\n"
		"    k = 100;\n"
		"  end when;\n"
		"end M;\n";
	Parsed p;
	parseInto(p, src);
	CHECK_EQ(p.parser.errCount, 0);
	HybridEnv env = makeHybridEnv(*p.model);
	env.values["x"] = 3.0; // continuous value at the event
	std::vector<WhenClause> clauses = collectWhenClauses(*p.model);
	string err;
	HybridEvent ev;
	// At the start time only sample() fires.
	CHECK(fireEvents(clauses, env, 0.0, true, ev, err));
	CHECK_EQ(ev.assignments.size(), (size_t)1);
	CHECK_EQ(ev.reinits.size(), (size_t)1);
	if (ev.reinits.size() == 1) {
		CHECK(ev.reinits[0].first == "x");
		CHECK_NEAR(ev.reinits[0].second, 2 * 3.0 + 0.0, 1e-12); // k is evaluated before the event: pre values
	}
	CHECK_NEAR(env.values["k"], 1.0, 0.0);
	double te = 0;
	CHECK(nextEventTime(clauses, env, 0.0, te, err));
	CHECK_NEAR(te, 1.0, 1e-12);
	CHECK(nextEventTime(clauses, env, 1.0, te, err));
	CHECK_NEAR(te, 2.0, 1e-12);

	// Two branches assigning the same variable at the same time is an error.
	CHECK(fireEvents(clauses, env, 2.0, false, ev, err) || !err.empty());
}

static void testStrictChecker() {
	string text;
	{
		Parsed p;
		parseInto(p, kSwitching);
		CHECK_EQ(check(p, text), 0);
	}
	{
		Parsed p;
		parseInto(p, "model M\n  Real x(start=1);\nequation\n  der(x) = -x;\n  when x < 0 then\n    reinit(x, 1);\n  end when;\nend M;\n");
		CHECK_EQ(p.parser.errCount, 0);
		CHECK(check(p, text) > 0);
		CHECK(contains(text, "when-condition must be `time >= e`"));
	}
	{
		Parsed p;
		parseInto(p, "model M\n  Real x(start=1);\nequation\n  der(x) = if x > 0 then -1 else 1;\nend M;\n");
		CHECK(check(p, text) > 0);
		CHECK(contains(text, "condition depends on continuous variable 'x'"));
	}
	{
		Parsed p;
		parseInto(p, "model M\n  Real x(start=1);\nequation\n  der(x) = -abs(x);\nend M;\n");
		CHECK(check(p, text) > 0);
		CHECK(contains(text, "function abs() is not supported in continuous equations"));
	}
	{
		Parsed p;
		parseInto(p, "model M\n  Real x(start=1);\nequation\n  der(x) = -time*x;\nend M;\n");
		CHECK(check(p, text) > 0);
		CHECK(contains(text, "'time' may only be used in when-conditions"));
	}
	{
		Parsed p;
		parseInto(p, "model M\n  discrete Real d(start=1);\n  Real x(start=1);\nequation\n  der(x) = -x;\n  d = 2;\nend M;\n");
		CHECK(check(p, text) > 0);
		CHECK(contains(text, "can only be assigned in a when-equation"));
	}
	{
		Parsed p;
		parseInto(p, "model M\n  Real x(start=1);\n  Real y(start=2);\nequation\n  der(x) = -x;\n  der(y) = -y;\n  when time >= 1 then\n    y = 3;\n  end when;\nend M;\n");
		CHECK(check(p, text) > 0);
		CHECK(contains(text, "only discrete variables can be assigned"));
	}
	{
		Parsed p;
		parseInto(p, "model M\n  Boolean b(start=true);\n  Real x(start=1);\nequation\n  if b then\n    der(x) = -x;\n  else\n    der(x) = x;\n  end if;\nend M;\n");
		CHECK_EQ(check(p, text), 0);
	}
	{
		Parsed p;
		parseInto(p, "model M\n  Boolean b(start=true);\n  Real x(start=1);\nequation\n  der(x) = -x;\n  if b then\n    when time >= 1 then\n      b = false;\n    end when;\n  end if;\nend M;\n");
		CHECK(check(p, text) > 0);
		CHECK(contains(text, "when-equations inside if-equations are not supported yet"));
	}
}

static const char *kBallModel =
	"model Ball\n"
	"  parameter Real g = 9.81;\n"
	"  Real h(start = 1);\n"
	"  Real v(start = 0);\n"
	"  Boolean up(start = false);\n"
	"equation\n"
	"  der(h) = v;\n"
	"  der(v) = if up then -g else -2*g;\n"
	"  when h < 0 then\n"
	"    reinit(v, -0.5*pre(v));\n"
	"    up = true;\n"
	"  end when;\n"
	"  when time >= 5 then\n"
	"    up = false;\n"
	"  end when;\n"
	"end Ball;\n";

static SolverFeatures stateFeatures() {
	SolverFeatures f = hybridFeatures();
	f.stateEvents = true;
	return f;
}

static void testGuards() {
	Parsed p;
	parseInto(p, kBallModel);
	CHECK_EQ(p.parser.errCount, 0);
	HybridEnv env = makeHybridEnv(*p.model);
	std::vector<Guard> guards = collectGuards(*p.model, env);
	// if-expression on the discrete `up` (static), when h < 0 (dynamic), when time >= 5 (a time event).
	CHECK_EQ(guards.size(), (size_t)3);
	if (guards.size() == 3) {
		CHECK(guards[0].kind == G_IFEXPR && !guards[0].dynamic);
		CHECK(guards[1].kind == G_WHEN && guards[1].dynamic);
		CHECK(guards[2].kind == G_WHEN && !guards[2].dynamic);
	}
	std::vector<WhenClause> clauses = collectWhenClauses(*p.model);
	CHECK_EQ(clauses.size(), (size_t)2);
	if (clauses.size() == 2) {
		CHECK_EQ(clauses[0].links[0].guardId, 1);
		CHECK_EQ(clauses[1].links[0].guardId, 2);
	}

	// evalGuards needs the continuous values.
	env.values["h"] = -0.5;
	env.values["v"] = -3.0;
	std::map<int, bool> truth;
	string err;
	CHECK(evalGuards(guards, env, truth, err));
	CHECK_EQ(truth.size(), (size_t)1);
	CHECK(truth.count(1) == 1 && truth[1]);
	env.values["h"] = 0.5;
	CHECK(evalGuards(guards, env, truth, err));
	CHECK(!truth[1]);

	// A when-branch with a dynamic guard fires only if its guard is in the fired set.
	std::set<int> fired;
	HybridEvent ev;
	CHECK(fireEvents(clauses, env, 1.0, false, ev, err, &fired, true));
	CHECK(ev.assignments.empty() && ev.reinits.empty());
	fired.insert(1);
	env.values["v"] = -4.0;
	CHECK(fireEvents(clauses, env, 1.0, false, ev, err, &fired, true));
	CHECK_EQ(ev.reinits.size(), (size_t)1);
	if (ev.reinits.size() == 1) {
		CHECK_NEAR(ev.reinits[0].second, 2.0, 1e-12);
	}
	CHECK_NEAR(env.values["up"], 1.0, 0.0);
	// Time events are not scheduled for state-event branches.
	double te = 0;
	CHECK(nextEventTime(clauses, env, 0.0, te, err));
	CHECK_NEAR(te, 5.0, 1e-12);
}

static void testDynamicTruthSpecialization() {
	const char *src =
		"model M\n"
		"  Real x(start = 0);\n"
		"equation\n"
		"  der(x) = if x < 1 then 1 else 0;\n"
		"end M;\n";
	for (int truth = 1; truth >= 0; truth--) {
		Parsed p;
		parseInto(p, src);
		HybridEnv env = makeHybridEnv(*p.model);
		std::vector<Guard> guards = collectGuards(*p.model, env);
		CHECK_EQ(guards.size(), (size_t)1);
		CHECK(guards.size() == 1 && guards[0].dynamic);
		std::map<int, bool> dyn;
		dyn[0] = truth != 0;
		string err;
		CHECK(specializeModel(*p.model, env, err, &dyn));
		const AstNode *eq = p.model->modelTrees.front()->pHead;
		CHECK_EQ(countNodes(eq, NULL, OP_IFEXPR), 0);
		const AstNode *rhs = eq->subs[1];
		CHECK(rhs != NULL && rhs->op == OP_CONST);
		if (rhs) {
			CHECK_NEAR(rhs->ty == TY_REAL ? rhs->value.d : (double)rhs->value.i[0], truth ? 1.0 : 0.0, 0.0);
		}
		// The guard ids survive copying through the index field of the if-expression node.
		CHECK_EQ(countNodes(eq, "x", OP_NONE), 1); // only der(x)
	}
}

static void testStateEventChecker() {
	string text;
	{
		Parsed p;
		parseInto(p, kBallModel);
		std::ostringstream os;
		CHECK_EQ(checkModel(*p.model, stateFeatures(), "test.mo", os), 0);
		// Without stateEvents the state-event condition is rejected.
		CHECK(check(p, text) > 0);
		CHECK(contains(text, "when-condition must be `time >= e`"));
	}
	{
		Parsed p;
		parseInto(p, "model M\n  Real x(start=1);\nequation\n  der(x) = if x > 0 then -1 else 1;\nend M;\n");
		std::ostringstream os;
		CHECK_EQ(checkModel(*p.model, stateFeatures(), "test.mo", os), 0);
	}
	{
		Parsed p;
		parseInto(p, "model M\n  Real x(start=1);\nequation\n  der(x) = -abs(x);\nend M;\n");
		std::ostringstream os;
		CHECK(checkModel(*p.model, stateFeatures(), "test.mo", os) > 0); // still not allowed in continuous equations
	}
	{
		Parsed p;
		parseInto(p, "model M\n  Real x(start=1);\n  Boolean b(start=false);\nequation\n  der(x) = -x;\n  when edge(b) then\n    reinit(x, 1);\n  end when;\nend M;\n");
		std::ostringstream os;
		CHECK_EQ(checkModel(*p.model, stateFeatures(), "test.mo", os), 0); // edge() in a when-condition
	}
	const char *bad[] = {
		// outside a when-condition
		"model M\n  Real x(start=1);\n  Boolean b(start=false);\n  Integer n(start=0);\nequation\n  der(x) = -x;\n  when time >= 1 then\n    n = if edge(b) then 1 else 0;\n  end when;\nend M;\n",
		// edge() of a continuous variable
		"model M\n  Real x(start=1);\nequation\n  der(x) = -x;\n  when edge(x) then\n    reinit(x, 1);\n  end when;\nend M;\n",
		// edge() needs a Boolean
		"model M\n  Real x(start=1);\n  Integer k(start=0);\nequation\n  der(x) = -x;\n  when edge(k) then\n    reinit(x, 1);\n  end when;\nend M;\n",
		// an expression instead of a variable
		"model M\n  Real x(start=1);\n  Boolean b(start=false);\nequation\n  der(x) = -x;\n  when change(not b) then\n    reinit(x, 1);\n  end when;\nend M;\n",
	};
	for (const char *src : bad) {
		Parsed p;
		parseInto(p, src);
		std::ostringstream os;
		CHECK(checkModel(*p.model, stateFeatures(), "test.mo", os) > 0);
	}
}

static void testNotHybrid() {
	Parsed p;
	parseInto(p, "model M\n  Real x(start=1);\nequation\n  der(x) = -x;\nend M;\n");
	CHECK(!isHybridModel(*p.model));
	CHECK(collectWhenClauses(*p.model).empty());
	HybridEnv env = makeHybridEnv(*p.model);
	double te = 0;
	string err;
	CHECK(!nextEventTime(collectWhenClauses(*p.model), env, 0.0, te, err));
	CHECK(err.empty());
}

int main() {
	testEnvAndEval();
	testSpecialize();
	testIfEquationSpecialization();
	testTimeEvents();
	testInitialSampleAndReinit();
	testStrictChecker();
	testGuards();
	testDynamicTruthSpecialization();
	testStateEventChecker();
	testNotHybrid();
	return sastest::finish("hybrid_model_test");
}
