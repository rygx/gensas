//
// Copyright (C) 2022, UChicago Argonne, LLC. All rights reserved.
//
// Software Name: Generic Semi-Analytical Simulation Tool (GenSAS)
// By: Argonne National Laboratory
// OPEN SOURCE LICENSE
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
//
// ******************************************************************************************************
// DISCLAIMER
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// ***************************************************************************************************
//

#include "sas/ModelChecker.h"
#include "sas/HybridModel.h"
#include <cstring>
#include <iostream>
#include <set>
#include <utility>

using namespace std;

namespace che {
	namespace core {

		namespace {

			class Checker {
			public:
				Checker(SasModel &model, const SolverFeatures &features, const string &filename, ostream &err) : model(model), features(features), filename(filename), err(err) {}

				int run() {
					for (const shared_ptr<IdPearl> &id : model.idTable) {
						if (id->idTy == ID_DISC && !features.discreteVariables) {
							report(id->line, "discrete variable '" + name(id.get()) + "' (Integer, Boolean or discrete) is not supported by the solver yet");
						}
					}
					for (const shared_ptr<AstTree> &tree : model.modelTrees) {
						walkStatement(tree->pHead, tree->line);
					}
					if (features.timeEventsOnly) {
						for (const shared_ptr<AstTree> &tree : model.modelTrees) {
							strictStatement(tree->pHead, tree->line, true);
						}
					}
					return errors;
				}

			private:
				SasModel &model;
				const SolverFeatures &features;
				const string &filename;
				ostream &err;
				int errors = 0;
				set<pair<int, string>> reported;

				static string name(const IdPearl *id) {
					return id->str == NULL ? string() : string((const char *)id->str);
				}

				void report(int line, const string &msg) {
					if (!reported.insert(make_pair(line, msg)).second) {
						return;
					}
					err << filename << ":" << line << ": error: " << msg << endl;
					errors++;
				}

				// A statement list is linked through `next`.
				void walkStatements(AstNode *p, int line) {
					for (; p != NULL; p = p->next) {
						walkStatement(p, line, false);
					}
				}

				// `followNext` is false for list elements (the list loop iterates).
				void walkStatement(AstNode *p, int line, bool followNext = false) {
					if (p == NULL) {
						return;
					}
					if (p->line > 0) {
						line = p->line;
					}
					if (p->nodeKind == NK_IfEquation) {
						if (!features.ifEquations) {
							report(line, "if-equations are not supported by the solver yet");
						}
						walkExpr(p->subs[0], line);
						if (p->subs[1] != NULL) {
							walkStatements(p->subs[1]->subs[0], line);
							walkStatements(p->subs[1]->subs[1], line);
						}
					} else if (p->nodeKind == NK_WhenEquation) {
						if (!features.whenEquations) {
							report(line, "when-equations are not supported by the solver yet");
						}
						walkExpr(p->subs[0], line);
						if (p->subs[1] != NULL) {
							walkStatements(p->subs[1]->subs[0], line);
							walkStatement(p->subs[1]->subs[1], line); // elsewhen chain
						}
					} else {
						walkExpr(p, line, false);
					}
					if (followNext) {
						walkStatement(p->next, line, true);
					}
				}

				static bool inList(const string &n, const char *const *list) {
					for (; *list != NULL; list++) {
						if (n == *list) {
							return true;
						}
					}
					return false;
				}

				void checkFunction(const string &n, int line) {
					static const char *const supported[] = {"der", "sin", "cos", "exp", "sqrt", "inv", "log", "atan", "asin", "acos", "tan", "sinh", "cosh", "tanh", "log10", "atan2", NULL};
					static const char *const piecewise[] = {"abs", "min", "max", "sign", "floor", NULL};
					static const char *const events[] = {"pre", "edge", "change", "sample", "initial", "noEvent", "reinit", NULL};
					static const char *const terminateAssert[] = {"terminate", "assert", NULL};
					if (inList(n, supported)) {
						return;
					}
					if (inList(n, piecewise)) {
						if (!features.piecewiseFunctions) {
							report(line, "function " + n + "() is not supported by the solver yet");
						}
					} else if (inList(n, events)) {
						if (!features.eventFunctions) {
							report(line, "function " + n + "() is not supported by the solver yet");
						}
					} else if (inList(n, terminateAssert)) {
						if (!features.terminateAndAssert) {
							report(line, "function " + n + "() is not supported by the solver yet");
						}
					} else {
						report(line, "unknown function '" + n + "'");
					}
				}

				// ---- Rules for the time-event simulator (features.timeEventsOnly) ----

				bool inWhenCondition = false; // checking a state-event when-condition (edge/change are allowed)

				enum ExprCtx { CTX_CONT,  // continuous equation
							   CTX_COND,  // evaluated from parameters and discrete variables only
							   CTX_EVENT, // evaluated at an event: also continuous values and `time`
				};

				shared_ptr<IdPearl> lookup(const AstNode *p) {
					const char *n = (const char *)p->value.p;
					return n == NULL ? nullptr : model.lookupIdFromModel(n, (int)strlen(n));
				}

				static bool isRelational(Op op) {
					return op == OP_EQUAL || op == OP_UNEQUAL || op == OP_GREAT || op == OP_LESS || op == OP_GREAT_EQ || op == OP_LESS_EQ || op == OP_AND || op == OP_OR || op == OP_NOT;
				}

				void strictExpr(const AstNode *p, int line, ExprCtx ctx) {
					if (p == NULL) {
						return;
					}
					if (p->op == OP_ID && p->ty != TY_FUNCTION && p->value.p != NULL) {
						string n = (const char *)p->value.p;
						if (n == "time") {
							if (ctx != CTX_EVENT) {
								report(line, "'time' may only be used in when-conditions and when-bodies");
							}
							return;
						}
						shared_ptr<IdPearl> id = lookup(p);
						if (id != nullptr && ctx == CTX_COND && (id->idTy == ID_VAR || id->idTy == ID_IVAR || id->idTy == ID_TMPVAR)) {
							report(line, "condition depends on continuous variable '" + n + "'; state events are not supported yet");
						}
						return;
					}
					if (p->op == OP_IFEXPR) {
						strictExpr(p->subs[0], line, (ctx == CTX_EVENT || features.stateEvents) ? CTX_EVENT : CTX_COND);
						if (p->subs[1] != NULL) {
							strictExpr(p->subs[1]->subs[0], line, ctx);
							strictExpr(p->subs[1]->subs[1], line, ctx);
						}
						return;
					}
					if (isRelational(p->op) && ctx == CTX_CONT) {
						report(line, "relational and logical operators can only be used in conditions");
					}
					if (p->op == OP_CALL) {
						const AstNode *callee = p->subs[0];
						string fn = (callee != NULL && callee->value.p != NULL) ? (const char *)callee->value.p : "";
						if (fn == "der") {
							if (ctx != CTX_CONT) {
								report(line, "der() cannot be used in conditions or event bodies");
							}
						} else if (fn == "abs" || fn == "min" || fn == "max" || fn == "sign" || fn == "floor" || fn == "pre" || fn == "noEvent") {
							if (ctx == CTX_CONT) {
								report(line, "function " + fn + "() is not supported in continuous equations yet");
							}
						} else if (fn == "edge" || fn == "change") {
							const AstNode *arg = p->subs[1];
							shared_ptr<IdPearl> id = (arg != NULL && arg->op == OP_ID && arg->ty != TY_FUNCTION) ? lookup(arg) : nullptr;
							if (!features.stateEvents) {
								report(line, "function " + fn + "() requires state events, which are not supported yet");
							} else if (!inWhenCondition) {
								report(line, "function " + fn + "() can only be used in when-conditions");
							} else if (arg == NULL || arg->next != NULL || id == nullptr || id->idTy != ID_DISC) {
								report(line, "the argument of " + fn + "() must be one discrete variable");
							} else if (fn == "edge" && id->ty != TY_BOOL) {
								report(line, "the argument of edge() must be a Boolean variable");
							}
							return;
						} else if (fn == "sample") {
							report(line, "sample() can only be used as a when-condition");
						} else if (fn == "initial") {
							report(line, "initial() can only be used as the whole condition of a when-equation");
						}
						for (const AstNode *a = p->subs[1]; a != NULL; a = a->next) {
							strictExpr(a, line, ctx);
						}
						return;
					}
					strictExpr(p->subs[0], line, ctx);
					strictExpr(p->subs[1], line, ctx);
				}

				// assert(cond, "message" [, level]) and terminate("message") statements.
				void strictAssertOrTerminate(const AstNode *s, int line, bool allowTerminate) {
					string e;
					if (isAssertCall(s)) {
						AssertCall ac;
						if (!parseAssertCall(s, ac, e)) {
							report(line, e);
						} else {
							strictExpr(ac.cond, line, CTX_EVENT);
						}
					} else if (!allowTerminate) {
						report(line, "terminate() can only be used in a when-equation");
					} else {
						string message;
						if (!parseTerminateCall(s, message, e)) {
							report(line, e);
						}
					}
				}

				void strictWhen(const AstNode *w, int line) {
					for (; w != NULL && w->nodeKind == NK_WhenEquation; w = w->subs[1] != NULL ? w->subs[1]->subs[1] : NULL) {
						if (w->line > 0) {
							line = w->line;
						}
						TimeCond tc;
						if (!classifyTimeCondition(w->subs[0], tc)) {
							if (!features.stateEvents) {
								report(line, "when-condition must be `time >= e` or `sample(start, interval)`; state events are not supported yet");
							} else {
								inWhenCondition = true;
								strictExpr(w->subs[0], line, CTX_EVENT);
								inWhenCondition = false;
							}
						} else {
							strictExpr(tc.a, line, CTX_COND);
							strictExpr(tc.b, line, CTX_COND);
						}
						if (w->subs[1] == NULL) {
							continue;
						}
						for (const AstNode *s = w->subs[1]->subs[0]; s != NULL; s = s->next) {
							if (s->op == OP_EQN && s->subs[0] != NULL && s->subs[0]->op == OP_ID) {
								shared_ptr<IdPearl> id = lookup(s->subs[0]);
								if (id != nullptr && id->idTy != ID_DISC) {
									report(line, "only discrete variables can be assigned in a when-equation; use reinit() for '" + string((const char *)s->subs[0]->value.p) + "'");
								}
								strictExpr(s->subs[1], line, CTX_EVENT);
							} else if (isAssertCall(s) || isTerminateCall(s)) {
								strictAssertOrTerminate(s, line, true);
							} else if (s->op == OP_CALL && s->subs[0] != NULL && s->subs[0]->value.p != NULL && string((const char *)s->subs[0]->value.p) == "reinit") {
								if (s->subs[1] != NULL && s->subs[1]->next != NULL) {
									strictExpr(s->subs[1]->next, line, CTX_EVENT);
								}
							} else {
								report(line, "unsupported statement in when-equation; expect `d = expr;`, `reinit(x, expr);`, `assert(...)` or `terminate(...)`");
							}
						}
					}
				}

				// `top` is true for statements outside any if-equation.
				void strictStatement(const AstNode *p, int line, bool top) {
					for (; p != NULL; p = p->next) {
						if (p->line > 0) {
							line = p->line;
						}
						if (p->nodeKind == NK_WhenEquation) {
							if (!top) {
								report(line, "when-equations inside if-equations are not supported yet");
							} else {
								strictWhen(p, line);
							}
						} else if (p->nodeKind == NK_IfEquation) {
							strictExpr(p->subs[0], line, features.stateEvents ? CTX_EVENT : CTX_COND);
							if (p->subs[1] != NULL) {
								strictStatement(p->subs[1]->subs[0], line, false);
								strictStatement(p->subs[1]->subs[1], line, false);
							}
						} else if (p->op == OP_EQN) {
							const AstNode *lhs = p->subs[0];
							if (lhs != NULL && lhs->op == OP_ID && lhs->ty != TY_FUNCTION) {
								shared_ptr<IdPearl> id = lookup(lhs);
								if (id != nullptr && id->idTy == ID_DISC) {
									report(line, "discrete variable '" + string((const char *)lhs->value.p) + "' can only be assigned in a when-equation");
									continue;
								}
							}
							strictExpr(lhs, line, CTX_CONT);
							strictExpr(p->subs[1], line, CTX_CONT);
						} else if (isAssertCall(p) || isTerminateCall(p)) {
							if (!top) {
								report(line, "assert() and terminate() are not supported inside if-equations");
							} else {
								strictAssertOrTerminate(p, line, false);
							}
						}
					}
				}

				// Expression nodes: subs[0], subs[1], and `next` (call argument chain).
				void walkExpr(AstNode *p, int line, bool followNext = true) {
					if (p == NULL) {
						return;
					}
					if (p->op == OP_IFEXPR && !features.ifExpressions) {
						report(line, "if-expressions are not supported by the solver yet");
					}
					if (p->op == OP_CALL) {
						AstNode *callee = p->subs[0];
						if (callee != NULL && callee->op == OP_ID && callee->value.p != NULL) {
							checkFunction((const char *)callee->value.p, line);
						}
						if (callee != NULL && callee->op == OP_ID && callee->value.p != NULL && string((const char *)callee->value.p) == "assert") {
							// the level argument AssertionLevel.warning is not an identifier of the model
							for (AstNode *a = p->subs[1]; a != NULL && a->op != OP_MEMBER; a = a->next) {
								walkExpr(a, line, false);
							}
						} else {
							walkExpr(p->subs[1], line); // arguments (chain)
						}
					} else if (p->op == OP_ID && p->ty != TY_FUNCTION && p->value.p != NULL) {
						const char *n = (const char *)p->value.p;
						shared_ptr<IdPearl> id = model.lookupIdFromModel(n, (int)strlen(n));
						if (id == nullptr || id->idTy == ID_NONE) {
							if (string(n) == "time") {
								if (!features.simulationTime) {
									report(line, "'time' is not supported by the solver yet");
								}
							} else {
								report(line, "undeclared identifier '" + string(n) + "'");
							}
						}
					} else {
						walkExpr(p->subs[0], line);
						walkExpr(p->subs[1], line);
					}
					if (followNext) {
						walkExpr(p->next, line);
					}
				}
			};

		} // namespace

		int checkModel(SasModel &model, const SolverFeatures &features, const string &filename, ostream &err) {
			return Checker(model, features, filename, err).run();
		}
	} // namespace core
} // namespace che
