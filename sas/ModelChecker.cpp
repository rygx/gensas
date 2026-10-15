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
					static const char *const supported[] = {"der", "sin", "cos", "exp", "sqrt", "inv", NULL};
					static const char *const piecewise[] = {"abs", "min", "max", "sign", "floor", NULL};
					static const char *const events[] = {"pre", "edge", "change", "sample", "noEvent", "reinit", NULL};
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
						walkExpr(p->subs[1], line); // arguments (chain)
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
