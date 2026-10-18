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
#include "sas/HybridModel.h"
#include "sas/MathFunctions.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <memory>
#include <sstream>

using namespace std;

namespace che {
	namespace core {

		namespace {

			string idName(const AstNode *p) {
				return (p != NULL && p->value.p != NULL) ? string((const char *)p->value.p) : string();
			}

			bool isCall(const AstNode *p, const char *name) {
				return p != NULL && p->op == OP_CALL && p->subs[0] != NULL && p->subs[0]->op == OP_ID && idName(p->subs[0]) == name;
			}

			bool isTimeId(const AstNode *p) {
				return p != NULL && p->op == OP_ID && p->ty != TY_FUNCTION && idName(p) == "time";
			}

			double tolOf(double t) {
				return 1e-9 * max(1.0, fabs(t));
			}

			void freeTree(AstNode *p) {
				if (p == NULL) {
					return;
				}
				freeTree(p->subs[0]);
				freeTree(p->subs[1]);
				freeTree(p->next);
				delete p;
			}

			string lineMsg(int line, const string &msg) {
				ostringstream os;
				if (line > 0) {
					os << "line " << line << ": ";
				}
				os << msg;
				return os.str();
			}
		} // namespace

		double HybridEnv::coerce(const string &name, double v) const {
			map<string, Type>::const_iterator it = types.find(name);
			if (it == types.end()) {
				return v;
			}
			if (it->second == TY_BOOL) {
				return v != 0.0 ? 1.0 : 0.0;
			}
			if (it->second == TY_INT) {
				return v < 0 ? ceil(v) : floor(v); // truncate toward zero
			}
			return v;
		}

		bool HybridEnv::eval(const AstNode *p, double &out, string &err) const {
			if (p == NULL) {
				err = "missing expression";
				return false;
			}
			switch (p->op) {
			case OP_CONST:
				out = (p->ty == TY_REAL) ? p->value.d : (double)p->value.i[0];
				return true;
			case OP_ID: {
				string n = idName(p);
				if (p->ty == TY_FUNCTION) {
					err = "function '" + n + "' used as a value";
					return false;
				}
				map<string, double>::const_iterator it = values.find(n);
				if (it != values.end()) {
					out = it->second;
					return true;
				}
				if (n == "time") {
					out = time;
					return true;
				}
				err = "no value for '" + n + "'";
				return false;
			}
			case OP_IFEXPR: {
				double c;
				if (!eval(p->subs[0], c, err)) {
					return false;
				}
				const AstNode *br = p->subs[1];
				if (br == NULL) {
					err = "malformed if-expression";
					return false;
				}
				return eval(c != 0.0 ? br->subs[0] : br->subs[1], out, err);
			}
			case OP_MINUS:
			case OP_POS:
			case OP_NOT: {
				double a;
				if (!eval(p->subs[0], a, err)) {
					return false;
				}
				out = (p->op == OP_MINUS) ? -a : (p->op == OP_POS ? a : (a == 0.0 ? 1.0 : 0.0));
				return true;
			}
			case OP_ADD:
			case OP_SUB:
			case OP_MUL:
			case OP_DIV:
			case OP_POW:
			case OP_EQUAL:
			case OP_UNEQUAL:
			case OP_GREAT:
			case OP_LESS:
			case OP_GREAT_EQ:
			case OP_LESS_EQ:
			case OP_AND:
			case OP_OR: {
				double a, b;
				if (!eval(p->subs[0], a, err) || !eval(p->subs[1], b, err)) {
					return false;
				}
				switch (p->op) {
				case OP_ADD:
					out = a + b;
					break;
				case OP_SUB:
					out = a - b;
					break;
				case OP_MUL:
					out = a * b;
					break;
				case OP_DIV:
					out = a / b;
					break;
				case OP_POW:
					out = pow(a, b);
					break;
				case OP_EQUAL:
					out = (a == b) ? 1.0 : 0.0;
					break;
				case OP_UNEQUAL:
					out = (a != b) ? 1.0 : 0.0;
					break;
				case OP_GREAT:
					out = (a > b) ? 1.0 : 0.0;
					break;
				case OP_LESS:
					out = (a < b) ? 1.0 : 0.0;
					break;
				case OP_GREAT_EQ:
					out = (a >= b) ? 1.0 : 0.0;
					break;
				case OP_LESS_EQ:
					out = (a <= b) ? 1.0 : 0.0;
					break;
				case OP_AND:
					out = (a != 0.0 && b != 0.0) ? 1.0 : 0.0;
					break;
				default:
					out = (a != 0.0 || b != 0.0) ? 1.0 : 0.0;
					break;
				}
				return true;
			}
			case OP_CALL: {
				string fn = idName(p->subs[0]);
				if ((fn == "edge" || fn == "change") && p->subs[1] != NULL && p->subs[1]->op == OP_ID && p->subs[1]->next == NULL) {
					string n = idName(p->subs[1]);
					map<string, double>::const_iterator now = values.find(n);
					map<string, double>::const_iterator before = eventPre.find(n);
					if (now == values.end()) {
						err = "no value for '" + n + "'";
						return false;
					}
					if (!atEvent || before == eventPre.end()) {
						out = 0.0;
					} else if (fn == "edge") {
						out = (now->second != 0.0 && before->second == 0.0) ? 1.0 : 0.0;
					} else {
						out = (now->second != before->second) ? 1.0 : 0.0;
					}
					return true;
				}
				vector<double> args;
				for (const AstNode *a = p->subs[1]; a != NULL; a = a->next) {
					double v;
					if (fn == "pre" && a->op == OP_ID) {
						string n = idName(a);
						map<string, double>::const_iterator it = pre.find(n);
						if (it == pre.end()) {
							it = values.find(n);
							if (it == values.end()) {
								err = "no value for '" + n + "'";
								return false;
							}
						}
						v = it->second;
					} else if (!eval(a, v, err)) {
						return false;
					}
					args.push_back(v);
				}
				size_t n = args.size();
				if (n == 1 && evalMathFunction(fn, args[0], out)) {
					return true;
				}
				if ((fn == "abs" || fn == "sign" || fn == "floor" || fn == "pre" || fn == "noEvent") && n == 1) {
					double a = args[0];
					if (fn == "abs") {
						out = fabs(a);
					} else if (fn == "sign") {
						out = (a > 0) ? 1.0 : (a < 0 ? -1.0 : 0.0);
					} else if (fn == "floor") {
						out = floor(a);
					} else {
						out = a; // pre, noEvent
					}
					return true;
				}
				if (fn == "atan2" && n == 2) {
					out = atan2(args[0], args[1]);
					return true;
				}
				if ((fn == "min" || fn == "max") && n == 2) {
					out = (fn == "min") ? min(args[0], args[1]) : max(args[0], args[1]);
					return true;
				}
				err = "function '" + fn + "' cannot be evaluated here";
				return false;
			}
			default:
				err = "expression cannot be evaluated here";
				return false;
			}
		}

		HybridEnv makeHybridEnv(const SasModel &model) {
			HybridEnv env;
			for (const shared_ptr<IdPearl> &id : model.idTable) {
				if (id->str == NULL) {
					continue;
				}
				string n = (const char *)id->str;
				if (id->idTy == ID_PARAM) {
					env.values[n] = (id->ty == TY_INT || id->ty == TY_BOOL) ? (double)id->value.i : id->value.d;
				} else if (id->idTy == ID_DISC) {
					double v = 0.0;
					if (id->initVal) {
						v = (id->ty == TY_INT || id->ty == TY_BOOL) ? (double)id->value.i : id->value.d;
					}
					env.discrete.insert(n);
					env.types[n] = id->ty;
					env.values[n] = v;
				}
			}
			return env;
		}

		bool isAssertCall(const AstNode *s) {
			return isCall(s, "assert");
		}

		bool isTerminateCall(const AstNode *s) {
			return isCall(s, "terminate");
		}

		bool parseAssertCall(const AstNode *call, AssertCall &out, string &err) {
			const AstNode *cond = call != NULL ? call->subs[1] : NULL;
			const AstNode *msg = cond != NULL ? cond->next : NULL;
			const AstNode *level = msg != NULL ? msg->next : NULL;
			if (cond == NULL || msg == NULL || msg->op != OP_STR || (level != NULL && level->next != NULL)) {
				err = "assert() needs a condition and a message string, assert(condition, \"message\")";
				return false;
			}
			out.cond = cond;
			out.message = idName(msg);
			out.warning = false;
			if (level != NULL) {
				const string name = level->op == OP_MEMBER && level->subs[0] != NULL ? idName(level->subs[0]) + "." + idName(level) : string();
				if (name == "AssertionLevel.warning") {
					out.warning = true;
				} else if (name != "AssertionLevel.error") {
					err = "the level of assert() must be AssertionLevel.error or AssertionLevel.warning";
					return false;
				}
			}
			return true;
		}

		bool parseTerminateCall(const AstNode *call, string &message, string &err) {
			const AstNode *msg = call != NULL ? call->subs[1] : NULL;
			if (msg == NULL || msg->op != OP_STR || msg->next != NULL) {
				err = "terminate() needs one message string, terminate(\"message\")";
				return false;
			}
			message = idName(msg);
			return true;
		}

		namespace {
			bool hasHybridNode(const AstNode *p) {
				for (; p != NULL; p = p->next) {
					if (p->nodeKind == NK_IfEquation || p->nodeKind == NK_WhenEquation || p->op == OP_IFEXPR || isTimeId(p) || isAssertCall(p) || isTerminateCall(p)) {
						return true;
					}
					if (hasHybridNode(p->subs[0]) || hasHybridNode(p->subs[1])) {
						return true;
					}
				}
				return false;
			}
		} // namespace

		namespace {
			AstNode *cloneChain(const AstNode *p);

			// Deep copy of one expression node (without its `next` sibling).
			AstNode *cloneNode(const AstNode *p) {
				if (p == NULL) {
					return NULL;
				}
				AstNode *c = new AstNode(*p);
				c->next = NULL;
				c->index = -1;
				c->subs[0] = cloneChain(p->subs[0]);
				c->subs[1] = cloneChain(p->subs[1]);
				return c;
			}

			AstNode *cloneChain(const AstNode *p) {
				AstNode *head = NULL;
				AstNode **tail = &head;
				for (; p != NULL; p = p->next) {
					*tail = cloneNode(p);
					tail = &(*tail)->next;
				}
				return head;
			}

			AstNode *mkNode(Op op, AstNode *a, AstNode *b = NULL) {
				AstNode *n = new AstNode();
				n->op = op;
				n->nodeKind = NK_Expression;
				n->subs[0] = a;
				n->subs[1] = b;
				return n;
			}

			AstNode *mkConst(double v) {
				AstNode *n = mkNode(OP_CONST, NULL);
				n->ty = TY_REAL;
				n->value.d = v;
				return n;
			}

			// if c then a else b
			AstNode *mkIf(AstNode *c, AstNode *a, AstNode *b) {
				return mkNode(OP_IFEXPR, c, mkNode(OP_BRANCH, a, b));
			}

			struct RewriteCtx {
				PiecewiseMode mode = PW_EVENTS;
				double eps = 1e-6;
				int count = 0;
				SasModel *model = NULL;
			};

			// sqrt(e*e + eps*eps); takes ownership of e.
			AstNode *smoothAbs(AstNode *e, RewriteCtx &ctx) {
				const double eps = ctx.eps;
				AstNode *square = mkNode(OP_MUL, e, cloneNode(e));
				AstNode *radicand = mkNode(OP_ADD, square, mkConst(eps * eps));
				AstNode *callee = mkNode(OP_ID, NULL);
				callee->ty = TY_FUNCTION;
				// The model's identifier table must know the function (the solver resolves names there).
				callee->value.p = (decltype(callee->value.p))ctx.model->lookupIdentifier("sqrt", 4)->str;
				return mkNode(OP_CALL, callee, radicand);
			}

			// A call of the function `name` with the single argument arg (takes ownership of arg).
			AstNode *mkCall(const string &name, AstNode *arg, RewriteCtx &ctx) {
				AstNode *callee = mkNode(OP_ID, NULL);
				callee->ty = TY_FUNCTION;
				callee->value.p = (decltype(callee->value.p))ctx.model->lookupIdentifier(name.c_str(), (int)name.size())->str;
				return mkNode(OP_CALL, callee, arg);
			}

			// tan, sinh, cosh, tanh and log10 in terms of the functions the solver knows; NULL for other calls.
			AstNode *rewriteDerivedFunction(const string &fn, const AstNode *a, RewriteCtx &ctx) {
				if (fn == "tan") {
					return mkNode(OP_DIV, mkCall("sin", cloneNode(a), ctx), mkCall("cos", cloneNode(a), ctx));
				}
				if (fn == "sinh" || fn == "cosh" || fn == "tanh") {
					AstNode *up = mkCall("exp", cloneNode(a), ctx);
					AstNode *down = mkCall("exp", mkNode(OP_MINUS, cloneNode(a)), ctx);
					if (fn == "tanh") {
						return mkNode(OP_DIV, mkNode(OP_SUB, up, down), mkNode(OP_ADD, cloneNode(up), cloneNode(down)));
					}
					return mkNode(OP_MUL, mkConst(0.5), mkNode(fn == "sinh" ? OP_SUB : OP_ADD, up, down));
				}
				if (fn == "log10") {
					return mkNode(OP_MUL, mkConst(0.4342944819032518), mkCall("log", cloneNode(a), ctx));
				}
				return NULL;
			}

			// The value of a constant expression (numbers and parameters combined with + - * / and unary minus).
			bool constantNumber(const AstNode *p, RewriteCtx &ctx, double &out) {
				if (p == NULL) {
					return false;
				}
				double a, b;
				switch (p->op) {
				case OP_CONST:
					out = p->ty == TY_INT ? (double)p->value.i[0] : p->value.d;
					return true;
				case OP_ID: {
					if (p->ty == TY_FUNCTION) {
						return false;
					}
					const char *n = (const char *)p->value.p;
					shared_ptr<IdPearl> id = ctx.model->lookupIdFromModel((char *)n, (int)strlen(n));
					if (id == nullptr || id->idTy != ID_PARAM || (id->ty != TY_INT && id->ty != TY_REAL)) {
						return false;
					}
					out = id->ty == TY_INT ? (double)id->value.i : id->value.d;
					return true;
				}
				case OP_MINUS:
					if (!constantNumber(p->subs[0], ctx, a)) {
						return false;
					}
					out = -a;
					return true;
				case OP_POS:
					return constantNumber(p->subs[0], ctx, out);
				case OP_ADD:
				case OP_SUB:
				case OP_MUL:
				case OP_DIV:
					if (!constantNumber(p->subs[0], ctx, a) || !constantNumber(p->subs[1], ctx, b)) {
						return false;
					}
					if (p->op == OP_DIV && b == 0.0) {
						return false;
					}
					if (p->op == OP_ADD) {
						out = a + b;
					} else if (p->op == OP_SUB) {
						out = a - b;
					} else if (p->op == OP_MUL) {
						out = a * b;
					} else {
						out = a / b;
					}
					return true;
				default:
					return false;
				}
			}

			// x^y in terms of products and the functions the solver knows: a small integer power is a
			// product (1/product for a negative one), x^0.5 is sqrt(x), and any other power is
			// exp(y*log(x)), which needs x > 0.
			AstNode *rewritePow(const AstNode *x, const AstNode *y, RewriteCtx &ctx) {
				const int maxIntPower = 16;
				double cx, n;
				if (constantNumber(x, ctx, cx) && constantNumber(y, ctx, n) && std::isfinite(std::pow(cx, n))) {
					return mkConst(std::pow(cx, n));
				}
				if (constantNumber(y, ctx, n)) {
					if (n == std::floor(n) && std::fabs(n) <= maxIntPower) {
						if (n == 0.0) {
							return mkConst(1.0);
						}
						AstNode *prod = cloneNode(x);
						for (int i = 1; i < (int)std::fabs(n); i++) {
							prod = mkNode(OP_MUL, prod, cloneNode(x));
						}
						return n < 0.0 ? mkNode(OP_DIV, mkConst(1.0), prod) : prod;
					}
					if (n == 0.5) {
						return mkCall("sqrt", cloneNode(x), ctx);
					}
				}
				double bx;
				AstNode *logx = (constantNumber(x, ctx, bx) && bx > 0.0) ? mkConst(std::log(bx)) : mkCall("log", cloneNode(x), ctx);
				return mkCall("exp", mkNode(OP_MUL, cloneNode(y), logx), ctx);
			}

			// atan2(y, x) with the half-angle formula atan2(y, x) = 2 atan(y / (sqrt(x^2 + y^2) + x)), which is smooth
			// for x >= 0; on the left half plane it is pi - atan2(y, -x) (y >= 0) or -pi - atan2(y, -x) (y < 0).
			// The switches at x = 0 (continuous) and y = 0 with x < 0 (the branch cut) are state events.
			AstNode *rewriteAtan2(const AstNode *y, const AstNode *x, RewriteCtx &ctx) {
				auto radius = [&]() {
					AstNode *sum = mkNode(OP_ADD, mkNode(OP_MUL, cloneNode(x), cloneNode(x)), mkNode(OP_MUL, cloneNode(y), cloneNode(y)));
					return mkCall("sqrt", sum, ctx);
				};
				auto half = [&](Op op) { // 2 atan(y / (r +- x))
					AstNode *den = mkNode(op, radius(), cloneNode(x));
					return mkNode(OP_MUL, mkConst(2.0), mkCall("atan", mkNode(OP_DIV, cloneNode(y), den), ctx));
				};
				const double pi = 3.14159265358979323846;
				AstNode *upper = mkNode(OP_SUB, mkConst(pi), half(OP_SUB));
				AstNode *lower = mkNode(OP_SUB, mkConst(-pi), half(OP_SUB));
				AstNode *left = mkIf(mkNode(OP_GREAT_EQ, cloneNode(y), mkConst(0.0)), upper, lower);
				return mkIf(mkNode(OP_GREAT_EQ, cloneNode(x), mkConst(0.0)), half(OP_ADD), left);
			}

			// Rewrites the expression chain starting at p; returns the new head.
			AstNode *rewriteChain(AstNode *p, RewriteCtx &ctx);

			AstNode *rewriteOne(AstNode *p, RewriteCtx &ctx) {
				if (p->op == OP_IFEXPR) {
					// The condition is evaluated directly; only the branches are equations.
					if (p->subs[1] != NULL) {
						p->subs[1]->subs[0] = rewriteChain(p->subs[1]->subs[0], ctx);
						p->subs[1]->subs[1] = rewriteChain(p->subs[1]->subs[1], ctx);
					}
					return p;
				}
				p->subs[0] = rewriteChain(p->subs[0], ctx);
				p->subs[1] = rewriteChain(p->subs[1], ctx);
				if (p->op == OP_POW && p->subs[0] != NULL && p->subs[1] != NULL) {
					AstNode *pw = rewritePow(p->subs[0], p->subs[1], ctx);
					ctx.count++;
					pw->next = p->next;
					pw->line = p->line;
					p->next = NULL;
					freeTree(p);
					return pw;
				}
				if (p->op != OP_CALL || p->subs[0] == NULL || p->subs[0]->op != OP_ID) {
					return p;
				}
				string fn = idName(p->subs[0]);
				const AstNode *a = p->subs[1];
				const AstNode *b = a != NULL ? a->next : NULL;
				AstNode *repl = NULL;
				// A function of constants is a constant; the solver series need a variable as the argument.
				double ca, cb, cv;
				if (a != NULL && b == NULL && constantNumber(a, ctx, ca) && evalMathFunction(fn, ca, cv) && std::isfinite(cv)) {
					repl = mkConst(cv);
				} else if (fn == "atan2" && a != NULL && b != NULL && b->next == NULL && constantNumber(a, ctx, ca) && constantNumber(b, ctx, cb)) {
					repl = mkConst(std::atan2(ca, cb));
				}
				if (repl == NULL && a != NULL && b == NULL) {
					repl = rewriteDerivedFunction(fn, a, ctx);
				}
				if (repl == NULL && fn == "atan2" && a != NULL && b != NULL && b->next == NULL) {
					repl = rewriteAtan2(a, b, ctx);
				}
				if (repl != NULL) {
					// done
				} else if (fn == "noEvent" && a != NULL && b == NULL) {
					// noEvent(e) only asks the tool not to generate events for e. Keeping the events is
					// always allowed and more accurate, so the call is dropped.
					repl = cloneNode(a);
				} else if (ctx.mode == PW_SMOOTH) {
					if ((fn == "abs" || fn == "sign") && a != NULL && b == NULL) {
						repl = fn == "abs" ? smoothAbs(cloneNode(a), ctx) : mkNode(OP_DIV, cloneNode(a), smoothAbs(cloneNode(a), ctx));
					} else if ((fn == "min" || fn == "max") && a != NULL && b != NULL && b->next == NULL) {
						AstNode *sum = mkNode(OP_ADD, cloneNode(a), cloneNode(b));
						AstNode *dev = smoothAbs(mkNode(OP_SUB, cloneNode(a), cloneNode(b)), ctx);
						repl = mkNode(OP_MUL, mkConst(0.5), mkNode(fn == "max" ? OP_ADD : OP_SUB, sum, dev));
					}
				} else if ((fn == "abs" || fn == "sign") && a != NULL && b == NULL) {
					if (fn == "abs") {
						repl = mkIf(mkNode(OP_GREAT_EQ, cloneNode(a), mkConst(0.0)), cloneNode(a), mkNode(OP_MINUS, cloneNode(a)));
					} else {
						repl = mkIf(mkNode(OP_GREAT, cloneNode(a), mkConst(0.0)), mkConst(1.0),
									mkIf(mkNode(OP_LESS, cloneNode(a), mkConst(0.0)), mkConst(-1.0), mkConst(0.0)));
					}
				} else if ((fn == "min" || fn == "max") && a != NULL && b != NULL && b->next == NULL) {
					repl = mkIf(mkNode(fn == "max" ? OP_GREAT_EQ : OP_LESS_EQ, cloneNode(a), cloneNode(b)), cloneNode(a), cloneNode(b));
				}
				if (repl == NULL) {
					return p;
				}
				ctx.count++;
				repl->next = p->next;
				repl->line = p->line;
				p->next = NULL;
				freeTree(p);
				return repl;
			}

			AstNode *rewriteChain(AstNode *p, RewriteCtx &ctx) {
				AstNode *head = NULL;
				AstNode **tail = &head;
				while (p != NULL) {
					AstNode *nxt = p->next;
					p->next = NULL;
					AstNode *r = rewriteOne(p, ctx);
					AstNode *last = r;
					while (last->next != NULL) {
						last = last->next;
					}
					*tail = r;
					tail = &last->next;
					p = nxt;
				}
				return head;
			}

			// Statements: equations are rewritten, if-equation branches are searched, when-equations skipped.
			AstNode *rewriteStatements(AstNode *p, RewriteCtx &ctx) {
				for (AstNode *s = p; s != NULL; s = s->next) {
					if (s->nodeKind == NK_WhenEquation) {
						continue;
					}
					if (s->nodeKind == NK_IfEquation) {
						if (s->subs[1] != NULL) {
							s->subs[1]->subs[0] = rewriteStatements(s->subs[1]->subs[0], ctx);
							s->subs[1]->subs[1] = rewriteStatements(s->subs[1]->subs[1], ctx);
						}
					} else if (s->op == OP_EQN) {
						AstNode *nxt = s->next;
						s->next = NULL;
						s->subs[0] = rewriteChain(s->subs[0], ctx);
						s->subs[1] = rewriteChain(s->subs[1], ctx);
						s->next = nxt;
					}
				}
				return p;
			}
		} // namespace

		int rewritePiecewiseFunctions(SasModel &model, PiecewiseMode mode, double eps) {
			RewriteCtx ctx;
			ctx.mode = mode;
			ctx.eps = eps;
			ctx.model = &model;
			for (const shared_ptr<AstTree> &tree : model.modelTrees) {
				tree->pHead = rewriteStatements(tree->pHead, ctx);
			}
			return ctx.count;
		}

		bool isHybridModel(const SasModel &model) {
			for (const shared_ptr<IdPearl> &id : model.idTable) {
				if (id->idTy == ID_DISC) {
					return true;
				}
			}
			for (const shared_ptr<AstTree> &tree : model.modelTrees) {
				if (hasHybridNode(tree->pHead)) {
					return true;
				}
			}
			return false;
		}

		bool classifyTimeCondition(const AstNode *cond, TimeCond &tc) {
			tc = TimeCond();
			if (cond == NULL) {
				return false;
			}
			if ((cond->op == OP_GREAT_EQ || cond->op == OP_GREAT) && isTimeId(cond->subs[0]) && cond->subs[1] != NULL) {
				tc.kind = TC_AT;
				tc.a = cond->subs[1];
				return true;
			}
			if ((cond->op == OP_LESS_EQ || cond->op == OP_LESS) && isTimeId(cond->subs[1]) && cond->subs[0] != NULL) {
				tc.kind = TC_AT;
				tc.a = cond->subs[0];
				return true;
			}
			if (isCall(cond, "initial") && cond->subs[1] == NULL) {
				tc.kind = TC_INITIAL;
				return true;
			}
			if (isCall(cond, "sample") && cond->subs[1] != NULL && cond->subs[1]->next != NULL && cond->subs[1]->next->next == NULL) {
				tc.kind = TC_SAMPLE;
				tc.a = cond->subs[1];
				tc.b = cond->subs[1]->next;
				return true;
			}
			return false;
		}

		vector<WhenClause> collectWhenClauses(const SasModel &model) {
			vector<WhenClause> clauses;
			for (const shared_ptr<AstTree> &tree : model.modelTrees) {
				const AstNode *w = tree->pHead;
				if (w == NULL || w->nodeKind != NK_WhenEquation) {
					continue;
				}
				WhenClause clause;
				while (w != NULL && w->nodeKind == NK_WhenEquation) {
					WhenLink link;
					link.cond = w->subs[0];
					link.line = w->line;
					link.guardId = w->index;
					const AstNode *br = w->subs[1];
					link.body = br != NULL ? br->subs[0] : NULL;
					clause.links.push_back(link);
					w = br != NULL ? br->subs[1] : NULL;
				}
				clauses.push_back(clause);
			}
			return clauses;
		}

		namespace {

			// Does the expression depend on `time` or on a variable the environment does not know
			// (that is, a continuous variable)?
			bool refsDynamic(const AstNode *p, const HybridEnv &env) {
				for (; p != NULL; p = p->next) {
					if (p->op == OP_ID && p->ty != TY_FUNCTION) {
						string n = idName(p);
						if (n == "time" || env.values.find(n) == env.values.end()) {
							return true;
						}
					}
					if (refsDynamic(p->subs[0], env) || refsDynamic(p->subs[1], env)) {
						return true;
					}
				}
				return false;
			}

			void numberExprGuards(AstNode *p, const HybridEnv &env, vector<Guard> &out) {
				for (; p != NULL; p = p->next) {
					if (p->op == OP_IFEXPR) {
						Guard g;
						g.id = (int)out.size();
						g.kind = G_IFEXPR;
						g.cond = p->subs[0];
						g.dynamic = refsDynamic(p->subs[0], env);
						p->index = g.id;
						out.push_back(g);
					}
					numberExprGuards(p->subs[0], env, out);
					numberExprGuards(p->subs[1], env, out);
				}
			}

			void numberStatementGuards(AstNode *p, const HybridEnv &env, vector<Guard> &out) {
				for (; p != NULL; p = p->next) {
					if (p->nodeKind == NK_IfEquation) {
						Guard g;
						g.id = (int)out.size();
						g.kind = G_IFEQ;
						g.cond = p->subs[0];
						g.line = p->line;
						g.dynamic = refsDynamic(p->subs[0], env);
						p->index = g.id;
						out.push_back(g);
						numberExprGuards(p->subs[0], env, out);
						if (p->subs[1] != NULL) {
							numberStatementGuards(p->subs[1]->subs[0], env, out);
							numberStatementGuards(p->subs[1]->subs[1], env, out);
						}
					} else if (p->nodeKind == NK_WhenEquation) {
						AstNode *w = p;
						while (w != NULL && w->nodeKind == NK_WhenEquation) {
							Guard g;
							TimeCond tc;
							g.id = (int)out.size();
							g.kind = G_WHEN;
							g.cond = w->subs[0];
							g.line = w->line;
							g.dynamic = !classifyTimeCondition(w->subs[0], tc);
							w->index = g.id;
							out.push_back(g);
							w = w->subs[1] != NULL ? w->subs[1]->subs[1] : NULL;
						}
					} else if (isAssertCall(p)) {
						AssertCall ac;
						string e;
						if (parseAssertCall(p, ac, e)) {
							Guard g;
							g.id = (int)out.size();
							g.kind = G_ASSERT;
							g.cond = ac.cond;
							g.line = p->line;
							g.dynamic = true;
							g.message = ac.message;
							g.warning = ac.warning;
							p->index = g.id;
							out.push_back(g);
						}
					} else {
						numberExprGuards(p, env, out);
					}
				}
			}
		} // namespace

		vector<Guard> collectGuards(SasModel &model, const HybridEnv &env) {
			vector<Guard> guards;
			for (const shared_ptr<AstTree> &tree : model.modelTrees) {
				numberStatementGuards(tree->pHead, env, guards);
			}
			return guards;
		}

		bool evalGuards(const vector<Guard> &guards, const HybridEnv &env, map<int, bool> &truth, string &err) {
			for (const Guard &g : guards) {
				if (!g.dynamic) {
					continue;
				}
				double v;
				string e;
				if (!env.eval(g.cond, v, e)) {
					err = lineMsg(g.line, e);
					return false;
				}
				truth[g.id] = (g.kind == G_ASSERT) ? (v == 0.0) : (v != 0.0);
			}
			return true;
		}

		namespace {
			void gatherRelations(const AstNode *p, std::vector<const AstNode *> &out) {
				if (p == NULL) {
					return;
				}
				if (p->op == OP_GREAT || p->op == OP_LESS || p->op == OP_GREAT_EQ || p->op == OP_LESS_EQ) {
					out.push_back(p);
				}
				gatherRelations(p->subs[0], out);
				gatherRelations(p->subs[1], out);
			}
		} // namespace

		vector<const AstNode *> collectRelations(const vector<Guard> &guards) {
			vector<const AstNode *> rels;
			for (const Guard &g : guards) {
				if (g.dynamic) {
					gatherRelations(g.cond, rels);
				}
			}
			return rels;
		}

		namespace {
			void gatherNames(const AstNode *p, set<string> &out) {
				for (; p != NULL; p = p->next) {
					if (p->op == OP_ID && p->ty != TY_FUNCTION) {
						out.insert(idName(p));
					}
					gatherNames(p->subs[0], out);
					gatherNames(p->subs[1], out);
				}
			}
		} // namespace

		set<string> collectGuardNames(const vector<Guard> &guards) {
			set<string> names;
			for (const Guard &g : guards) {
				if (g.dynamic) {
					gatherNames(g.cond, names);
				}
			}
			return names;
		}

		bool relationMargin(const AstNode *relation, const HybridEnv &env, double &margin, string &err) {
			double a, b;
			if (!env.eval(relation->subs[0], a, err) || !env.eval(relation->subs[1], b, err)) {
				return false;
			}
			margin = a - b;
			return true;
		}

		namespace {

			// Next time strictly after t at which the link's condition becomes true.
			bool linkNextTime(const WhenLink &link, HybridEnv &env, double t, double &te, bool &found, string &err) {
				found = false;
				TimeCond tc;
				if (!classifyTimeCondition(link.cond, tc)) {
					return true; // a state-event condition: detected by the simulator, not scheduled
				}
				if (tc.kind == TC_INITIAL) {
					return true; // fires at the start time only
				}
				env.time = t;
				double a;
				if (!env.eval(tc.a, a, err)) {
					err = lineMsg(link.line, err);
					return false;
				}
				if (tc.kind == TC_AT) {
					if (a > t + tolOf(t)) {
						te = a;
						found = true;
					}
					return true;
				}
				double dt;
				if (!env.eval(tc.b, dt, err)) {
					err = lineMsg(link.line, err);
					return false;
				}
				if (!(dt > 0.0)) {
					err = lineMsg(link.line, "sample() interval must be positive");
					return false;
				}
				double k = 0.0;
				if (t >= a - tolOf(t)) {
					k = floor((t - a) / dt + 1e-9) + 1.0;
				}
				te = a + k * dt;
				found = true;
				return true;
			}

			// Does the link fire at time T? Time-event conditions are tested against T; other conditions
			// fire if their guard became true (is in `fired`).
			bool linkTriggers(const WhenLink &link, HybridEnv &env, double T, bool initialOnly, bool noTimeEvents, const set<int> *fired, bool &fires, string &err) {
				fires = false;
				TimeCond tc;
				if (!classifyTimeCondition(link.cond, tc)) {
					fires = fired != NULL && link.guardId >= 0 && fired->count(link.guardId) > 0;
					return true;
				}
				if (tc.kind == TC_INITIAL) {
					fires = initialOnly && !noTimeEvents;
					return true;
				}
				if (noTimeEvents || (initialOnly && tc.kind != TC_SAMPLE)) {
					return true;
				}
				env.time = T;
				double a;
				if (!env.eval(tc.a, a, err)) {
					err = lineMsg(link.line, err);
					return false;
				}
				if (tc.kind == TC_AT) {
					fires = fabs(a - T) <= tolOf(T);
					return true;
				}
				double dt;
				if (!env.eval(tc.b, dt, err)) {
					err = lineMsg(link.line, err);
					return false;
				}
				if (!(dt > 0.0)) {
					err = lineMsg(link.line, "sample() interval must be positive");
					return false;
				}
				double k = floor((T - a) / dt + 0.5);
				fires = (k >= 0.0) && fabs(a + k * dt - T) <= tolOf(T);
				return true;
			}
		} // namespace

		bool nextEventTime(const vector<WhenClause> &clauses, HybridEnv &env, double t, double &te, string &err) {
			bool any = false;
			double best = 0.0;
			for (const WhenClause &clause : clauses) {
				for (const WhenLink &link : clause.links) {
					double c;
					bool found;
					if (!linkNextTime(link, env, t, c, found, err)) {
						return false;
					}
					if (found && (!any || c < best)) {
						best = c;
						any = true;
					}
				}
			}
			te = best;
			return any;
		}

		bool fireEvents(const vector<WhenClause> &clauses, HybridEnv &env, double T, bool initialOnly, HybridEvent &ev, string &err, const set<int> *fired, bool noTimeEvents) {
			ev = HybridEvent();
			ev.time = T;
			env.time = T;
			env.pre.clear();
			for (const string &n : env.discrete) {
				env.pre[n] = env.values[n];
			}
			map<string, double> assigned;
			for (const WhenClause &clause : clauses) {
				for (const WhenLink &link : clause.links) {
					bool fires;
					if (!linkTriggers(link, env, T, initialOnly, noTimeEvents, fired, fires, err)) {
						return false;
					}
					if (!fires) {
						continue;
					}
					if (ev.line == 0) {
						ev.line = link.line;
					}
					for (const AstNode *s = link.body; s != NULL; s = s->next) {
						if (s->op == OP_EQN && s->subs[0] != NULL && s->subs[0]->op == OP_ID) {
							string n = idName(s->subs[0]);
							double v;
							string e;
							if (!env.eval(s->subs[1], v, e)) {
								err = lineMsg(s->line > 0 ? s->line : link.line, e);
								return false;
							}
							if (!env.discrete.count(n)) {
								err = lineMsg(link.line, "'" + n + "' is not a discrete variable; use reinit() to change a continuous variable");
								return false;
							}
							if (assigned.count(n)) {
								err = lineMsg(link.line, "discrete variable '" + n + "' is assigned by more than one when-branch at time " + to_string(T));
								return false;
							}
							assigned[n] = env.coerce(n, v);
							ev.assignments.push_back(make_pair(n, assigned[n]));
						} else if (isCall(s, "reinit") && s->subs[1] != NULL && s->subs[1]->op == OP_ID && s->subs[1]->next != NULL) {
							string n = idName(s->subs[1]);
							double v;
							string e;
							if (!env.eval(s->subs[1]->next, v, e)) {
								err = lineMsg(link.line, e);
								return false;
							}
							ev.reinits.push_back(make_pair(n, v));
						} else if (isAssertCall(s)) {
							AssertCall ac;
							string e;
							double v;
							if (!parseAssertCall(s, ac, e) || !env.eval(ac.cond, v, e)) {
								err = lineMsg(link.line, e);
								return false;
							}
							if (v == 0.0) {
								ostringstream os;
								os << "assertion violated at t=" << T << ": " << ac.message;
								if (!ac.warning) {
									err = lineMsg(link.line, os.str());
									return false;
								}
								ev.warnings.push_back(os.str());
							}
						} else if (isTerminateCall(s)) {
							string e;
							if (!parseTerminateCall(s, ev.message, e)) {
								err = lineMsg(link.line, e);
								return false;
							}
							ev.terminate = true;
						} else {
							err = lineMsg(link.line, "unsupported statement in when-equation; expect `d = expr;`, `reinit(x, expr);`, `assert(...)` or `terminate(...)`");
							return false;
						}
					}
					break; // only the first triggering branch of a when/elsewhen chain fires
				}
			}
			for (const pair<const string, double> &a : assigned) {
				env.values[a.first] = a.second;
			}
			return true;
		}

		namespace {

			// Specializes the expression chain starting at p (following `next`); takes ownership.
			AstNode *specExpr(AstNode *p, const HybridEnv &env, const map<int, bool> *dyn, string &err);

			AstNode *specOne(AstNode *p, const HybridEnv &env, const map<int, bool> *dyn, string &err) {
				if (p->op == OP_ID && p->ty != TY_FUNCTION) {
					string n = idName(p);
					if (env.discrete.count(n)) {
						map<string, double>::const_iterator it = env.values.find(n);
						AstNode *c = new AstNode();
						c->op = OP_CONST;
						c->ty = TY_REAL;
						c->nodeKind = p->nodeKind;
						c->value.d = it != env.values.end() ? it->second : 0.0;
						delete p;
						return c;
					}
					return p;
				}
				if (p->op == OP_IFEXPR) {
					double c;
					string e;
					AstNode *br = p->subs[1];
					bool known = false;
					if (dyn != NULL && p->index >= 0) {
						map<int, bool>::const_iterator it = dyn->find(p->index);
						if (it != dyn->end()) {
							c = it->second ? 1.0 : 0.0;
							known = true;
						}
					}
					if (br == NULL || (!known && !env.eval(p->subs[0], c, e))) {
						if (err.empty()) {
							err = br == NULL ? "malformed if-expression" : e;
						}
						return p;
					}
					AstNode *chosen;
					if (c != 0.0) {
						chosen = br->subs[0];
						br->subs[0] = NULL;
					} else {
						chosen = br->subs[1];
						br->subs[1] = NULL;
					}
					freeTree(p);
					return specExpr(chosen, env, dyn, err);
				}
				if (p->op == OP_CALL) {
					p->subs[1] = specExpr(p->subs[1], env, dyn, err);
					return p;
				}
				p->subs[0] = specExpr(p->subs[0], env, dyn, err);
				p->subs[1] = specExpr(p->subs[1], env, dyn, err);
				return p;
			}

			AstNode *specExpr(AstNode *p, const HybridEnv &env, const map<int, bool> *dyn, string &err) {
				if (p == NULL) {
					return NULL;
				}
				AstNode *nxt = p->next;
				p->next = NULL;
				AstNode *r = specOne(p, env, dyn, err);
				if (r != NULL) {
					r->next = specExpr(nxt, env, dyn, err);
				} else {
					freeTree(nxt);
				}
				return r;
			}

			void expandStatement(AstNode *st, const HybridEnv &env, const map<int, bool> *dyn, list<shared_ptr<AstTree>> &out, string &err) {
				if (st == NULL) {
					return;
				}
				AstNode *nxt = st->next;
				st->next = NULL;
				if (st->nodeKind == NK_WhenEquation || isAssertCall(st) || isTerminateCall(st)) {
					freeTree(st); // handled by the simulator, not part of the continuous model
				} else if (st->nodeKind == NK_IfEquation) {
					double c;
					string e;
					AstNode *br = st->subs[1];
					AstNode *chosen = NULL;
					bool known = false;
					if (dyn != NULL && st->index >= 0) {
						map<int, bool>::const_iterator it = dyn->find(st->index);
						if (it != dyn->end()) {
							c = it->second ? 1.0 : 0.0;
							known = true;
						}
					}
					if (br == NULL || (!known && !env.eval(st->subs[0], c, e))) {
						if (err.empty()) {
							err = lineMsg(st->line, br == NULL ? "malformed if-equation" : e);
						}
					} else if (c != 0.0) {
						chosen = br->subs[0];
						br->subs[0] = NULL;
					} else {
						chosen = br->subs[1];
						br->subs[1] = NULL;
					}
					freeTree(st);
					expandStatement(chosen, env, dyn, out, err);
				} else {
					st->subs[0] = specExpr(st->subs[0], env, dyn, err);
					st->subs[1] = specExpr(st->subs[1], env, dyn, err);
					shared_ptr<AstTree> tree = make_shared<AstTree>();
					tree->pHead = st;
					tree->line = st->line;
					out.push_back(tree);
				}
				expandStatement(nxt, env, dyn, out, err);
			}
		} // namespace

		bool specializeModel(SasModel &model, const HybridEnv &env, string &err, const map<int, bool> *dynamicTruth) {
			err.clear();
			list<shared_ptr<AstTree>> old;
			old.swap(model.modelTrees);
			for (const shared_ptr<AstTree> &tree : old) {
				AstNode *head = tree->pHead;
				tree->pHead = NULL;
				expandStatement(head, env, dynamicTruth, model.modelTrees, err);
			}
			return err.empty();
		}

	} // namespace core
} // namespace che
