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

#include "sas/MathFunctions.h"
#include "sas/SasExpr.h"
#include "sas/SasLexico.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

// Parser for the supported subset of flat Modelica:
//
//   model <Name> ["description"]
//     { declaration }
//     { "initial" "equation" {statement} | "equation" {statement} }
//   end <Name>;
//
//   declaration := { parameter | constant | discrete } (Real | Integer | Boolean) <id>
//                  [ "(" attribute "=" expr { "," attribute "=" expr } ")" ]
//                  [ "=" expr ]  ["description"] ";"
//   statement   := expr "=" expr ";"                      (simple equation)
//               |  "reinit" "(" id "," expr ")" ";"       (inside when)
//               |  "terminate" | "assert" call ";"
//               |  if-equation | when-equation
//
// The parser only builds the AST (see SasModel.h for the node layouts). Whether the
// solver can run a construct is decided by the checker (ModelChecker.h), so that
// unsupported syntax is reported with a clear message instead of failing in the solver.

namespace che {
	namespace core {

		namespace {

			// ---------------------------------------------------------------------
			// Constant expressions (parameter values, start values, initial equations)
			// ---------------------------------------------------------------------

			bool evalConst(const AstNode *p, SasModel *model, double &out);

			bool evalConstCall(const AstNode *p, SasModel *model, double &out) {
				if (p->subs[0] == NULL || p->subs[0]->op != OP_ID || p->subs[0]->value.p == NULL || p->subs[1] == NULL) {
					return false;
				}
				const string name = (const char *)p->subs[0]->value.p;
				double a = 0.0, b = 0.0;
				if (!evalConst(p->subs[1], model, a)) {
					return false;
				}
				const AstNode *second = p->subs[1]->next;
				if (evalMathFunction(name, a, out)) {
					// sin, cos, exp, sqrt, inv, log, atan, asin, acos, tan, sinh, cosh, tanh, log10
				} else if (name == "abs") {
					out = fabs(a);
				} else if (name == "atan2") {
					if (second == NULL || !evalConst(second, model, b)) {
						return false;
					}
					out = atan2(a, b);
				} else if (name == "min" || name == "max" || name == "pow") {
					if (second == NULL || !evalConst(second, model, b)) {
						return false;
					}
					out = (name == "min") ? fmin(a, b) : (name == "max") ? fmax(a, b)
																		 : pow(a, b);
				} else {
					return false;
				}
				return true;
			}

			// Evaluates numbers, earlier parameters and arithmetic on them.
			bool evalConst(const AstNode *p, SasModel *model, double &out) {
				if (p == NULL) {
					return false;
				}
				double a = 0.0, b = 0.0;
				switch (p->op) {
				case OP_CONST:
					out = (p->ty == TY_INT || p->ty == TY_BOOL) ? (double)p->value.i[0] : p->value.d;
					return true;
				case OP_ID: {
					if (p->ty == TY_FUNCTION || p->value.p == NULL) {
						return false;
					}
					const char *name = (const char *)p->value.p;
					shared_ptr<IdPearl> id = model->lookupIdFromModel(name, (int)strlen(name));
					if (id == nullptr || id->idTy != ID_PARAM || !id->initVal) {
						return false;
					}
					out = (id->ty == TY_INT || id->ty == TY_BOOL) ? (double)id->value.i : id->value.d;
					return true;
				}
				case OP_POS:
					return evalConst(p->subs[0], model, out);
				case OP_MINUS:
					if (!evalConst(p->subs[0], model, a)) {
						return false;
					}
					out = -a;
					return true;
				case OP_ADD:
				case OP_SUB:
				case OP_MUL:
				case OP_DIV:
				case OP_POW:
					if (!evalConst(p->subs[0], model, a) || !evalConst(p->subs[1], model, b)) {
						return false;
					}
					out = (p->op == OP_ADD) ? a + b : (p->op == OP_SUB) ? a - b
												  : (p->op == OP_MUL)	? a * b
												  : (p->op == OP_DIV)	? a / b
																		: pow(a, b);
					return true;
				case OP_CALL:
					return evalConstCall(p, model, out);
				default:
					return false;
				}
			}

			// Stores a constant into an IdPearl according to its declared type.
			void setValue(IdPearl *id, double v) {
				if (id->ty == TY_INT || id->ty == TY_BOOL) {
					id->value.i = (int)v;
				} else {
					id->value.d = v;
				}
				id->initVal = true;
			}

			// ---------------------------------------------------------------------
			// Helpers
			// ---------------------------------------------------------------------

			// Error recovery: skip to just past the next ';', but never past `end`.
			void skipStatement(SasModelParser *parser) {
				while (parser->currentToken != TK_SEMICOLON && parser->currentToken != TK_EOF && parser->currentToken != TK_K_END) {
					getNextToken(parser);
				}
				if (parser->currentToken == TK_SEMICOLON) {
					getNextToken(parser);
				}
			}

			string calleeName(const AstNode *call) {
				if (call != NULL && call->op == OP_CALL && call->subs[0] != NULL && call->subs[0]->op == OP_ID && call->subs[0]->value.p != NULL) {
					return (const char *)call->subs[0]->value.p;
				}
				return "";
			}

			// ---------------------------------------------------------------------
			// Arrays and for-equations
			//
			// An array is a set of ordinary variables named like "x[2,3]". Subscripts must be constant
			// (numbers, parameters, loop variables) and are resolved while parsing, and for-equations are
			// unrolled, so the solver never sees either.
			// ---------------------------------------------------------------------

			const long kMaxArrayElements = 100000;
			const long kMaxLoopIterations = 1000000;

			void freeAst(AstNode *p) {
				while (p != NULL) {
					AstNode *next = p->next;
					freeAst(p->subs[0]);
					freeAst(p->subs[1]);
					delete p;
					p = next;
				}
			}

			struct LoopVar {
				string name;
				long value;
			};

			struct Elab {
				SasModelParser *parser;
				int line;
				vector<LoopVar> env;
				long iterations = 0;
				Elab(SasModelParser *p, int l) : parser(p), line(l) {}
			};

			AstNode *elabChain(Elab &ctx, const AstNode *p);

			AstNode *mkIntConst(long v) {
				AstNode *n = new AstNode();
				n->op = OP_CONST;
				n->ty = TY_INT;
				n->nodeKind = NK_Expression;
				n->value.i[0] = (int)v;
				n->value.i[1] = 0;
				return n;
			}

			string elementName(const string &base, const vector<int> &idx) {
				string s = base + "[";
				for (size_t k = 0; k < idx.size(); k++) {
					s += (k > 0 ? "," : "") + std::to_string(idx[k]);
				}
				return s + "]";
			}

			// All index tuples of an array, last subscript varying fastest.
			vector<vector<int>> allIndices(const vector<int> &dims) {
				vector<vector<int>> out;
				vector<int> idx(dims.size(), 1);
				while (true) {
					out.push_back(idx);
					size_t k = dims.size();
					while (k > 0 && idx[k - 1] == dims[k - 1]) {
						idx[k - 1] = 1;
						k--;
					}
					if (k == 0) {
						break;
					}
					idx[k - 1]++;
				}
				return out;
			}

			bool evalInteger(Elab &ctx, const AstNode *e, long &out) {
				double v = 0.0;
				if (!evalConst(e, ctx.parser->currentModel.get(), v) || v != std::floor(v) || std::fabs(v) > 1e9) {
					return false;
				}
				out = (long)v;
				return true;
			}

			AstNode *elabNode(Elab &ctx, const AstNode *p);

			AstNode *elabIndex(Elab &ctx, const AstNode *p) {
				SasModelParser *parser = ctx.parser;
				const AstNode *base = p->subs[0];
				if (base == NULL || base->op != OP_ID || base->ty == TY_FUNCTION || base->value.p == NULL) {
					parser->error("Line %d: only declared arrays can have subscripts.", ctx.line);
					return mkIntConst(0);
				}
				const string name = (const char *)base->value.p;
				map<string, vector<int>>::const_iterator it = parser->arrayDims.find(name);
				if (it == parser->arrayDims.end()) {
					parser->error("Line %d: '%s' is not an array.", ctx.line, name.c_str());
					return mkIntConst(0);
				}
				const vector<int> &dims = it->second;
				vector<int> idx;
				for (const AstNode *s = p->subs[1]; s != NULL; s = s->next) {
					AstNode *e = elabNode(ctx, s);
					long v = 0;
					const bool ok = evalInteger(ctx, e, v);
					freeAst(e);
					if (!ok) {
						parser->error("Line %d: the subscript of '%s' must be a constant integer expression (numbers, parameters, loop variables).", ctx.line, name.c_str());
						return mkIntConst(0);
					}
					idx.push_back((int)v);
				}
				if (idx.size() != dims.size()) {
					parser->error("Line %d: '%s' has %d dimension(s) and needs one subscript for each.", ctx.line, name.c_str(), (int)dims.size());
					return mkIntConst(0);
				}
				for (size_t k = 0; k < idx.size(); k++) {
					if (idx[k] < 1 || idx[k] > dims[k]) {
						parser->error("Line %d: subscript %d of '%s' is out of range 1..%d.", ctx.line, idx[k], name.c_str(), dims[k]);
						return mkIntConst(0);
					}
				}
				const string full = elementName(name, idx);
				shared_ptr<IdPearl> id = parser->currentModel->lookupIdFromModel(full.c_str(), (int)full.size());
				if (id == nullptr) {
					parser->error("Line %d: '%s' is not declared.", ctx.line, full.c_str());
					return mkIntConst(0);
				}
				AstNode *n = new AstNode();
				n->op = OP_ID;
				n->nodeKind = base->nodeKind;
				n->line = p->line;
				n->value.p = (decltype(n->value.p))id->str;
				return n;
			}

			// size(x, k), sum(x) and prod(x) on a declared array; NULL for other calls.
			AstNode *elabArrayCall(Elab &ctx, const AstNode *p) {
				const string fn = calleeName(p);
				if (fn != "size" && fn != "sum" && fn != "prod") {
					return NULL;
				}
				const AstNode *arg = p->subs[1];
				if (arg == NULL || arg->op != OP_ID || arg->ty == TY_FUNCTION || arg->value.p == NULL) {
					return NULL;
				}
				map<string, vector<int>>::const_iterator it = ctx.parser->arrayDims.find((const char *)arg->value.p);
				if (it == ctx.parser->arrayDims.end()) {
					return NULL;
				}
				const vector<int> &dims = it->second;
				if (fn == "size") {
					long k = 0;
					AstNode *e = arg->next != NULL ? elabNode(ctx, arg->next) : NULL;
					const bool ok = e != NULL && arg->next->next == NULL && evalInteger(ctx, e, k) && k >= 1 && k <= (long)dims.size();
					freeAst(e);
					if (!ok) {
						ctx.parser->error("Line %d: size() needs an array and a constant dimension number, for example size(x, 1).", ctx.line);
						return mkIntConst(1);
					}
					return mkIntConst(dims[k - 1]);
				}
				if (arg->next != NULL) {
					return NULL;
				}
				AstNode *acc = NULL;
				for (const vector<int> &idx : allIndices(dims)) {
					const string full = elementName((const char *)arg->value.p, idx);
					shared_ptr<IdPearl> id = ctx.parser->currentModel->lookupIdFromModel(full.c_str(), (int)full.size());
					if (id == nullptr) {
						continue;
					}
					AstNode *n = new AstNode();
					n->op = OP_ID;
					n->nodeKind = arg->nodeKind;
					n->value.p = (decltype(n->value.p))id->str;
					if (acc == NULL) {
						acc = n;
					} else {
						AstNode *both = new AstNode();
						both->op = fn == "sum" ? OP_ADD : OP_MUL;
						both->nodeKind = NK_Expression;
						both->subs[0] = acc;
						both->subs[1] = n;
						acc = both;
					}
				}
				return acc != NULL ? acc : mkIntConst(fn == "sum" ? 0 : 1);
			}

			// Unrolls a for-equation: the body once per value of the iterator.
			AstNode *expandFor(Elab &ctx, const AstNode *p) {
				SasModelParser *parser = ctx.parser;
				const string iter = p->value.p != NULL ? (const char *)p->value.p : "";
				vector<long> bounds;
				for (const AstNode *r = p->subs[0]; r != NULL; r = r->next) {
					AstNode *e = elabNode(ctx, r);
					long v = 0;
					const bool ok = evalInteger(ctx, e, v);
					freeAst(e);
					if (!ok) {
						parser->error("Line %d: the range of the for-loop over '%s' must be constant integer expressions.", ctx.line, iter.c_str());
						return NULL;
					}
					bounds.push_back(v);
				}
				if (bounds.size() < 2 || bounds.size() > 3) {
					parser->error("Line %d: the range of the for-loop over '%s' must be `first:last` or `first:step:last`.", ctx.line, iter.c_str());
					return NULL;
				}
				const long lo = bounds[0];
				const long step = bounds.size() == 3 ? bounds[1] : 1;
				const long hi = bounds.back();
				if (step == 0) {
					parser->error("Line %d: the step of the for-loop over '%s' is zero.", ctx.line, iter.c_str());
					return NULL;
				}
				long count = 0;
				if (step > 0 && hi >= lo) {
					count = (hi - lo) / step + 1;
				} else if (step < 0 && lo >= hi) {
					count = (lo - hi) / (-step) + 1;
				}
				if (ctx.iterations + count > kMaxLoopIterations) {
					parser->error("Line %d: the for-loops unroll to more than %ld iterations.", ctx.line, kMaxLoopIterations);
					return NULL;
				}
				ctx.iterations += count;
				AstNode *head = NULL;
				AstNode **tail = &head;
				for (long k = 0; k < count; k++) {
					ctx.env.push_back(LoopVar{iter, lo + k * step});
					AstNode *body = elabChain(ctx, p->subs[1]);
					ctx.env.pop_back();
					*tail = body;
					while (*tail != NULL) {
						tail = &(*tail)->next;
					}
				}
				return head;
			}

			// Copy of one node (without `next`) with loop variables replaced by their values and subscripts resolved.
			AstNode *elabNode(Elab &ctx, const AstNode *p) {
				if (p->op == OP_INDEX) {
					return elabIndex(ctx, p);
				}
				if (p->op == OP_ID && p->ty != TY_FUNCTION && p->value.p != NULL) {
					const string name = (const char *)p->value.p;
					for (size_t k = ctx.env.size(); k > 0; k--) {
						if (ctx.env[k - 1].name == name) {
							return mkIntConst(ctx.env[k - 1].value);
						}
					}
					if (ctx.parser->arrayDims.count(name) > 0) {
						ctx.parser->error("Line %d: array '%s' is used without subscripts; whole-array expressions are not supported, use a for-equation.", ctx.line, name.c_str());
					}
				}
				if (p->op == OP_CALL) {
					AstNode *r = elabArrayCall(ctx, p);
					if (r != NULL) {
						return r;
					}
				}
				AstNode *c = new AstNode(*p);
				c->next = NULL;
				c->index = -1;
				c->subs[0] = elabChain(ctx, p->subs[0]);
				c->subs[1] = elabChain(ctx, p->subs[1]);
				return c;
			}

			// Elaborates a chain of expressions or statements; a for-equation becomes the statements it unrolls to.
			AstNode *elabChain(Elab &ctx, const AstNode *p) {
				AstNode *head = NULL;
				AstNode **tail = &head;
				for (; p != NULL; p = p->next) {
					AstNode *r = p->nodeKind == NK_ForEquation ? expandFor(ctx, p) : elabNode(ctx, p);
					*tail = r;
					while (*tail != NULL) {
						tail = &(*tail)->next;
					}
				}
				return head;
			}

			// A value in a declaration (parameter value, start value, binding) with subscripts resolved.
			AstNode *parseDeclValue(SasModelParser *parser) {
				Elab ctx(parser, parser->coord.line);
				AstNode *raw = parseValueExpr(parser);
				AstNode *value = elabChain(ctx, raw);
				freeAst(raw);
				return value;
			}

			// The elements of an array literal, nested literals flattened (a scalar gives itself).
			void flattenValues(const AstNode *v, vector<const AstNode *> &out) {
				if (v != NULL && v->op == OP_ARRAY) {
					for (const AstNode *e = v->subs[0]; e != NULL; e = e->next) {
						flattenValues(e, out);
					}
				} else {
					out.push_back(v);
				}
			}

			// One value per element for an attribute: an array literal of the right size, or one value for all.
			bool valuesPerElement(SasModelParser *parser, const AstNode *value, size_t n, const string &what, vector<const AstNode *> &per) {
				vector<const AstNode *> vals;
				flattenValues(value, vals);
				if (vals.size() == n) {
					per = vals;
				} else if (vals.size() == 1) {
					per.assign(n, vals[0]);
				} else {
					parser->error("%s has %d value(s) but the array has %d element(s).", what.c_str(), (int)vals.size(), (int)n);
					return false;
				}
				return true;
			}

			// ---------------------------------------------------------------------
			// Declarations
			// ---------------------------------------------------------------------

			// Parses  ( [each] name = value {, [each] name = value} )  after a variable name.
			// `start` becomes the initial value of every element (an array literal gives one value per
			// element); `unit` is stored; other attributes are parsed and ignored with a warning.
			void parseModification(SasModelParser *parser, const vector<shared_ptr<IdPearl>> &elems, const string &name) {
				// Attribute names (start, fixed, ...) are not model symbols: scan them
				// without entering them into the id table.
				parser->noInject = true;
				getNextToken(parser); // consume '('
				parser->noInject = false;
				while (true) {
					string attr;
					if (parser->currentToken == TK_K_EACH) {
						parser->noInject = true;
						getNextToken(parser);
						parser->noInject = false;
					}
					if (parser->currentToken == TK_K_UNIT) {
						attr = "unit";
					} else if (parser->currentToken == TK_ID && parser->value.p != NULL) {
						attr = (const char *)parser->value.p;
					} else {
						parser->error("Expect attribute name in modification.");
						break;
					}
					getNextToken(parser);
					expect(parser, TK_EQN);
					if (attr == "unit") {
						if (parser->currentToken == TK_STRING) {
							for (const shared_ptr<IdPearl> &e : elems) {
								e->setUnit((uchar *)parser->value.p);
							}
							getNextToken(parser);
						} else {
							parser->error("Expect a string for the unit.");
						}
					} else {
						AstNode *value = parseDeclValue(parser);
						if (attr == "start") {
							vector<const AstNode *> per;
							if (valuesPerElement(parser, value, elems.size(), "Attribute start of '" + name + "'", per)) {
								for (size_t k = 0; k < elems.size(); k++) {
									double v = 0.0;
									if (evalConst(per[k], parser->currentModel.get(), v)) {
										setValue(elems[k].get(), v);
									} else {
										parser->error("Attribute start of '%s' must be a constant expression.", name.c_str());
										break;
									}
								}
							}
						} else if (attr != "fixed") {
							parser->warn("Attribute '%s' of '%s' is ignored.", attr.c_str(), name.c_str());
						}
						freeAst(value);
					}
					if (parser->currentToken != TK_COMMA) {
						break;
					}
					parser->noInject = true;
					getNextToken(parser); // consume ',' and scan the next attribute name
					parser->noInject = false;
				}
				expect(parser, TK_RPAREN);
			}

			bool startsDeclaration(int tok) {
				return tok == TK_K_PARA || tok == TK_K_CONSTANT || tok == TK_K_DISCRETE || tok == TK_K_INPUT || tok == TK_K_OUTPUT || tok == TK_K_REAL || tok == TK_K_INTEGER || tok == TK_K_BOOLEAN;
			}

			// Skips `annotation(...)`, which does not change the equations. The tokens inside are scanned without
			// entering their names into the id table.
			void skipAnnotation(SasModelParser *parser) {
				parser->noInject = true;
				getNextToken(parser); // consume `annotation`
				int depth = 0;
				while (parser->currentToken != TK_EOF) {
					if (parser->currentToken == TK_LPAREN) {
						depth++;
					} else if (parser->currentToken == TK_RPAREN) {
						depth--;
						if (depth <= 0) {
							parser->noInject = false;
							getNextToken(parser);
							return;
						}
					} else if (depth == 0) {
						break; // `annotation` without parentheses
					}
					getNextToken(parser);
				}
				parser->noInject = false;
			}

			// Array dimensions  [n, m]  after a type or a name.
			void parseDims(SasModelParser *parser, vector<int> &dims) {
				getNextToken(parser); // consume '['
				while (true) {
					AstNode *e = parseDeclValue(parser);
					double v = 0.0;
					if (!evalConst(e, parser->currentModel.get(), v) || v != std::floor(v) || v < 1 || v > kMaxArrayElements) {
						parser->error("Line %d: array sizes must be constant positive integers (numbers or earlier parameters).", parser->coord.line);
						v = 1;
					}
					freeAst(e);
					dims.push_back((int)v);
					if (parser->currentToken != TK_COMMA) {
						break;
					}
					getNextToken(parser);
				}
				expect(parser, TK_RBRACKET);
			}

			// One variable (or array) of a declaration: name, dimensions, modification, binding and description.
			void parseDeclarator(SasModelParser *parser, Type ty, bool isParam, bool isDiscrete, bool isInput, const vector<int> &typeDims) {
				shared_ptr<IdPearl> var = parser->currentId;
				const string name = (const char *)var->str;
				const int line = parser->coord.line;
				getNextToken(parser);
				vector<int> dims = typeDims;
				if (parser->currentToken == TK_LBRACKET) {
					parseDims(parser, dims);
				}
				const IdType idTy = isParam ? ID_PARAM : ((isDiscrete || ty != TY_REAL) ? ID_DISC : ID_VAR);

				// The variable itself, or one variable per array element.
				vector<shared_ptr<IdPearl>> elems;
				if (var->idTy != ID_NONE || parser->arrayDims.count(name) > 0) {
					parser->error("Redeclaration of '%s'.", name.c_str());
				}
				if (dims.empty()) {
					elems.push_back(var);
				} else {
					long total = 1;
					for (int d : dims) {
						total *= d;
					}
					if (total > kMaxArrayElements) {
						parser->error("Line %d: array '%s' has more than %ld elements.", line, name.c_str(), kMaxArrayElements);
						dims.assign(1, 1);
					}
					parser->arrayDims[name] = dims;
					for (const vector<int> &idx : allIndices(dims)) {
						const string full = elementName(name, idx);
						shared_ptr<IdPearl> el = parser->currentModel->lookupIdentifier(full.c_str(), (int)full.size());
						elems.push_back(el);
					}
				}
				for (const shared_ptr<IdPearl> &el : elems) {
					el->ty = ty;
					el->line = line;
					el->idTy = idTy;
				}

				if (parser->currentToken == TK_LPAREN) {
					parseModification(parser, elems, name);
				}

				if (parser->currentToken == TK_EQN) {
					getNextToken(parser);
					AstNode *value = parseDeclValue(parser);
					vector<const AstNode *> vals;
					flattenValues(value, vals);
					if (vals.size() != elems.size()) {
						parser->error("Line %d: the binding of '%s' has %d value(s) but it has %d element(s).", line, name.c_str(), (int)vals.size(), (int)elems.size());
					} else if (!isParam && idTy != ID_VAR) {
						parser->error("Binding equations on discrete variable '%s' are not supported; assign it in a when-equation.", name.c_str());
					} else {
						for (size_t k = 0; k < elems.size(); k++) {
							double v = 0.0;
							if (isParam) {
								if (evalConst(vals[k], parser->currentModel.get(), v)) {
									setValue(elems[k].get(), v);
								} else {
									parser->error("Value of parameter '%s' must be a constant expression.", name.c_str());
									break;
								}
							} else {
								// Real y = expr; is the equation y = expr.
								AstNode *lhs = new AstNode();
								lhs->op = OP_ID;
								lhs->value.p = (decltype(lhs->value.p))elems[k]->str;
								Elab ctx(parser, line);
								AstNode *eqn = new AstNode();
								eqn->op = OP_EQN;
								eqn->subs[0] = lhs;
								eqn->subs[1] = elabNode(ctx, vals[k]);
								shared_ptr<AstTree> tree = make_shared<AstTree>();
								tree->pHead = eqn;
								tree->line = line;
								parser->currentModel->modelTrees.push_back(tree);
							}
						}
					}
					freeAst(value);
				} else if (isParam) {
					parser->error("Parameter '%s' has no value.", name.c_str());
				} else if (isInput) {
					parser->error("Input '%s' has no value; give it a binding, for example `input Real u = 1;`.", name.c_str());
				}

				if (parser->currentToken == TK_STRING) {
					for (const shared_ptr<IdPearl> &el : elems) {
						el->setAnnotation((uchar *)parser->value.p);
					}
					getNextToken(parser);
				}
				if (parser->currentToken == TK_K_ANNOTATION) {
					skipAnnotation(parser);
				}
			}

			// One declaration statement: prefixes, a type and one or more declarators separated by commas.
			void parseDeclaration(SasModelParser *parser) {
				bool isParam = false, isDiscrete = false, isInput = false;
				while (parser->currentToken == TK_K_PARA || parser->currentToken == TK_K_CONSTANT || parser->currentToken == TK_K_DISCRETE ||
					   parser->currentToken == TK_K_INPUT || parser->currentToken == TK_K_OUTPUT) {
					if (parser->currentToken == TK_K_DISCRETE) {
						isDiscrete = true;
					} else if (parser->currentToken == TK_K_INPUT) {
						isInput = true;
					} else if (parser->currentToken != TK_K_OUTPUT) { // an output is an ordinary variable
						isParam = true;
					}
					getNextToken(parser);
				}
				Type ty = TY_REAL;
				if (parser->currentToken == TK_K_REAL) {
					ty = TY_REAL;
				} else if (parser->currentToken == TK_K_INTEGER) {
					ty = TY_INT;
				} else if (parser->currentToken == TK_K_BOOLEAN) {
					ty = TY_BOOL;
				} else {
					parser->error("Expect a type (Real, Integer or Boolean).");
					skipStatement(parser);
					return;
				}
				getNextToken(parser);
				vector<int> typeDims; // Real[3] x;
				if (parser->currentToken == TK_LBRACKET) {
					parseDims(parser, typeDims);
				}
				while (true) {
					if (parser->currentToken != TK_ID) {
						parser->error("actual token %s: expect identifier.", Lexer::getInstance().tokenStrings[parser->currentToken]);
						skipStatement(parser);
						return;
					}
					parseDeclarator(parser, ty, isParam, isDiscrete, isInput, typeDims);
					if (parser->currentToken != TK_COMMA) {
						break;
					}
					getNextToken(parser);
				}
				expect(parser, TK_SEMICOLON);
			}

			// ---------------------------------------------------------------------
			// Equations
			// ---------------------------------------------------------------------

			enum StmtContext {
				CTX_TOP,	 // equation section
				CTX_INITIAL, // initial equation section
				CTX_IF,		 // inside an if-equation (outside any when)
				CTX_WHEN	 // inside a when-equation
			};

			AstNode *parseStatement(SasModelParser *parser, StmtContext ctx);

			bool endsStatementList(int tok) {
				return tok == TK_EOF || tok == TK_K_END || tok == TK_K_ELSE || tok == TK_K_ELSEIF || tok == TK_K_ELSEWHEN || tok == TK_K_EQN || tok == TK_K_INIT || tok == TK_K_ANNOTATION;
			}

			// Statements are linked through AstNode::next.
			AstNode *parseStatementList(SasModelParser *parser, StmtContext ctx) {
				AstNode *head = NULL, *tail = NULL;
				while (!endsStatementList(parser->currentToken)) {
					uchar *before = parser->reader.cursur;
					AstNode *st = parseStatement(parser, ctx);
					if (st != NULL) {
						if (head == NULL) {
							head = st;
						} else {
							tail->next = st;
						}
						tail = st;
					}
					if (parser->reader.cursur == before) { // no progress: avoid looping forever
						getNextToken(parser);
					}
				}
				return head;
			}

			void expectEnd(SasModelParser *parser, int what) {
				expect(parser, TK_K_END);
				expect(parser, what);
				expect(parser, TK_SEMICOLON);
			}

			// Handles `x = <const>` / `der(x) = <const>` in an initial equation section.
			void applyInitialEquation(SasModelParser *parser, AstNode *eq, int line) {
				SasModel *model = parser->currentModel.get();
				AstNode *lhs = eq->subs[0];
				bool isDer = false;
				AstNode *target = lhs;
				if (lhs != NULL && lhs->op == OP_CALL && calleeName(lhs) == "der") {
					isDer = true;
					target = lhs->subs[1];
				}
				double v = 0.0;
				shared_ptr<IdPearl> id = nullptr;
				if (target != NULL && target->op == OP_ID && target->ty != TY_FUNCTION && target->value.p != NULL) {
					const char *name = (const char *)target->value.p;
					id = model->lookupIdFromModel(name, (int)strlen(name));
				}
				if (id == nullptr || (id->idTy != ID_VAR && id->idTy != ID_DISC)) {
					parser->error("Line %d: unsupported initial equation; expect `x = <constant>` or `der(x) = <constant>` for a declared variable.", line);
				} else if (!evalConst(eq->subs[1], model, v)) {
					parser->error("Line %d: initial value of '%s' must be a constant expression.", line, (const char *)id->str);
				} else if (isDer) {
					id->derVal.d = v;
					id->initDer = true;
				} else {
					setValue(id.get(), v);
				}
			}

			AstNode *parseIfEquation(SasModelParser *parser, StmtContext ctx, bool isElseIf) {
				AstNode *node = new AstNode();
				node->nodeKind = NK_IfEquation;
				node->line = parser->coord.line; // line of the `if` / `elseif`
				getNextToken(parser);			 // consume `if` / `elseif`
				node->subs[0] = parseValueExpr(parser);
				expect(parser, TK_K_THEN);
				AstNode *branch = new AstNode();
				branch->op = OP_BRANCH;
				branch->nodeKind = NK_Equation;
				branch->subs[0] = parseStatementList(parser, ctx);
				if (parser->currentToken == TK_K_ELSEIF) {
					// `elseif` is desugared into a nested if in the else branch.
					branch->subs[1] = parseIfEquation(parser, ctx, true);
				} else if (parser->currentToken == TK_K_ELSE) {
					getNextToken(parser);
					branch->subs[1] = parseStatementList(parser, ctx);
				}
				node->subs[1] = branch;
				if (!isElseIf) {
					expectEnd(parser, TK_K_IF);
				}
				return node;
			}

			AstNode *parseWhenEquation(SasModelParser *parser, StmtContext ctx, bool isElseWhen) {
				if (ctx == CTX_WHEN) {
					parser->error("when-equations cannot be nested.");
				}
				AstNode *node = new AstNode();
				node->nodeKind = NK_WhenEquation;
				node->line = parser->coord.line; // line of the `when` / `elsewhen`
				getNextToken(parser);			 // consume `when` / `elsewhen`
				node->subs[0] = parseValueExpr(parser);
				expect(parser, TK_K_THEN);
				AstNode *branch = new AstNode();
				branch->op = OP_BRANCH;
				branch->nodeKind = NK_Equation;
				branch->subs[0] = parseStatementList(parser, CTX_WHEN);
				if (parser->currentToken == TK_K_ELSEWHEN) {
					branch->subs[1] = parseWhenEquation(parser, ctx, true);
				}
				node->subs[1] = branch;
				if (!isElseWhen) {
					expectEnd(parser, TK_K_WHEN);
				}
				return node;
			}

			// for i in 1:n loop ... end for;  and  for i in 1:n, j in 1:m loop ... end for;
			// The result is a NK_ForEquation node (nested for several iterators); elaboration unrolls it.
			AstNode *parseForEquation(SasModelParser *parser, StmtContext ctx) {
				const int line = parser->coord.line;
				getNextToken(parser); // consume `for`
				vector<pair<uchar *, AstNode *>> iterators;
				while (true) {
					if (parser->currentToken != TK_ID) {
						parser->error("Line %d: expect a loop variable after `for`.", line);
						skipStatement(parser);
						return NULL;
					}
					uchar *name = parser->currentId->str;
					getNextToken(parser);
					expect(parser, TK_K_IN);
					AstNode *first = parseValueExpr(parser);
					AstNode *last = first;
					int parts = 1;
					while (parser->currentToken == TK_COLON && parts < 3) {
						getNextToken(parser);
						AstNode *e = parseValueExpr(parser);
						last->next = e;
						last = e;
						parts++;
					}
					if (parts < 2) {
						parser->error("Line %d: the range of a for-loop must be `first:last` or `first:step:last`.", line);
					}
					iterators.push_back(make_pair(name, first));
					if (parser->currentToken != TK_COMMA) {
						break;
					}
					getNextToken(parser);
				}
				expect(parser, TK_K_LOOP);
				AstNode *node = parseStatementList(parser, ctx);
				expectEnd(parser, TK_K_FOR);
				for (size_t k = iterators.size(); k > 0; k--) {
					AstNode *loop = new AstNode();
					loop->nodeKind = NK_ForEquation;
					loop->line = line;
					loop->value.p = (decltype(loop->value.p))iterators[k - 1].first;
					loop->subs[0] = iterators[k - 1].second;
					loop->subs[1] = node;
					node = loop;
				}
				return node;
			}

			// Returns the statement node, or NULL when nothing is to be kept (error recovery).
			AstNode *parseStatement(SasModelParser *parser, StmtContext ctx) {
				const int line = parser->coord.line;
				if (parser->currentToken == TK_K_FOR) {
					return parseForEquation(parser, ctx);
				}
				if (parser->currentToken == TK_K_IF) {
					if (ctx == CTX_INITIAL) {
						parser->error("Line %d: if-equations are not supported in initial equation sections.", line);
					}
					AstNode *node = parseIfEquation(parser, ctx == CTX_TOP ? CTX_IF : ctx, false);
					return node;
				}
				if (parser->currentToken == TK_K_WHEN) {
					if (ctx == CTX_INITIAL) {
						parser->error("Line %d: when-equations are not supported in initial equation sections.", line);
					}
					return parseWhenEquation(parser, ctx, false);
				}
				if (parser->currentToken == TK_SEMICOLON) {
					parser->error("Line %d: empty statement.", line);
					getNextToken(parser);
					return NULL;
				}
				AstNode *expr = parseExpr(parser);
				bool valid = false;
				if (expr != NULL && expr->op == OP_EQN) {
					valid = true;
				} else if (expr != NULL && expr->op == OP_CALL) {
					const string name = calleeName(expr);
					if (name == "reinit") {
						valid = true;
						if (ctx != CTX_WHEN) {
							parser->error("Line %d: reinit() is only allowed inside a when-equation.", line);
						}
					} else if (name == "terminate" || name == "assert") {
						valid = true;
					}
				}
				if (!valid) {
					parser->error("Line %d: expect an equation `lhs = rhs;`.", line);
					skipStatement(parser);
					return NULL;
				}
				expect(parser, TK_SEMICOLON);
				if (ctx == CTX_INITIAL && expr->op != OP_EQN) {
					parser->error("Line %d: unsupported statement in initial equation section.", line);
					return NULL;
				}
				return expr; // initial equations are applied by parseEquationSection, after subscripts are resolved
			}

			void parseEquationSection(SasModelParser *parser, StmtContext ctx) {
				while (!endsStatementList(parser->currentToken)) {
					uchar *before = parser->reader.cursur;
					const int line = parser->coord.line;
					AstNode *st = parseStatement(parser, ctx);
					if (st != NULL) {
						// Subscripts are resolved and for-equations unrolled: one statement may become many.
						Elab elab(parser, line);
						AstNode *chain = elabChain(elab, st);
						freeAst(st);
						while (chain != NULL) {
							AstNode *s = chain;
							chain = chain->next;
							s->next = NULL;
							if (ctx == CTX_INITIAL) {
								applyInitialEquation(parser, s, line);
								freeAst(s);
							} else {
								shared_ptr<AstTree> tree = make_shared<AstTree>();
								tree->pHead = s;
								tree->line = line;
								parser->currentModel->modelTrees.push_back(tree);
							}
						}
					}
					if (parser->reader.cursur == before) {
						getNextToken(parser);
					}
				}
			}

		} // namespace

		void parseSingleSasModel(SasModelParser *parser) {
			getNextToken(parser);
			if (parser->currentToken != TK_K_MDL) {
				parser->error("actual token %s: expect \"model\".", Lexer::getInstance().tokenStrings[parser->currentToken]);
			}
			getNextToken(parser);
			if (parser->currentToken != TK_ID) {
				parser->error("actual token %s: expect identifier.", Lexer::getInstance().tokenStrings[parser->currentToken]);
				parser->currentId = parser->currentModel->lookupIdentifier("Unknown", 7);
			}
			parser->currentModel->modelInfo = parser->currentId;
			parser->currentId->idTy = ID_MODEL;

			getNextToken(parser);
			if (parser->currentToken == TK_STRING) {
				parser->currentModel->modelInfo->setAnnotation((uchar *)parser->value.p);
				getNextToken(parser);
			}

			while (startsDeclaration(parser->currentToken) || parser->currentToken == TK_K_ANNOTATION) {
				if (parser->currentToken == TK_K_ANNOTATION) {
					skipAnnotation(parser);
					expect(parser, TK_SEMICOLON);
				} else {
					parseDeclaration(parser);
				}
			}

			while (parser->currentToken != TK_K_END && parser->currentToken != TK_EOF) {
				if (parser->currentToken == TK_K_ANNOTATION) {
					skipAnnotation(parser);
					expect(parser, TK_SEMICOLON);
				} else if (parser->currentToken == TK_K_INIT) {
					getNextToken(parser);
					expect(parser, TK_K_EQN);
					parseEquationSection(parser, CTX_INITIAL);
				} else if (parser->currentToken == TK_K_EQN) {
					getNextToken(parser);
					parseEquationSection(parser, CTX_TOP);
				} else {
					if (startsDeclaration(parser->currentToken)) {
						parser->error("Declarations must come before the equation sections.");
					} else if (parser->currentToken == TK_K_EXTENDS) {
						parser->error("Inheritance (extends) is not supported; flatten the model first.");
					} else if (parser->currentToken == TK_K_ALGORITHM) {
						parser->error("Algorithm sections are not supported; write the statements as equations.");
					} else if (parser->currentToken == TK_K_FUNCTION || parser->currentToken == TK_K_PACKAGE || parser->currentToken == TK_K_CLASS) {
						parser->error("Nested classes and user-defined functions are not supported.");
					} else {
						parser->error("actual token %s: expect `initial equation`, `equation` or `end`.", Lexer::getInstance().tokenStrings[parser->currentToken]);
					}
					skipStatement(parser);
				}
			}

			if (parser->currentToken == TK_K_END) {
				getNextToken(parser);
				if (parser->currentToken == TK_ID && parser->currentId != parser->currentModel->modelInfo) {
					parser->warn("End name does not match the model name.");
				}
				getNextToken(parser);
				expect(parser, TK_SEMICOLON);
			} else {
				parser->error("Missing end statement of model.");
			}
		}
	} // namespace core
} // namespace che
