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
#include <string>

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
			// Declarations
			// ---------------------------------------------------------------------

			// Parses  ( name = value {, name = value} )  after a variable name.
			// `start` becomes the initial value; `unit` is stored; other attributes are
			// parsed and ignored with a warning.
			void parseModification(SasModelParser *parser, const shared_ptr<IdPearl> &var) {
				// Attribute names (start, fixed, ...) are not model symbols: scan them
				// without entering them into the id table.
				parser->noInject = true;
				getNextToken(parser); // consume '('
				parser->noInject = false;
				while (true) {
					string attr;
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
							var->setUnit((uchar *)parser->value.p);
							getNextToken(parser);
						} else {
							parser->error("Expect a string for the unit.");
						}
					} else {
						AstNode *value = parseValueExpr(parser);
						if (attr == "start") {
							double v = 0.0;
							if (evalConst(value, parser->currentModel.get(), v)) {
								setValue(var.get(), v);
							} else {
								parser->error("Attribute start of '%s' must be a constant expression.", (const char *)var->str);
							}
						} else if (attr != "fixed") {
							parser->warn("Attribute '%s' of '%s' is ignored.", attr.c_str(), (const char *)var->str);
						}
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

			// One variable of a declaration: name, modification, binding and description string.
			void parseDeclarator(SasModelParser *parser, Type ty, bool isParam, bool isDiscrete, bool isInput) {
				shared_ptr<IdPearl> var = parser->currentId;
				const int line = parser->coord.line;
				if (var->idTy != ID_NONE) {
					parser->error("Redeclaration of '%s'.", (const char *)var->str);
				}
				var->ty = ty;
				var->line = line;
				var->idTy = isParam ? ID_PARAM : ((isDiscrete || ty != TY_REAL) ? ID_DISC : ID_VAR);
				getNextToken(parser);

				if (parser->currentToken == TK_LPAREN) {
					parseModification(parser, var);
				}

				if (parser->currentToken == TK_EQN) {
					getNextToken(parser);
					AstNode *value = parseValueExpr(parser);
					double v = 0.0;
					if (!isParam) {
						if (var->idTy == ID_VAR) {
							// Real y = expr; is the equation y = expr.
							AstNode *lhs = new AstNode();
							lhs->op = OP_ID;
							lhs->value.p = (decltype(lhs->value.p))var->str;
							AstNode *eqn = new AstNode();
							eqn->op = OP_EQN;
							eqn->subs[0] = lhs;
							eqn->subs[1] = value;
							shared_ptr<AstTree> tree = make_shared<AstTree>();
							tree->pHead = eqn;
							tree->line = line;
							parser->currentModel->modelTrees.push_back(tree);
						} else {
							parser->error("Binding equations on discrete variable '%s' are not supported; assign it in a when-equation.", (const char *)var->str);
						}
					} else if (evalConst(value, parser->currentModel.get(), v)) {
						setValue(var.get(), v);
					} else {
						parser->error("Value of parameter '%s' must be a constant expression.", (const char *)var->str);
					}
				} else if (isParam) {
					parser->error("Parameter '%s' has no value.", (const char *)var->str);
				} else if (isInput) {
					parser->error("Input '%s' has no value; give it a binding, for example `input Real u = 1;`.", (const char *)var->str);
				}

				if (parser->currentToken == TK_STRING) {
					var->setAnnotation((uchar *)parser->value.p);
					getNextToken(parser);
				}
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
				while (true) {
					if (parser->currentToken != TK_ID) {
						parser->error("actual token %s: expect identifier.", Lexer::getInstance().tokenStrings[parser->currentToken]);
						skipStatement(parser);
						return;
					}
					parseDeclarator(parser, ty, isParam, isDiscrete, isInput);
					if (parser->currentToken != TK_COMMA) {
						break;
					}
					getNextToken(parser);
				}
				if (parser->currentToken == TK_K_ANNOTATION) {
					skipAnnotation(parser);
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

			// Returns the statement node, or NULL when nothing is to be kept (error
			// recovery, or an initial equation that was applied directly).
			AstNode *parseStatement(SasModelParser *parser, StmtContext ctx) {
				const int line = parser->coord.line;
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
				if (ctx == CTX_INITIAL) {
					if (expr->op == OP_EQN) {
						applyInitialEquation(parser, expr, line);
					} else {
						parser->error("Line %d: unsupported statement in initial equation section.", line);
					}
					return NULL;
				}
				return expr;
			}

			void parseEquationSection(SasModelParser *parser, StmtContext ctx) {
				while (!endsStatementList(parser->currentToken)) {
					uchar *before = parser->reader.cursur;
					const int line = parser->coord.line;
					AstNode *st = parseStatement(parser, ctx);
					if (st != NULL) {
						shared_ptr<AstTree> tree = make_shared<AstTree>();
						tree->pHead = st;
						tree->line = line;
						parser->currentModel->modelTrees.push_back(tree);
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
