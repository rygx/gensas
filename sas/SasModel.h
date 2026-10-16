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
#ifndef SAS_MODEL_H
#define SAS_MODEL_H

#include "sas/SasConfig.h"
#include <list>
#include <string>
#include <vector>
#include <iostream>
#include <memory>

namespace che {
	namespace core {
		typedef union {
			void *p;
			float f;
			double d;
			int i[2];
		} TokenValue;

		typedef union {
			float f;
			double d;
			int i;
		} NumValue;

		typedef TokenValue IdValue;

		typedef struct {
			char *filename = NULL;
			int line = 0;
			int col = 0;
		} Coord;

		// A message about the model text, for editors: where it is (1-based line and column; 0 if unknown)
		// and what is wrong. `source` says which stage found it: "parser" or "checker".
		struct Diagnostic {
			enum Severity { Error,
							Warning };
			Severity severity = Error;
			int line = 0;
			int column = 0;
			std::string message;
			std::string source;
		};

		typedef enum {
#define TOKEN(k, n) k,
#include "sas/token.txt"
#undef TOKEN
			NUM_TK,
		} TokenKind;

		typedef enum {
			NK_NA,
			NK_Expression,
			NK_Equation,
			NK_Decl,
			NK_Condition,
			NK_Domain,
			NK_Token,
			NK_Assignment,
			NK_ModelBlock,
			NK_ModelDef,
			NK_EndModel,
			NK_VarStatement,
			NK_ParameterStatement,
			NK_InitialBlock,
			NK_InitialCondition,
			NK_EquationBlock,
			NK_IfEquation,	 // subs[0]=condition; subs[1]=OP_BRANCH{ subs[0]=then-statements, subs[1]=else-statements }.
							 // Statement lists are linked through `next`; `elseif` is desugared to a nested if in the else list.
			NK_WhenEquation, // subs[0]=condition; subs[1]=OP_BRANCH{ subs[0]=body-statements, subs[1]=next `elsewhen` NK_WhenEquation or NULL }
			NK_ForEquation,	 // value.p=iterator name; subs[0]=range (lo, hi) or (lo, step, hi) linked through `next`; subs[1]=body statements.
							 // Only exists while parsing: for-equations are unrolled before the parser returns.
		} NodeKind;

		typedef enum {
			TY_NULL,
			TY_REAL,
			TY_INT,
			TY_STR,
			TY_FUNCTION,
			TY_BOOL,
		} Type;

		typedef enum {
			ID_NONE,
			ID_PARAM,
			ID_VAR,
			ID_IVAR,
			ID_TMPVAR,
			ID_MODEL,
			ID_FUNC,
			ID_DISC // discrete (Integer/Boolean/discrete Real) variable
		} IdType;

		enum Op {
#define OPINFO(op, prec, name, func, opcode) op,
#include "sas/opinfo.txt"
#undef OPINFO
		};

		typedef struct keyword {
			const char *name;
			int len;
			int tok;
		} func;

		class Token {
		public:
			TokenKind kind;
			TokenValue value;
		};

		class AstNode {
		public:
			AstNode *next = NULL;
			TokenValue value = {};
			NodeKind nodeKind = NK_NA;
			Op op = OP_NONE;
			Type ty = TY_NULL;
			int isfunc = 0;
			int index = -1;
			int pairIdx = -1;
			int line = 0; // source line, set on if/when statement nodes (for diagnostics)
			AstNode *subs[2] = {NULL, NULL};
		};

		typedef enum {
			EQN_NONE,
			EQN_ALGE,
			EQN_DIFF,
			EQN_COMP
		} EqnType;

		class AstTree {
		public:
			AstNode *pHead = NULL;
			int depth = 0;
			EqnType eqnType = EQN_NONE;
			int eqnIdx = -1;
			int line = 0; // source line of the statement (for diagnostics)

			AstTree() {
				pHead = NULL;
				depth = 0;
				eqnType = EQN_NONE;
				eqnIdx = -1;
			}

			virtual ~AstTree();

			void printTree(std::ostream &ost);

			void refreshDepth();
		};

		class IdPearl {
		public:
			unsigned char *str = NULL;
			int size = 0;
			unsigned char *annotation = NULL;
			unsigned char *unit = NULL;
			Type ty = TY_NULL;
			IdType idTy = ID_NONE;
			NumValue value = {};
			NumValue derVal = {};
			IdPearl *pair = NULL;
			bool initVal = false;
			bool initDer = false;
			int index = -1;
			int line = 0;			 // source line of the declaration (for diagnostics)
			bool isConstant = false; // declared `constant` (a parameter that cannot be changed)
			bool isInput = false;	 // declared `input`

			IdPearl(const unsigned char *s, int sz, const unsigned char *ann, const unsigned char *un, Type ty, IdType idTy);

			IdPearl(const unsigned char *s, int sz);

			IdPearl(const IdPearl &);

			void setName(const unsigned char *s, int sz);

			void setAnnotation(const unsigned char *ann);

			void setUnit(const unsigned char *un);

			virtual ~IdPearl();

		private:
			void resetName();

			void resetAnnotation();

			void resetUnit();
		};

		typedef IdPearl StrPearl;

		class SasModel {
		public:
			int errCount = 0;
			int warnCount = 0;
			std::shared_ptr<IdPearl> modelInfo = nullptr;
			std::list<std::shared_ptr<AstTree>> modelTrees;
			std::list<std::shared_ptr<IdPearl>> idTable;
			std::vector<std::list<std::shared_ptr<IdPearl>>> idNecklaces;

			std::list<std::shared_ptr<SasModel>> subModels;
			// String literals scanned for this model; AST nodes point into these.
			std::list<std::shared_ptr<StrPearl>> strTable;
			// Where error and warning messages go (default: std::cerr; null discards them).
			std::ostream *diag = &std::cerr;

			SasModel();

			void error(const char *formatMsg, ...);

			void warn(const char *formatMsg, ...);

			std::shared_ptr<IdPearl> lookupIdFromModel(const char *iden, int len);

			std::shared_ptr<IdPearl> injectNewId(const char *iden, int len);

			std::shared_ptr<IdPearl> lookupIdentifier(const char *iden, int len);

			void printAllModelTrees();

			void printIdTable();

			virtual ~SasModel();
		};

		// Owns the top-level models created for one parse/run. There is no process-wide
		// instance: every parser holds its own pool, so independent runs can proceed
		// concurrently without sharing mutable state.
		class ModelPool {
		public:
			std::list<std::shared_ptr<SasModel>> topModels;

			ModelPool() = default;

			ModelPool(const ModelPool &) = delete;
			ModelPool &operator=(const ModelPool &) = delete;

			std::shared_ptr<SasModel> createNewTopModel();

			virtual ~ModelPool() = default;

			static std::shared_ptr<SasModel> findModel(std::list<std::shared_ptr<SasModel>> &, AstNode *);
		};

	} // namespace core
} // namespace che

#endif