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
#ifndef SAS_LEXICO_H
#define SAS_LEXICO_H
#include "sas/SasInput.h"
#include "sas/SasModel.h"
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

using namespace che::io;

namespace che {
	namespace core {

		struct UserFunction;

		class SasModelParser {
		public:
			BufferedFileReader reader;
			int errCount = 0;
			int warnCount = 0;
			TokenValue value;
			Coord coord;
			Coord prevCoord;
			Coord peekCoord;
			TokenValue peekValue;
			int currentToken;
			unsigned char *peekPoint;
			std::shared_ptr<IdPearl> currentId = NULL;
			std::shared_ptr<SasModel> currentModel = NULL;
			// Models created while parsing with this parser (one pool per parser, no global state).
			ModelPool pool;
			// String literals scanned while no model is attached to the parser.
			std::list<std::shared_ptr<StrPearl>> strayStrings;
			// When set, identifiers are scanned without being added to the model's id table
			// (used for attribute names such as start/fixed that are not model symbols).
			bool noInject = false;
			// Declared arrays: name -> dimension sizes. Each element is an ordinary variable named like "x[2,3]".
			std::map<std::string, std::vector<int>> arrayDims;
			// Functions defined in the source, inlined at their calls once the model has been read.
			std::map<std::string, std::shared_ptr<UserFunction>> functions;
			// Values that replace those written in the model: parameter values and start values, by name
			// ("k", or "x[2]" for an array element; an array name sets all its elements). The parser
			// applies them while it reads the declarations, so parameters computed from them follow.
			std::map<std::string, double> parameterOverrides, startOverrides;
			// The keys of the two maps that matched a declaration.
			std::set<std::string> usedOverrides;
			// Where error and warning messages go (default: std::cerr; null discards them).
			std::ostream *diag = &std::cerr;
			// Every error and warning, also when `diag` is null.
			std::vector<Diagnostic> diagnostics;

			SasModelParser();

			void parseString(const char *str);

			void parseFromFile(const char *filename);

			// void error(const char* filename, int line, const char* errorMsg);

			// void warn(const char* filename, int line, const char* errorMsg);

			void error(const char *formatMsg, ...);

			void warn(const char *formatMsg, ...);

			int bufferMargin();

			bool hasReachedBufferSize();

			void beginPeekToken();

			void endPeekToken();

			virtual ~SasModelParser();
		};

		typedef int (*scanner)(SasModelParser *);

		class Lexer {
		public:
			scanner scanners[256];
			const char *tokenStrings[NUM_TK] = {
#define TOKEN(k, n) n,
#include "sas/token.txt"
#undef TOKEN
			};

			const char *eqnNames[4]{
				"None", "Algebraic equation", "Differential equation", "Composite equation"};

			static Lexer &getInstance() {
				static Lexer instance;
				return instance;
			}

		private:
			Lexer();
		};

		int getNextToken(SasModelParser *parser);

		void expect(SasModelParser *parser, int tok);
	} // namespace core
} // namespace che

#endif