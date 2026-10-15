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
#ifndef SAS_HYBRID_MODEL_H
#define SAS_HYBRID_MODEL_H

#include "sas/SasModel.h"
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace che {
	namespace core {

		// Values by name used to evaluate event and discrete expressions.
		class HybridEnv {
		public:
			std::map<std::string, double> values; // parameters, discrete variables and (at events) continuous variables
			std::map<std::string, double> pre;	  // discrete values just before the current event, for pre()
			std::map<std::string, Type> types;	  // declared types of the discrete variables
			std::set<std::string> discrete;		  // names of the discrete variables
			double time = 0.0;

			// Evaluates an expression. Booleans are 1.0 / 0.0. Returns false and sets `err` on failure.
			bool eval(const AstNode *p, double &out, std::string &err) const;

			// Converts a value to the declared type of discrete variable `name`.
			double coerce(const std::string &name, double v) const;
		};

		// Builds an environment from a parsed model: parameters and the initial values of the discrete variables.
		HybridEnv makeHybridEnv(const SasModel &model);

		// True if the model uses discrete variables, if/when statements, if-expressions or `time`.
		bool isHybridModel(const SasModel &model);

		// Time-event condition forms:
		//   TC_AT:     time >= a   (also >, and the mirrored a <= time, a < time)
		//   TC_SAMPLE: sample(a, b), true at a, a + b, a + 2b, ...
		enum TimeCondKind { TC_NONE,
							TC_AT,
							TC_SAMPLE };

		struct TimeCond {
			TimeCondKind kind = TC_NONE;
			const AstNode *a = NULL;
			const AstNode *b = NULL;
		};

		// Recognizes a time-event condition; returns false (kind TC_NONE) for any other expression.
		bool classifyTimeCondition(const AstNode *cond, TimeCond &tc);

		// One `when` / `elsewhen` branch.
		struct WhenLink {
			const AstNode *cond = NULL;
			const AstNode *body = NULL; // statement list
			int line = 0;
		};

		// A `when ... elsewhen ... end when` statement.
		struct WhenClause {
			std::vector<WhenLink> links;
		};

		// The top-level when-statements of a model.
		std::vector<WhenClause> collectWhenClauses(const SasModel &model);

		// Smallest event time strictly after `t` over all clauses; returns false if there is none.
		// Returns false with a non-empty `err` on evaluation errors.
		bool nextEventTime(const std::vector<WhenClause> &clauses, HybridEnv &env, double t, double &te, std::string &err);

		struct HybridEvent {
			double time = 0.0;
			int line = 0;
			std::vector<std::pair<std::string, double>> assignments; // discrete variables set by the event
			std::vector<std::pair<std::string, double>> reinits;	 // continuous variables reinitialized by the event
		};

		// Fires the clauses that trigger at time T: runs the body of the first triggering branch of each
		// clause, evaluating all right-hand sides with the values before the event, then applies the
		// discrete assignments to `env.values`. Continuous reinits are returned in `ev` and are not
		// applied to `env`. With `initialOnly`, only sample() conditions fire (used at the start time).
		// Returns false with `err` set on errors; `ev` is empty when nothing fired.
		bool fireEvents(const std::vector<WhenClause> &clauses, HybridEnv &env, double T, bool initialOnly, HybridEvent &ev, std::string &err);

		// Specializes `model` (a private copy) for the discrete state in `env`: when-statements are
		// removed, if-equations and if-expressions are resolved, and discrete variables are replaced by
		// their current values. What remains is a plain continuous model.
		bool specializeModel(SasModel &model, const HybridEnv &env, std::string &err);

	} // namespace core
} // namespace che

#endif
