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
			// Discrete values at the start of the event instant being handled; edge() and change()
			// compare against them and are false outside event handling (`atEvent` is false).
			std::map<std::string, double> eventPre;
			bool atEvent = false;

			// Evaluates an expression. Booleans are 1.0 / 0.0. Returns false and sets `err` on failure.
			bool eval(const AstNode *p, double &out, std::string &err) const;

			// Converts a value to the declared type of discrete variable `name`.
			double coerce(const std::string &name, double v) const;
		};

		// Builds an environment from a parsed model: parameters and the initial values of the discrete variables.
		HybridEnv makeHybridEnv(const SasModel &model);

		// How abs(), min(), max() and sign() are treated in continuous equations.
		enum PiecewiseMode {
			PW_EVENTS, // exact: switching points are state events (default)
			PW_SMOOTH  // smooth approximations with parameter eps, no events
		};

		// Rewrites abs(), min(), max() and sign() calls in continuous equations. Conditions and
		// when-bodies are left alone (they are evaluated directly). Returns the number of calls
		// rewritten. Call before isHybridModel() and checkModel().
		//  - PW_EVENTS: into if-expressions (abs(e) = if e >= 0 then e else -e, ...), so that the
		//    switching points become state events and every mode is smooth.
		//  - PW_SMOOTH: into smooth functions without events, with s(e) = sqrt(e*e + eps*eps):
		//    abs(e) = s(e), max(a,b) = (a + b + s(a-b)) / 2, min(a,b) = (a + b - s(a-b)) / 2,
		//    sign(e) = e / s(e). The error is of the order eps; very small eps makes the series
		//    steps near the kink small.
		int rewritePiecewiseFunctions(SasModel &model, PiecewiseMode mode = PW_EVENTS, double eps = 1e-6);

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
			int guardId = -1; // id of the branch's condition, see collectGuards()
		};

		// A `when ... elsewhen ... end when` statement.
		struct WhenClause {
			std::vector<WhenLink> links;
		};

		// The top-level when-statements of a model.
		std::vector<WhenClause> collectWhenClauses(const SasModel &model);

		// A condition that can change its truth value during a simulation: the condition of a when
		// branch, an if-expression or an if-equation. collectGuards() numbers them (the id is stored in
		// AstNode::index of the when/if node, so it survives copying the AST) and says which ones depend
		// on continuous variables or `time` and therefore need event detection ("dynamic").
		enum GuardKind { G_IFEXPR,
						 G_IFEQ,
						 G_WHEN };

		struct Guard {
			int id = -1;
			GuardKind kind = G_IFEXPR;
			const AstNode *cond = NULL;
			int line = 0;
			bool dynamic = false;
		};

		// Numbers the guards of `model` and classifies them. `env` must hold the parameters and discrete
		// variables only (continuous variables are the names it does not know). The recognized time-event
		// conditions of when branches (`time >= e`, `sample(...)`) are not dynamic guards: they are
		// scheduled exactly by nextEventTime().
		std::vector<Guard> collectGuards(SasModel &model, const HybridEnv &env);

		// Truth values of the dynamic guards under the values in `env` (including `env.time`).
		bool evalGuards(const std::vector<Guard> &guards, const HybridEnv &env, std::map<int, bool> &truth, std::string &err);

		// Relational operators (<, <=, >, >=) in the conditions of the dynamic guards. The simulator
		// watches their margins (left minus right side) to find crossings that a guard evaluated at
		// two sample times cannot see.
		std::vector<const AstNode *> collectRelations(const std::vector<Guard> &guards);

		// Margin (left side minus right side) of a relation from collectRelations().
		bool relationMargin(const AstNode *relation, const HybridEnv &env, double &margin, std::string &err);

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
		// applied to `env`. A branch triggers if its time-event condition holds at T (only sample()
		// conditions with `initialOnly`, none with `noTimeEvents`) or if its guard id is in `fired`
		// (the guard became true). Returns false with `err` set on errors; `ev` is empty when nothing
		// fired.
		bool fireEvents(const std::vector<WhenClause> &clauses, HybridEnv &env, double T, bool initialOnly, HybridEvent &ev, std::string &err, const std::set<int> *fired = NULL, bool noTimeEvents = false);

		// Specializes `model` (a private copy) for the discrete state in `env`: when-statements are
		// removed, if-equations and if-expressions are resolved, and discrete variables are replaced by
		// their current values. What remains is a plain continuous model. Conditions of dynamic guards
		// (see collectGuards) take their truth value from `dynamicTruth` instead of being evaluated.
		bool specializeModel(SasModel &model, const HybridEnv &env, std::string &err, const std::map<int, bool> *dynamicTruth = NULL);

	} // namespace core
} // namespace che

#endif
