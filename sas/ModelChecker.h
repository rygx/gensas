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

#ifndef SAS_MODEL_CHECKER_H
#define SAS_MODEL_CHECKER_H

#include "sas/SasModel.h"
#include <iosfwd>
#include <string>

namespace che {
	namespace core {

		// Language features the parser accepts but that only some solver versions can run.
		// A feature is reported as an error by checkModel() unless its flag is set. The
		// default (all false) describes the current solver: continuous Real variables,
		// parameters, `der()` equations and algebraic equations.
		struct SolverFeatures {
			bool discreteVariables = false;	 // Integer, Boolean, discrete Real
			bool ifExpressions = false;		 // if c then a else b
			bool ifEquations = false;		 // if c then <equations> else ... end if
			bool whenEquations = false;		 // when c then ... end when
			bool piecewiseFunctions = false; // abs, min, max, sign, floor
			bool eventFunctions = false;	 // pre, edge, change, sample, noEvent, reinit
			bool simulationTime = false;	 // `time`
			bool terminateAndAssert = false; // terminate(), assert()
			// Restricts hybrid constructs to what the time-event simulator can run: when-conditions are
			// `time >= e` or `sample(start, interval)`, if-conditions depend only on parameters and
			// discrete variables (no state events), and discrete-only operators stay out of
			// continuous equations. Use together with the flags above.
			bool timeEventsOnly = false;
		};

		// Checks a parsed model: undeclared identifiers, unknown functions, and use of
		// features that `features` does not enable. Diagnostics are written to `err` as
		// "<filename>:<line>: error: <message>". Returns the number of errors.
		int checkModel(SasModel &model, const SolverFeatures &features, const std::string &filename, std::ostream &err);

	} // namespace core
} // namespace che

#endif
