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
#ifndef SAS_HYBRID_SIMULATOR_H
#define SAS_HYBRID_SIMULATOR_H

#include "sas/HybridModel.h"
#include "sas/SasComputation.h"
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace che {
	namespace core {

		// Simulates a flat hybrid model with time and state events. Between events the model is
		// specialized for the current discrete state (see specializeModel) and solved with the
		// semi-analytical solver; at each event the when-equations run, discrete variables and continuous
		// states are updated, and the next segment starts from the new state.
		//
		// Time events (`time >= e`, `sample(...)`) are scheduled exactly. State events are conditions on
		// continuous variables or `time`; they are detected by sampling each solution segment and are
		// located by bisection on the Pade solution. A condition that flips true -> false needs no
		// action; a when-condition that flips false -> true fires, and an if-condition that flips
		// switches the mode of the model for the next segment. Crossings that start and end within one
		// sample interval of a segment are not detected (see setGuardSamples).
		//
		// The model must have passed checkModel() with `timeEventsOnly` (and `stateEvents`) set.
		class HybridSimulator {
		public:
			explicit HybridSimulator(std::shared_ptr<SasModel> model);

			// Runs from options.startTime to options.endTime. Returns the solutions (owned by the caller;
			// consecutive segments share their boundary time) or NULL on error, see error().
			SasSolutionSet *run(const SasComputationOptions &options);

			const std::string &error() const { return err; }

			// True if options.progress stopped the run (run() then returned NULL with error() "cancelled").
			bool cancelled() const { return wasCancelled; }

			// Events that fired, in time order. Events at the start time come first.
			const std::vector<HybridEvent> &events() const { return firedEvents; }

			// Values of the discrete and continuous variables at the end time, by name.
			const std::map<std::string, double> &finalValues() const { return finalVals; }

			// Warnings of the last run: assert() with AssertionLevel.warning that failed.
			const std::vector<std::string> &warnings() const { return warningList; }

			// True if terminate() ended the last run; the solutions then stop at terminationTime().
			bool terminated() const { return wasTerminated; }
			double terminationTime() const { return terminateTime; }
			const std::string &terminationMessage() const { return terminateMessage; }

			// Number of solver segments (between events) that were run.
			int segmentCount() const { return nSegments; }

			// Number of segments that reused the flattened model of an earlier segment (SasComputationOptions::modeCache).
			int modeCacheHits() const { return nCacheHits; }

			// Number of samples per solution segment used to detect state events (default 48).
			void setGuardSamples(int n) { guardSamples = n < 2 ? 2 : n; }

		private:
			std::shared_ptr<SasModel> model;
			std::string err;
			std::vector<HybridEvent> firedEvents;
			std::map<std::string, double> finalVals;
			int nSegments = 0;
			int nCacheHits = 0;
			bool wasCancelled = false;
			std::vector<std::string> warningList;
			bool wasTerminated = false;
			double terminateTime = 0.0;
			std::string terminateMessage;
			int guardSamples = 48;

			bool fail(const std::string &msg);
		};

	} // namespace core
} // namespace che

#endif
