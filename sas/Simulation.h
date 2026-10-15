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
#ifndef SAS_SIMULATION_H
#define SAS_SIMULATION_H

#include "sas/HybridModel.h"
#include "sas/SasComputation.h"
#include "sas/SasResult.h"

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace che {
	namespace core {

		// Everything needed to run one Modelica simulation. Filled from a job file (parseModelicaJob) or
		// directly by the caller. Limits on the numeric options are applied by Simulation (see
		// applyJobLimits); out-of-range values are clamped and reported as warnings.
		struct ModelicaJob {
			// The model: a file, or the Modelica text itself when `modelFile` is empty.
			std::string modelFile;
			std::string modelSource;

			SasComputationOptions options;
			PiecewiseMode piecewise = PW_EVENTS;
			double smoothEps = 1e-6;
			double outputStep = 0.01; // spacing of the samples in the result

			// Output files requested by a job file (already resolved against its directory). Simulation
			// does not write them; the caller does, see Simulation::writeResult.
			std::string jsonOutput;
			std::string csvOutput;
			std::string matOutput;

			// Print the model trees and solver diagnostics to stdout (debugging; not for concurrent runs,
			// the tree printer redirects std::cout).
			bool verbose = false;

			ModelicaJob();
		};

		// Parses a job file (JSON, see doc/JobFormat.md). Relative paths are resolved against `baseDir`
		// (the directory of the job file; may be empty). Unknown keys, wrong types and unsupported
		// versions or kinds are errors whose message names the key path.
		bool parseModelicaJob(const std::string &json, const std::string &baseDir, ModelicaJob &job, std::string &error);

		// Applies "dotted.key=value" assignments to the job text (value is JSON; a value that is not valid
		// JSON is taken as a string), e.g. "solver.endTime=5" or "model.file=other.mo".
		bool applyJobOverrides(std::string &json, const std::vector<std::string> &assignments, std::string &error);

		// Clamps the numeric options of `job` to the supported ranges; one warning per change is appended.
		void applyJobLimits(ModelicaJob &job, std::vector<std::string> &warnings);

		// Reports progress after each solution segment (time reached, end time); return true to cancel.
		typedef std::function<bool(double time, double endTime)> ProgressFn;

		// One simulation: parse, check, solve, sample. Independent Simulation objects may run on
		// different threads at the same time; one object must not be used from several threads.
		class Simulation {
		public:
			explicit Simulation(const ModelicaJob &job);

			// Runs the job. Returns false on error (see error()) or if `progress` asked to cancel
			// (see cancelled()). Never throws: failures inside the solver are reported through error().
			bool run(const ProgressFn &progress = ProgressFn());

			const std::string &error() const {
				return err;
			}

			bool cancelled() const {
				return wasCancelled;
			}

			// Messages about adjusted options and similar; valid right after construction.
			const std::vector<std::string> &warnings() const {
				return warns;
			}

			// The effective job (after limits were applied).
			const ModelicaJob &job() const {
				return jobData;
			}

			// The sampled result, valid after a successful run().
			const SasResult &result() const {
				return res;
			}

			// Wall-clock seconds spent in the last run().
			double elapsedSeconds() const {
				return elapsed;
			}

			// Writes the result as "json", "csv" or "mat". Returns false and sets `error` on failure.
			bool writeResult(const std::string &format, const std::string &path, std::string &error) const;

		private:
			bool runImpl(const ProgressFn &progress);

			ModelicaJob jobData;
			std::vector<std::string> warns;
			std::string err;
			bool wasCancelled = false;
			double elapsed = 0.0;
			SasResult res;
		};

	} // namespace core
} // namespace che

#endif
