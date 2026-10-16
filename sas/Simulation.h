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

		class SasModelParser;
		struct SolverFeatures;

		// What an editor needs to know about a model, found by Simulation::check() without running it.
		struct ModelVariable {
			std::string name;
			std::string kind; // "state" (der() is used), "algebraic" or "discrete"
			std::string type; // "Real", "Integer" or "Boolean"
			std::string unit;
			std::string description;
			bool hasStart = false; // a start value was given
			double start = 0.0;
			bool input = false;
			int line = 0;
		};

		struct ModelParameter {
			std::string name;
			std::string type;
			std::string unit;
			std::string description;
			double value = 0.0;
			bool constant = false; // declared `constant`: it cannot be changed
			int line = 0;
		};

		struct ModelInfo {
			std::string name;
			std::string description;
			bool hybrid = false; // has events (discrete variables, if/when, time)
			std::vector<ModelVariable> variables;
			std::vector<ModelParameter> parameters;
		};

		// The model info as JSON text: {"name", "description", "hybrid", "variables": [...], "parameters": [...]}.
		std::string modelInfoJson(const ModelInfo &info);

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

			// Messages about adjusted options (known right after construction) and about the last run.
			const std::vector<std::string> &warnings() const {
				return warns;
			}

			// Parses and checks the model without running it. Returns true if run() would get past this
			// stage. Errors and warnings about the text are in diagnostics(), the first error as text in
			// error().
			bool check();

			// Variables and parameters of the model, valid after a check() or run() that read the model
			// (even when it then failed a later check).
			bool hasModelInfo() const {
				return infoValid;
			}
			const ModelInfo &modelInfo() const {
				return info;
			}

			// Errors and warnings about the model text found by the parser and the checker in the last
			// run() or check(), with line numbers (empty if the model could not even be read).
			const std::vector<Diagnostic> &diagnostics() const {
				return diags;
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
			// Reads and checks the model; `hybrid` and `features` are what the solver stage needs.
			bool frontEnd(SasModelParser &parser, bool &hybrid, SolverFeatures &features);
			// The numbers of states and algebraic variables must match those of the equations.
			bool checkDaeShape(const SasComputationModel &comp);

			ModelicaJob jobData;
			std::vector<std::string> warns;
			std::vector<Diagnostic> diags;
			ModelInfo info;
			bool infoValid = false;
			size_t nJobWarnings = 0;
			std::string err;
			bool wasCancelled = false;
			double elapsed = 0.0;
			SasResult res;
		};

	} // namespace core
} // namespace che

#endif
