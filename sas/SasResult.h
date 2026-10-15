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
#ifndef SAS_RESULT_H
#define SAS_RESULT_H

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace che {
	namespace core {

		class SasSolutionSet;

		// An event of a hybrid simulation: the discrete variables it assigned and the continuous
		// variables it reinitialized.
		struct SasEventRecord {
			double time = 0.0;
			std::vector<std::pair<std::string, double>> assignments;
			std::vector<std::pair<std::string, double>> reinits;
		};

		// Values of the discrete variables after the events at each time in `times`; the first entry is
		// the state at the start time. values[k][i] belongs to names[i].
		struct SasDiscreteTrace {
			std::vector<std::string> names;
			std::vector<double> times;
			std::vector<std::vector<double>> values;
		};

		// Version of the result format written by all backends (see doc/ResultFormat.md).
		constexpr int kResultFormatVersion = 1;

		// Simulation result sampled on a time grid, independent of the output file format.
		// An event time appears twice in `time`: the sample before and the sample after the event.
		struct SasResult {
			double startTime = 0.0;
			double endTime = 0.0;
			double interval = 0.0;					   // requested sample spacing
			std::vector<double> time;				   // sample times
			std::vector<std::string> names;			   // continuous variables
			std::vector<std::vector<double>> values;   // values[i][j]: names[i] at time[j]
			std::vector<std::string> discreteNames;	   // hybrid models only
			std::vector<std::vector<double>> discrete; // discrete[i][j]: discreteNames[i] at time[j]
			std::vector<SasEventRecord> events;		   // hybrid models only

			bool empty() const {
				return time.empty();
			}
		};

		// Samples the solution every `interval`. Variables are matched by name across the solution
		// segments, so hybrid models whose equation structure depends on the mode are supported.
		// Returns an empty result when the set has no solutions.
		SasResult buildResult(SasSolutionSet &set, double interval);

		class ResultWriter {
		public:
			virtual ~ResultWriter() {
			}

			// Writes `result` to `path`. On failure returns false and sets `error`.
			virtual bool write(const SasResult &result, const std::string &path, std::string &error) const = 0;

			// Format name, as used by makeResultWriter.
			virtual std::string format() const = 0;
		};

		// Formats available in this library: "json" and "csv". Returns null for other names
		// (the MAT writer is in sas_result_mat_lib, see SasResultMat.h).
		std::unique_ptr<ResultWriter> makeResultWriter(const std::string &format);

	} // namespace core
} // namespace che

#endif
