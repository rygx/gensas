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
#include "sas/SasResultMat.h"

#include <vector>
#include "io/Mat5Writer.h"

using namespace std;

namespace che {
	namespace core {

		namespace {
			// Row-major vector of rows to the column-major layout of MAT files.
			vector<double> columnMajor(const vector<vector<double>> &rows, size_t cols) {
				vector<double> data(rows.size() * cols);
				for (size_t r = 0; r < rows.size(); r++) {
					for (size_t c = 0; c < cols; c++) {
						data[c * rows.size() + r] = rows[r][c];
					}
				}
				return data;
			}

			class MatResultWriter : public ResultWriter {
			public:
				std::string format() const override {
					return "mat";
				}

				bool write(const SasResult &result, const std::string &path, std::string &error) const override {
					if (result.empty() || result.names.empty()) {
						error = "no results to write";
						return false;
					}
					che::io::Mat5Writer w;
					const size_t n = result.time.size();
					double version = kResultFormatVersion;
					bool ok = w.open(path) && w.writeMatrix("gensas_version", 1, 1, &version) && w.writeMatrix("t", 1, n, result.time.data());
					if (ok) {
						vector<double> s = columnMajor(result.values, n);
						ok = w.writeMatrix("s", result.names.size(), n, s.data()) && w.writeChars("names", result.names);
					}
					if (ok && !result.discreteNames.empty()) {
						vector<double> d = columnMajor(result.discrete, n);
						ok = w.writeMatrix("d", result.discreteNames.size(), n, d.data()) && w.writeChars("dnames", result.discreteNames);
					}
					if (ok && !result.events.empty()) {
						vector<double> et;
						for (const SasEventRecord &e : result.events) {
							et.push_back(e.time);
						}
						ok = w.writeMatrix("event_t", 1, et.size(), et.data());
					}
					ok = w.close() && ok;
					if (!ok) {
						error = w.error();
					}
					return ok;
				}
			};
		} // namespace

		std::unique_ptr<ResultWriter> makeMatResultWriter() {
			return std::unique_ptr<ResultWriter>(new MatResultWriter());
		}

	} // namespace core
} // namespace che
