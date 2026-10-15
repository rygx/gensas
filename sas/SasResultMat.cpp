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

#include <algorithm>
#include <vector>
#include "matio.h"

using namespace std;

namespace che {
	namespace core {

		namespace {
			bool writeVar(mat_t *fp, const char *name, matio_classes cls, matio_types type, size_t rows, size_t cols, void *data) {
				size_t dims[2] = {rows, cols};
				matvar_t *var = Mat_VarCreate(name, cls, type, 2, dims, data, 0);
				if (var == NULL) {
					return false;
				}
				int rc = Mat_VarWrite(fp, var, MAT_COMPRESSION_NONE);
				Mat_VarFree(var);
				return rc == 0;
			}

			bool writeMatrix(mat_t *fp, const char *name, const vector<vector<double>> &rows, size_t cols) {
				// MAT files are column-major.
				vector<double> data(rows.size() * cols);
				for (size_t r = 0; r < rows.size(); r++) {
					for (size_t c = 0; c < cols; c++) {
						data[c * rows.size() + r] = rows[r][c];
					}
				}
				return writeVar(fp, name, MAT_C_DOUBLE, MAT_T_DOUBLE, rows.size(), cols, data.data());
			}

			// Names as a char matrix, one name per row, blank padded.
			bool writeNames(mat_t *fp, const char *name, const vector<string> &names) {
				size_t width = 1;
				for (const string &n : names) {
					width = max(width, n.size());
				}
				vector<char> data(names.size() * width, ' ');
				for (size_t r = 0; r < names.size(); r++) {
					for (size_t c = 0; c < names[r].size(); c++) {
						data[c * names.size() + r] = names[r][c];
					}
				}
				return writeVar(fp, name, MAT_C_CHAR, MAT_T_UTF8, names.size(), width, data.data());
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
					mat_t *fp = Mat_CreateVer(path.c_str(), NULL, MAT_FT_MAT5);
					if (fp == NULL) {
						error = "cannot create MAT file \"" + path + "\"";
						return false;
					}
					const size_t n = result.time.size();
					bool ok = true;
					double version = kResultFormatVersion;
					ok = ok && writeVar(fp, "gensas_version", MAT_C_DOUBLE, MAT_T_DOUBLE, 1, 1, &version);
					vector<double> t = result.time;
					ok = ok && writeVar(fp, "t", MAT_C_DOUBLE, MAT_T_DOUBLE, 1, n, t.data());
					ok = ok && writeMatrix(fp, "s", result.values, n);
					ok = ok && writeNames(fp, "names", result.names);
					if (!result.discreteNames.empty()) {
						ok = ok && writeMatrix(fp, "d", result.discrete, n);
						ok = ok && writeNames(fp, "dnames", result.discreteNames);
					}
					if (!result.events.empty()) {
						vector<double> et;
						for (const SasEventRecord &e : result.events) {
							et.push_back(e.time);
						}
						ok = ok && writeVar(fp, "event_t", MAT_C_DOUBLE, MAT_T_DOUBLE, 1, et.size(), et.data());
					}
					Mat_Close(fp);
					if (!ok) {
						error = "error while writing MAT file \"" + path + "\"";
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
