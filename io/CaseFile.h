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
#ifndef IO_CASE_FILE_H
#define IO_CASE_FILE_H

#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace che {
	namespace io {

		// A real double matrix, stored column by column.
		struct Matrix {
			int rank = 2;
			size_t dims[2] = {0, 0}; // rows, columns
			std::vector<double> values;
		};

		// Named matrices of a power system case (see kPsatMatrixNames).
		class MatrixSet {
		public:
			void set(const std::string &name, Matrix m) {
				items[name] = std::move(m);
			}

			// Null if the case does not contain the matrix.
			const Matrix *find(const std::string &name) const {
				std::map<std::string, Matrix>::const_iterator it = items.find(name);
				return it == items.end() ? NULL : &it->second;
			}

			const std::map<std::string, Matrix> &all() const {
				return items;
			}

		private:
			std::map<std::string, Matrix> items;
		};

		// The PSAT data matrices the solver reads (bus is mandatory, the others are optional).
		extern const char *const kPsatMatrixNames[];
		extern const size_t kPsatMatrixCount;

		// GenSAS case file, version 1 (see doc/CaseFormat.md). Two encodings of the same content:
		//  - binary: "GSCASE01", the header length, a JSON header and the matrices as little-endian
		//    float64 in column-major order (fast, for large cases);
		//  - JSON text, with the matrices as arrays of rows (readable, for small cases).
		// readCaseFile recognizes the encoding from the first bytes.
		bool readCaseFile(const std::string &path, MatrixSet &out, std::string &error);

		bool writeCaseFile(const std::string &path, const MatrixSet &in, std::string &error);

		// JSON text encoding.
		bool writeCaseJson(const std::string &path, const MatrixSet &in, std::string &error);

	} // namespace io
} // namespace che

#endif
