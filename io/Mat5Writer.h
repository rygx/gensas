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
#ifndef IO_MAT5_WRITER_H
#define IO_MAT5_WRITER_H

#include <cstddef>
#include <fstream>
#include <string>
#include <vector>

namespace che {
	namespace io {

		// Writes MATLAB level 5 MAT files (uncompressed) without any third-party library. MATLAB, Octave
		// and SciPy read them. Supported variables: real double matrices and char matrices. A single
		// variable must be smaller than 2 GB (a limit of the format).
		class Mat5Writer {
		public:
			Mat5Writer();
			~Mat5Writer();

			// Creates the file. Returns false and sets error() on failure.
			bool open(const std::string &path);

			// rows x cols matrix; `colMajor` holds the elements column by column.
			bool writeMatrix(const std::string &name, size_t rows, size_t cols, const double *colMajor);

			// One name per row, blank padded to the longest.
			bool writeChars(const std::string &name, const std::vector<std::string> &rows);

			// Flushes and closes the file; false if anything failed since open().
			bool close();

			const std::string &error() const {
				return err;
			}

		private:
			bool fail(const std::string &msg);
			bool writeElementTag(unsigned type, size_t bytes);
			bool writePadded(const void *data, size_t bytes);
			bool writeVariable(const std::string &name, unsigned arrayClass, size_t rows, size_t cols, unsigned dataType, const void *data, size_t dataBytes);

			std::ofstream out;
			std::string err;
			bool ok = false;
		};

	} // namespace io
} // namespace che

#endif
