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
#include "io/Mat5Writer.h"

#include <algorithm>
#include <cstdint>
#include <cstring>

using namespace std;

namespace che {
	namespace io {

		namespace {
			// Data types and array classes of the MAT-file level 5 format.
			const unsigned miINT8 = 1;
			const unsigned miUINT32 = 6;
			const unsigned miINT32 = 5;
			const unsigned miDOUBLE = 9;
			const unsigned miMATRIX = 14;
			const unsigned miUTF8 = 16;
			const unsigned mxCHAR_CLASS = 4;
			const unsigned mxDOUBLE_CLASS = 6;
			const size_t kMaxVariableBytes = 0x7fffffffu;

			size_t padTo8(size_t n) {
				return (n + 7) & ~(size_t)7;
			}
		} // namespace

		Mat5Writer::Mat5Writer() {
		}

		Mat5Writer::~Mat5Writer() {
			if (out.is_open()) {
				out.close();
			}
		}

		bool Mat5Writer::fail(const string &msg) {
			if (ok) {
				err = msg;
			}
			ok = false;
			return false;
		}

		bool Mat5Writer::open(const string &path) {
			err.clear();
			out.open(path.c_str(), ios::binary | ios::trunc);
			if (!out) {
				ok = true; // let fail() record the message
				return fail("cannot create MAT file \"" + path + "\"");
			}
			ok = true;
			// 128-byte header: 116 bytes of text, 8 bytes subsystem offset, version, byte-order mark.
			char header[128];
			memset(header, ' ', sizeof(header));
			const char text[] = "MATLAB 5.0 MAT-file, written by GenSAS";
			memcpy(header, text, sizeof(text) - 1);
			memset(header + 116, 0, 8);
			uint16_t version = 0x0100;
			uint16_t indicator = 0x4D49; // 'MI' in the byte order of this machine tells readers how to decode
			memcpy(header + 124, &version, 2);
			memcpy(header + 126, &indicator, 2);
			out.write(header, sizeof(header));
			return out ? true : fail("error writing MAT file \"" + path + "\"");
		}

		bool Mat5Writer::writeElementTag(unsigned type, size_t bytes) {
			uint32_t tag[2] = {(uint32_t)type, (uint32_t)bytes};
			out.write((const char *)tag, sizeof(tag));
			return (bool)out;
		}

		// Writes `bytes` of data followed by zero padding to a multiple of 8 bytes.
		bool Mat5Writer::writePadded(const void *data, size_t bytes) {
			out.write((const char *)data, (streamsize)bytes);
			static const char zeros[8] = {0, 0, 0, 0, 0, 0, 0, 0};
			out.write(zeros, (streamsize)(padTo8(bytes) - bytes));
			return (bool)out;
		}

		bool Mat5Writer::writeVariable(const string &name, unsigned arrayClass, size_t rows, size_t cols, unsigned dataType, const void *data, size_t dataBytes) {
			if (!ok || !out.is_open()) {
				return false;
			}
			if (name.empty() || name.size() > 63) {
				return fail("MAT variable name \"" + name + "\" must have 1 to 63 characters");
			}
			if (dataBytes > kMaxVariableBytes || rows > 0x7fffffffu || cols > 0x7fffffffu) {
				return fail("MAT variable \"" + name + "\" is too large for the level 5 format");
			}
			// Contents of the miMATRIX element: array flags, dimensions, name and data, each with an
			// 8-byte tag and padded to a multiple of 8 bytes.
			const size_t contents = (8 + 8) + (8 + 8) + (8 + padTo8(name.size())) + (8 + padTo8(dataBytes));
			if (contents > 0xffffffffu) {
				return fail("MAT variable \"" + name + "\" is too large for the level 5 format");
			}
			bool good = writeElementTag(miMATRIX, contents);
			uint32_t flags[2] = {(uint32_t)arrayClass, 0};
			good = good && writeElementTag(miUINT32, sizeof(flags));
			out.write((const char *)flags, sizeof(flags));
			int32_t dims[2] = {(int32_t)rows, (int32_t)cols};
			good = good && writeElementTag(miINT32, sizeof(dims));
			out.write((const char *)dims, sizeof(dims));
			good = good && writeElementTag(miINT8, name.size());
			good = good && writePadded(name.data(), name.size());
			good = good && writeElementTag(dataType, dataBytes);
			good = good && writePadded(data, dataBytes);
			if (!good || !out) {
				return fail("error while writing MAT variable \"" + name + "\"");
			}
			return true;
		}

		bool Mat5Writer::writeMatrix(const string &name, size_t rows, size_t cols, const double *colMajor) {
			return writeVariable(name, mxDOUBLE_CLASS, rows, cols, miDOUBLE, colMajor, rows * cols * sizeof(double));
		}

		bool Mat5Writer::writeChars(const string &name, const vector<string> &rows) {
			size_t width = 0;
			for (const string &r : rows) {
				width = max(width, r.size());
			}
			vector<char> data(rows.size() * width, ' ');
			for (size_t r = 0; r < rows.size(); r++) {
				for (size_t c = 0; c < rows[r].size(); c++) {
					data[c * rows.size() + r] = rows[r][c];
				}
			}
			return writeVariable(name, mxCHAR_CLASS, rows.size(), width, miUTF8, data.data(), data.size());
		}

		bool Mat5Writer::close() {
			if (out.is_open()) {
				out.flush();
				if (!out) {
					fail("error while writing the MAT file");
				}
				out.close();
			}
			return ok;
		}

	} // namespace io
} // namespace che
