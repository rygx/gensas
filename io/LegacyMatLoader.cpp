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
#include "io/LegacyMatLoader.h"

#include "matio.h"

using namespace std;

namespace che {
	namespace io {

		bool loadPsatMatFile(const string &path, MatrixSet &out, string &error) {
			mat_t *fp = Mat_Open(path.c_str(), MAT_ACC_RDONLY);
			if (fp == NULL) {
				error = "cannot open MAT file \"" + path + "\"";
				return false;
			}
			for (size_t i = 0; i < kPsatMatrixCount; i++) {
				matvar_t *v = Mat_VarRead(fp, kPsatMatrixNames[i]);
				if (v == NULL) {
					continue;
				}
				if (v->rank != 2 || v->class_type != MAT_C_DOUBLE || v->isComplex) {
					error = string("\"") + kPsatMatrixNames[i] + "\" is not a real double matrix";
					Mat_VarFree(v);
					Mat_Close(fp);
					return false;
				}
				Matrix m;
				m.dims[0] = v->dims[0];
				m.dims[1] = v->dims[1];
				const double *data = (const double *)v->data;
				m.values.assign(data, data + v->dims[0] * v->dims[1]);
				out.set(kPsatMatrixNames[i], std::move(m));
				Mat_VarFree(v);
			}
			Mat_Close(fp);
			return true;
		}

	} // namespace io
} // namespace che
