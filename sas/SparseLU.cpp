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
#include "sas/SparseLU.h"

#include "slu_ddefs.h"

namespace che {
	namespace core {

		struct SparseLU::Impl {
			int n = 0;
			bool ok = false;
			SuperMatrix L, U;
			int *permR = nullptr;
			int *permC = nullptr;
			SuperLUStat_t stat;

			void release() {
				if (ok) {
					Destroy_SuperNode_Matrix(&L);
					Destroy_CompCol_Matrix(&U);
					ok = false;
				}
				if (permR != nullptr) {
					SUPERLU_FREE(permR);
					permR = nullptr;
				}
				if (permC != nullptr) {
					SUPERLU_FREE(permC);
					permC = nullptr;
				}
			}
		};

		SparseLU::SparseLU() : impl(new Impl()) {}

		SparseLU::~SparseLU() {
			impl->release();
		}

		bool SparseLU::factored() const {
			return impl->ok;
		}

		bool SparseLU::factor(int n, const int *colPtr, const int *rowIdx, const double *values) {
			Impl &s = *impl;
			s.release();
			s.n = n;
			if (n <= 0) {
				return false;
			}
			const int nnz = colPtr[n];

			// SuperLU takes over the arrays of the matrix (Destroy_CompCol_Matrix frees them).
			double *a = doubleMalloc(nnz);
			int_t *asub = intMalloc(nnz);
			int_t *xa = intMalloc(n + 1);
			if (a == nullptr || asub == nullptr || xa == nullptr) {
				return false;
			}
			for (int k = 0; k < nnz; k++) {
				a[k] = values[k];
				asub[k] = rowIdx[k];
			}
			for (int j = 0; j <= n; j++) {
				xa[j] = colPtr[j];
			}
			SuperMatrix A, AC;
			dCreate_CompCol_Matrix(&A, n, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

			s.permR = (int *)SUPERLU_MALLOC(n * sizeof(int));
			s.permC = (int *)SUPERLU_MALLOC(n * sizeof(int));
			int *etree = (int *)SUPERLU_MALLOC(n * sizeof(int));
			superlu_options_t options;
			set_default_options(&options);
			options.ColPerm = COLAMD;
			get_perm_c(options.ColPerm, &A, s.permC);
			StatInit(&s.stat);
			sp_preorder(&options, &A, s.permC, etree, &AC);
			GlobalLU_t glu;
			int_t info = 0;
			dgstrf(&options, &AC, sp_ienv(2), sp_ienv(1), etree, nullptr, 0, s.permC, s.permR, &s.L, &s.U, &glu, &s.stat, &info);
			SUPERLU_FREE(etree);
			Destroy_CompCol_Permuted(&AC);
			Destroy_CompCol_Matrix(&A);
			StatFree(&s.stat);
			if (info != 0) {
				// info > 0: a zero pivot at that column (L and U are not valid then); info < 0: bad argument
				if (info > 0) {
					if (info <= n) {
						Destroy_SuperNode_Matrix(&s.L);
						Destroy_CompCol_Matrix(&s.U);
					}
				}
				if (s.permR != nullptr) {
					SUPERLU_FREE(s.permR);
					s.permR = nullptr;
				}
				if (s.permC != nullptr) {
					SUPERLU_FREE(s.permC);
					s.permC = nullptr;
				}
				return false;
			}
			s.ok = true;
			return true;
		}

		bool SparseLU::solve(double *b) {
			Impl &s = *impl;
			if (!s.ok) {
				return false;
			}
			SuperMatrix B;
			dCreate_Dense_Matrix(&B, s.n, 1, b, s.n, SLU_DN, SLU_D, SLU_GE);
			StatInit(&s.stat);
			int_t info = 0;
			dgstrs(NOTRANS, &s.L, &s.U, s.permC, s.permR, &B, &s.stat, &info);
			StatFree(&s.stat);
			Destroy_SuperMatrix_Store(&B);
			return info == 0;
		}

	} // namespace core
} // namespace che
