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
#ifndef _Che_SparseLU_H
#define _Che_SparseLU_H

#include <memory>
#include <vector>

namespace che {
	namespace core {

		// A sparse LU factorization of a square matrix (SuperLU, column ordering COLAMD, partial pivoting),
		// factored once and solved for any number of right-hand sides. The solver needs this for the
		// matrix of the algebraic equations, which is the same at every order of a series.
		//
		// This header does not include SuperLU or Armadillo (their headers do not mix well), so the matrix
		// is given in compressed column form.
		class SparseLU {
		public:
			SparseLU();
			~SparseLU();
			SparseLU(const SparseLU &) = delete;
			SparseLU &operator=(const SparseLU &) = delete;

			// Factors the n x n matrix with the column pointers colPtr (n + 1 entries), the row indices
			// rowIdx and the values of the nonzeros. Returns false if the matrix is singular (a zero pivot).
			bool factor(int n, const int *colPtr, const int *rowIdx, const double *values);

			// Solves A x = b; b (n entries) is overwritten with x. Returns false if there is no valid
			// factorization.
			bool solve(double *b);

			bool factored() const;

		private:
			struct Impl;
			std::unique_ptr<Impl> impl;
		};

	} // namespace core
} // namespace che

#endif
