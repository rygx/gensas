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
#include "util/AbstractCheCalculator.h"
#include "util/CheCompUtil.h"

#include <vector>

using namespace che::util;
using namespace che::io;

namespace che {
	namespace core {
		static vec getPowerSeriesValue(const mat &coeff, double alpha) {
			vec sol = coeff.col(coeff.n_cols - 1);
			for (int i = coeff.n_cols - 2; i >= 0; i--) {
				sol = alpha * sol + coeff.col(i);
			}
			return sol;
		}

		static mat solveToepLU(const mat &ct, const mat &y) {
			int d = ct.n_rows;
			int n = y.n_cols;

			mat x(d, n, fill::zeros);

			for (int i = 0; i < d; i++) {
				mat toep(n, n);
				for (int j = -n + 1; j < n; j++) {
					toep.diag(j).fill(ct(i, n - 1 - j));
				}
				// toep.print("toep:");
				try {
					x.row(i) = solve(toep, y.row(i).t(), solve_opts::no_approx).t();
				} catch (const std::runtime_error &e) {
					x.row(i).fill(datum::nan);
				}
			}
			return x;
		}

		// Solves one Toeplitz system per row with Levinson's recursion. The rows are independent, so each one is
		// done on plain arrays; the same arithmetic on matrix expressions allocated a dozen temporaries per
		// step of the recursion and dominated the run time of small models.
		static mat solveToepLevinson(const mat &ct, const mat &y) {
			const int d = ct.n_rows;
			const int n = y.n_cols;

			mat x(d, n, fill::zeros);
			if (n <= 0 || d <= 0) {
				return x;
			}
			std::vector<double> f(n), b(n), temp(n), cr(ct.n_cols), xr(n), yr(n);

			for (int r = 0; r < d; r++) {
				for (int j = 0; j < (int)ct.n_cols; j++) {
					cr[j] = ct(r, j);
				}
				for (int j = 0; j < n; j++) {
					yr[j] = y(r, j);
					f[j] = 0.0;
					b[j] = 0.0;
					temp[j] = 0.0;
					xr[j] = 0.0;
				}
				f[0] = 1.0 / cr[n - 1];
				b[n - 1] = f[0];
				xr[0] = yr[0] * f[0];

				for (int i = 1; i < n; i++) {
					double epsf = 0.0, epsb = 0.0, epsx = 0.0;
					for (int k = 0; k < i; k++) {
						epsf += f[k] * cr[n + i - 1 - k];
						epsb += b[n - i + k] * cr[n - 2 - k];
						epsx += xr[k] * cr[n + i - 1 - k];
					}
					const double alphaf = 1.0 / (1.0 - epsf * epsb);
					const double betaf = -epsf * alphaf;
					const double alphab = -epsb * alphaf;
					const double betab = alphaf;

					for (int k = 0; k < i; k++) {
						temp[k] = f[k];
					}
					temp[i] = 0.0;
					for (int k = 0; k <= i; k++) {
						const double bk = b[n - i - 1 + k];
						f[k] = temp[k] * alphaf + bk * betaf;
						b[n - i - 1 + k] = temp[k] * alphab + bk * betab;
					}
					const double ex = yr[i] - epsx;
					for (int k = 0; k <= i; k++) {
						xr[k] += b[n - i - 1 + k] * ex;
					}
				}
				for (int j = 0; j < n; j++) {
					x(r, j) = xr[j];
				}
			}

			return x;
		}

		CheCompOptions::CheCompOptions(int nLvl, double maxAlpha, double alphaTol, double segLen, double diffTol, double diffTolMax) {
			this->nLvl = nLvl;
			this->maxAlpha = maxAlpha;
			this->alphaTol = alphaTol;
			this->segLen = segLen;
			this->diffTol = diffTol;
			this->diffTolMax = diffTolMax;
		}

		CheSolution::CheSolution(int nState, int nLvl) : solution(nState, nLvl, fill::zeros) {
			this->type = CHESOL_NONE;
			this->nState = nState;
			this->nLvl = nLvl;
		}

		CheSolution::~CheSolution() {}

		CheSolutionPowerSeries::CheSolutionPowerSeries(int nState, int nLvl) : CheSolution(nState, nLvl) {
			type = CHESOL_PS;
		}

		vec CheSolutionPowerSeries::getSolValue(double alpha) {
			return getPowerSeriesValue(solution, alpha);
		}

		CheSolutionPowerSeries::~CheSolutionPowerSeries() {}

		CheSolutionPade::CheSolutionPade(int nState, int num, int den, PadeSolverType sol) : CheSolution(nState, num + den) {
			this->num = num;
			this->den = den;
			ready = false;
			solver = sol;
			this->type = CHESOL_PADE;
		}

		CheSolutionPade::CheSolutionPade(int nState, int nLvl, PadeSolverType sol)
			: CheSolutionPade(nState, nLvl - nLvl / 2, nLvl / 2, sol) {}

		vec CheSolutionPade::getSolValue(double alpha) {
			if (!ready) {
				if (solver == PADE_LEVINSON) {
					ready = genPadeCoeffLevinson();
				} else {
					ready = genPadeCoeffLU();
				}
			}
			if (!ready) {
				return getPowerSeriesValue(solution, alpha);
			} else {
				// Horner's rule on each row: the numerator, and the denominator 1 + d0*a + d1*a^2 + ...
				// The loop over the rows is the inner one, over contiguous memory (column-major storage).
				vec sol(nState);
				vec den(nState);
				const int nn = numerator.n_cols;
				const int nd = denomenator.n_cols;
				double *s = sol.memptr();
				double *dv = den.memptr();
				const double *np = numerator.memptr();
				const double *dp = denomenator.memptr();
				for (int r = 0; r < nState; r++) {
					s[r] = nn > 0 ? np[r + (size_t)(nn - 1) * nState] : 0.0;
					dv[r] = nd > 0 ? dp[r + (size_t)(nd - 1) * nState] : 1.0;
				}
				for (int k = nn - 2; k >= 0; k--) {
					const double *col = numerator.colptr(k);
					for (int r = 0; r < nState; r++) {
						s[r] = alpha * s[r] + col[r];
					}
				}
				for (int k = nd - 2; k >= 0; k--) {
					const double *col = denomenator.colptr(k);
					for (int r = 0; r < nState; r++) {
						dv[r] = alpha * dv[r] + col[r];
					}
				}
				bool finite = true;
				for (int r = 0; r < nState; r++) {
					if (nd > 0) {
						dv[r] = alpha * dv[r] + 1.0;
					}
					s[r] /= dv[r];
					finite = finite && std::isfinite(s[r]);
				}
				if (!finite) {
					uvec nfrows = find_nonfinite(sol);
					sol.rows(nfrows) = getPowerSeriesValue(solution.rows(nfrows), alpha);
					numerator.rows(nfrows) = solution.submat(nfrows, regspace<uvec>(0, this->num - 1));
					denomenator.rows(nfrows).fill(0.0);
				}
				return sol;
			}
		}

		// Derivative of the power series sum c_k a^k, rows 0..nRows-1, by Horner's rule on k c_k.
		static vec powerSeriesDerivative(const mat &c, int nRows, double alpha) {
			vec d(nRows, fill::zeros);
			for (int k = static_cast<int>(c.n_cols) - 1; k >= 1; k--) {
				for (int r = 0; r < nRows; r++) {
					d[r] = d[r] * alpha + k * c(r, k);
				}
			}
			return d;
		}

		vec CheSolutionPowerSeries::getDerValue(double alpha, int nRows) {
			if (nRows < 0 || nRows > nState) {
				nRows = nState;
			}
			return powerSeriesDerivative(solution, nRows, alpha);
		}

		vec CheSolutionPade::getDerValue(double alpha, int nRows) {
			if (nRows < 0 || nRows > nState) {
				nRows = nState;
			}
			if (!ready) {
				if (solver == PADE_LEVINSON) {
					ready = genPadeCoeffLevinson();
				} else {
					ready = genPadeCoeffLU();
				}
			}
			if (!ready) {
				return powerSeriesDerivative(solution, nRows, alpha);
			}
			// P = N / D with D = 1 + a Q(a); P' = (N' D - N D') / D^2. Horner's rule for each polynomial and its
			// derivative, the loop over the rows innermost on contiguous memory.
			const int nn = numerator.n_cols;
			const int nd = denomenator.n_cols;
			vec out(nRows);
			vec n(nRows, fill::zeros), dn(nRows, fill::zeros), q(nRows, fill::zeros), dq(nRows, fill::zeros);
			for (int k = nn - 1; k >= 0; k--) {
				const double *col = numerator.colptr(k);
				for (int r = 0; r < nRows; r++) {
					dn[r] = dn[r] * alpha + n[r];
					n[r] = n[r] * alpha + col[r];
				}
			}
			for (int k = nd - 1; k >= 0; k--) {
				const double *col = denomenator.colptr(k);
				for (int r = 0; r < nRows; r++) {
					dq[r] = dq[r] * alpha + q[r];
					q[r] = q[r] * alpha + col[r];
				}
			}
			bool finite = true;
			for (int r = 0; r < nRows; r++) {
				const double den = 1.0 + alpha * q[r];
				const double dden = q[r] + alpha * dq[r];
				out[r] = (dn[r] * den - n[r] * dden) / (den * den);
				finite = finite && std::isfinite(out[r]);
			}
			if (!finite) {
				const vec ps = powerSeriesDerivative(solution, nRows, alpha);
				for (int r = 0; r < nRows; r++) {
					if (!std::isfinite(out[r])) {
						out[r] = ps[r];
					}
				}
			}
			return out;
		}

		bool CheSolutionPade::genPadeCoeffLU() {
			mat augTc = join_rows(mat(nState, den, fill::zeros), solution);
			int nT = 2 * den - 1;
			mat ttxc = augTc.cols(augTc.n_cols - nT - 1, augTc.n_cols - 2);
			mat ytxc = augTc.tail_cols(den);

			this->denomenator = solveToepLU(ttxc, ytxc);
			this->numerator = augTc.cols(den, num + den - 1);

			for (int i = 0; i < den; i++) {
				this->numerator += augTc.cols(den - i - 1, num + den - i - 2).each_col() % this->denomenator.col(i);
			}
			return true;
		}

		bool CheSolutionPade::genPadeCoeffLevinson() {
			mat augTc = join_rows(mat(nState, den, fill::zeros), solution);
			int nT = 2 * den - 1;
			mat ttxc = augTc.cols(augTc.n_cols - nT - 1, augTc.n_cols - 2);
			mat ytxc = augTc.tail_cols(den);

			this->denomenator = -solveToepLevinson(ttxc, ytxc);
			this->numerator = augTc.cols(den, num + den - 1);

			for (int i = 0; i < den; i++) {
				for (int j = 0; j < num; j++) {
					for (int r = 0; r < nState; r++) {
						this->numerator(r, j) += augTc(r, den - i - 1 + j) * this->denomenator(r, i);
					}
				}
			}
			return true;
		}

		CheSolutionPade::~CheSolutionPade() {}

		CheSingleEmbedSystem::CheSingleEmbedSystem(const CheState &init, const chedata::PsatDataSet &baseSys, double startAlpha = 0)
			: initState(init), baseSys(baseSys) {
			this->startAlpha = startAlpha;
			yMatrix = CheCompUtil::getCheYMatrix(baseSys);
		}

		CheSolution *CheSolutionFactory::makeInitCheSol(int type, int nState, int nLvl) {
			CheSolution *pSol = NULL;
			if (type == CHESOL_PS) {
				pSol = new CheSolutionPowerSeries(nState, nLvl);
			} else if (type == CHESOL_PADE) {
				pSol = new CheSolutionPade(nState, nLvl, PADE_LEVINSON);
			}
			return pSol;
		}

		CheSolution *CheSolutionFactory::makeCopyCheSol(CheSolution *sol) {
			if (sol == NULL)
				return NULL;
			if (sol->type == CHESOL_PS) {
				CheSolutionPowerSeries *pSol = new CheSolutionPowerSeries(sol->nState, sol->nLvl);
				CheSolutionPowerSeries *pOri = dynamic_cast<CheSolutionPowerSeries *>(sol);
				if (pOri != nullptr) {
					pSol->solution = pOri->solution;
				}
				return pSol;
			} else if (sol->type == CHESOL_PADE) {
				CheSolutionPade *pSol = new CheSolutionPade(sol->nState, sol->nLvl, PADE_LEVINSON);
				CheSolutionPade *pOri = dynamic_cast<CheSolutionPade *>(sol);
				if (pOri != nullptr) {
					pSol->num = pOri->num;
					pSol->den = pOri->num;
					pSol->numerator = pOri->numerator;
					pSol->denomenator = pOri->denomenator;
					pSol->ready = pOri->ready;
				}
				return pSol;
			}
			return NULL;
		}

		mat getDer(const mat &c) {
			mat d(c.n_rows, c.n_cols, fill::zeros);
			for (int i = 1; i < c.n_cols; i++) {
				d.col(i - 1) = i * c.col(i);
			}
			return d;
		}

		CheSolution *CheSolutionFactory::makeDerivativeSol(CheSolution *sol, int nRows) {
			if (sol == NULL)
				return NULL;
			if (nRows < 0 || nRows > sol->nState) {
				nRows = sol->nState;
			}
			if (sol->type == CHESOL_PS) {
				CheSolutionPowerSeries *pSol = new CheSolutionPowerSeries(nRows, sol->nLvl);
				CheSolutionPowerSeries *pOri = dynamic_cast<CheSolutionPowerSeries *>(sol);
				if (pOri != nullptr) {
					pSol->solution = getDer(pOri->solution.head_rows(nRows));
				}
				return pSol;
			} else if (sol->type == CHESOL_PADE) {
				CheSolutionPade *pSol = new CheSolutionPade(nRows, sol->nLvl, PADE_LEVINSON);
				CheSolutionPade *pOri = dynamic_cast<CheSolutionPade *>(sol);
				if (pOri != nullptr) {
					pSol->solution = getDer(pOri->solution.head_rows(nRows));
					pSol->num = pOri->num;
					pSol->den = pOri->num;
					pSol->ready = false;
				}
				return pSol;
			}
			return NULL;
		}

		AbstractCheCalculator::AbstractCheCalculator(
			const chedata::PsatDataSet &sys, const CheCompOptions &compOpt)
			: baseSys(regulateIsland(sys)), compOpt(compOpt) {
			this->cheList = list<CheSingleEmbedSystem *>();
			this->solList = list<CheSolution *>();
			this->reachesMaxAlpha = false;
		}

		chedata::PsatDataSet AbstractCheCalculator::regulateIsland(const chedata::PsatDataSet &sys) {
			if (!sys.isFormatted) {
				chedata::PsatDataSet newSys(sys);
				newSys.renumberBuses();
				return newSys;
			} else {
				return sys;
			}
		}

		int AbstractCheCalculator::calc() {
			double alpha = 0;
			CheSingleEmbedSystem *pinit = getInitSystem(baseSys);
			cheList.push_back(pinit);
			CheSingleEmbedSystem *currSys;
			while ((currSys = cheList.back()) != NULL) {
				CheSolution *cheSol = getCheSolution();
				solList.push_back(cheSol);
				if (cheSol == NULL) {
					break;
				}
				cheList.push_back(getNewStage());
				if (reachesMaxAlpha) {
					break;
				}
			}
			return 0;
		}

		CheState AbstractCheCalculator::exportResult() {
			return cheList.back()->initState;
		}

		AbstractCheCalculator::~AbstractCheCalculator() {
			for (auto &&che : cheList) {
				if (che != NULL)
					delete che;
			}
			cheList.clear();
			for (auto &&sol : solList) {
				if (sol != NULL)
					delete sol;
			}
			solList.clear();
		}
	} // namespace core
} // namespace che