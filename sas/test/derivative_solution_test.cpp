// CheSolutionFactory::makeDerivativeSol with a number of rows: the derivative of the first rows of a
// solution equals the first rows of the derivative of the whole solution.
#include "sas/test/test_util.h"
#include "util/AbstractCheCalculator.h"

#include <cmath>
#include <memory>

using namespace che::core;

// The series of exp(a t) + b t^2 in every row: smooth, so the Pade approximants are accurate.
static void fillSeries(CheSolution *sol, int nState, int nLvl) {
	sol->solution.set_size(nState, nLvl);
	for (int r = 0; r < nState; r++) {
		const double a = -0.7 + 0.4 * r;
		double term = 1.0;
		for (int k = 0; k < nLvl; k++) {
			sol->solution(r, k) = term + (k == 2 ? 0.3 * (r + 1) : 0.0);
			term *= a / (k + 1);
		}
	}
}

static void check(int type) {
	const int nState = 5, nLvl = 15, nRows = 3;
	std::unique_ptr<CheSolution> sol(CheSolutionFactory::makeInitCheSol(type, nState, nLvl));
	fillSeries(sol.get(), nState, nLvl);
	std::unique_ptr<CheSolution> all(CheSolutionFactory::makeDerivativeSol(sol.get()));
	std::unique_ptr<CheSolution> head(CheSolutionFactory::makeDerivativeSol(sol.get(), nRows));
	CHECK(all != nullptr && head != nullptr);
	if (all == nullptr || head == nullptr) {
		return;
	}
	CHECK_EQ(all->nState, nState);
	CHECK_EQ(head->nState, nRows);
	for (double alpha : {0.0, 0.05, 0.3, 0.9}) {
		const arma::vec a = all->getSolValue(alpha);
		const arma::vec h = head->getSolValue(alpha);
		CHECK_EQ((int)a.n_rows, nState);
		CHECK_EQ((int)h.n_rows, nRows);
		for (int r = 0; r < nRows && r < (int)h.n_rows; r++) {
			CHECK_EQ(a(r), h(r));
		}
	}
	// the derivative of exp(a t) + 0.3 (r + 1) t^2 at t = 0.3 for the first row: a exp(a t) + 0.6 t
	const double d0 = head->getSolValue(0.3)(0);
	CHECK_NEAR(d0, -0.7 * std::exp(-0.7 * 0.3) + 0.6 * 0.3, 1e-6);
	// a number of rows out of range means all rows
	std::unique_ptr<CheSolution> big(CheSolutionFactory::makeDerivativeSol(sol.get(), 99));
	CHECK_EQ(big->nState, nState);
}

int main() {
	check(CHESOL_PADE);
	check(CHESOL_PS);
	return sastest::finish("derivative_solution_test");
}
