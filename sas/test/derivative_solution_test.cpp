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

// getDerValue is the analytic derivative of the evaluated solution (for a Pade approximant, the derivative of
// the approximant itself): it matches a central difference of getSolValue.
static void checkAnalyticDerivative(int type) {
	const int nState = 5, nLvl = 15;
	std::unique_ptr<CheSolution> sol(CheSolutionFactory::makeInitCheSol(type, nState, nLvl));
	fillSeries(sol.get(), nState, nLvl);
	const double h = 1e-5;
	for (double alpha : {0.0, 0.05, 0.3, 0.9}) {
		const arma::vec d = sol->getDerValue(alpha, nState);
		const arma::vec hi = sol->getSolValue(alpha + h);
		const arma::vec lo = sol->getSolValue(alpha - h);
		CHECK_EQ((int)d.n_rows, nState);
		for (int r = 0; r < nState && r < (int)d.n_rows; r++) {
			CHECK_NEAR(d(r), (hi(r) - lo(r)) / (2 * h), 1e-6 * (1.0 + std::fabs(d(r))));
		}
	}
	// a smaller number of rows gives the first rows
	const arma::vec all = sol->getDerValue(0.3, nState);
	const arma::vec head = sol->getDerValue(0.3, 3);
	CHECK_EQ((int)head.n_rows, 3);
	for (int r = 0; r < 3 && r < (int)head.n_rows; r++) {
		CHECK_EQ(all(r), head(r));
	}
}

// getSolValueRows gives the entries of getSolValue for the listed rows, in the order of the list (repeats allowed), also
// for a row of a Pade approximant that is not finite at alpha (it falls back to the power series, for good).
static void checkRows(int type) {
	const int nState = 6, nLvl = 15;
	auto make = [&]() {
		std::unique_ptr<CheSolution> s(CheSolutionFactory::makeInitCheSol(type, nState, nLvl));
		fillSeries(s.get(), nState, nLvl);
		for (int k = 0; k < nLvl; k++) {
			s->solution(4, k) = 1.0; // 1 / (1 - t): the Pade approximant has its pole at t = 1
		}
		return s;
	};
	std::unique_ptr<CheSolution> sol = make(), ref = make();
	const std::vector<int> rows = {5, 0, 4, 2, 2};
	// (alpha = 1 is the pole of row 4: after it the row stays on its power series, and so the values at 0.3 again)
	for (double alpha : {0.0, 0.3, 0.9, 1.0, 0.3, 1.0}) {
		const arma::vec some = sol->getSolValueRows(alpha, rows);
		const arma::vec all = ref->getSolValue(alpha);
		CHECK_EQ((int)some.n_rows, (int)rows.size());
		for (size_t i = 0; i < rows.size() && i < some.n_rows; i++) {
			// (not bitwise: a compiler may fuse the multiply and add of the Horner sums differently in the two loops)
			CHECK_NEAR(some(i), all(rows[i]), 1e-12 * (1.0 + std::fabs(all(rows[i]))));
		}
	}
	CHECK_EQ((int)sol->getSolValueRows(0.5, std::vector<int>()).n_rows, 0);
}

int main() {
	checkRows(CHESOL_PADE);
	checkRows(CHESOL_PS);
	check(CHESOL_PADE);
	check(CHESOL_PS);
	checkAnalyticDerivative(CHESOL_PADE);
	checkAnalyticDerivative(CHESOL_PS);
	return sastest::finish("derivative_solution_test");
}
