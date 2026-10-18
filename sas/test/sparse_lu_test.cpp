// SparseLU: a sparse factorization that is solved for several right-hand sides.
#include "sas/SparseLU.h"
#include "sas/test/test_util.h"

#include <cmath>
#include <random>
#include <vector>

using che::core::SparseLU;

namespace {
	struct Csc {
		int n = 0;
		std::vector<int> colPtr, rowIdx;
		std::vector<double> val;
	};

	// A random sparse matrix with a nonzero diagonal that is not dominant: the pivoting has to work.
	Csc randomMatrix(int n, double density, unsigned seed) {
		std::mt19937 gen(seed);
		std::uniform_real_distribution<double> u(-1.0, 1.0), p(0.0, 1.0);
		Csc a;
		a.n = n;
		a.colPtr.push_back(0);
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				if (i == j) {
					a.rowIdx.push_back(i);
					a.val.push_back(0.3 * u(gen) + (u(gen) > 0 ? 0.5 : -0.5));
				} else if (p(gen) < density) {
					a.rowIdx.push_back(i);
					a.val.push_back(u(gen));
				}
			}
			a.colPtr.push_back((int)a.rowIdx.size());
		}
		return a;
	}

	std::vector<double> multiply(const Csc &a, const std::vector<double> &x) {
		std::vector<double> y(a.n, 0.0);
		for (int j = 0; j < a.n; j++) {
			for (int k = a.colPtr[j]; k < a.colPtr[j + 1]; k++) {
				y[a.rowIdx[k]] += a.val[k] * x[j];
			}
		}
		return y;
	}
} // namespace

static void testSolves() {
	for (int n : {1, 2, 7, 60, 300}) {
		Csc a = randomMatrix(n, n < 20 ? 0.5 : 0.04, 100 + n);
		SparseLU lu;
		CHECK(!lu.factored());
		CHECK(lu.factor(a.n, a.colPtr.data(), a.rowIdx.data(), a.val.data()));
		CHECK(lu.factored());
		std::mt19937 gen(5);
		std::uniform_real_distribution<double> u(-1.0, 1.0);
		// the same factorization for several right-hand sides
		for (int rhs = 0; rhs < 4; rhs++) {
			std::vector<double> x(n);
			for (double &v : x) {
				v = u(gen);
			}
			std::vector<double> b = multiply(a, x);
			CHECK(lu.solve(b.data()));
			double err = 0.0, scale = 0.0;
			for (int i = 0; i < n; i++) {
				err = std::fmax(err, std::fabs(b[i] - x[i]));
				scale = std::fmax(scale, std::fabs(x[i]));
			}
			CHECK(err <= 1e-8 * scale);
		}
	}
}

static void testRefactorAndSingular() {
	SparseLU lu;
	Csc a = randomMatrix(30, 0.1, 1);
	CHECK(lu.factor(a.n, a.colPtr.data(), a.rowIdx.data(), a.val.data()));
	// a matrix with an empty column is singular; the object is then not usable until it is factored again
	Csc s = randomMatrix(30, 0.1, 2);
	s.val.erase(s.val.begin() + s.colPtr[4], s.val.begin() + s.colPtr[5]);
	s.rowIdx.erase(s.rowIdx.begin() + s.colPtr[4], s.rowIdx.begin() + s.colPtr[5]);
	const int removed = s.colPtr[5] - s.colPtr[4];
	for (size_t j = 5; j < s.colPtr.size(); j++) {
		s.colPtr[j] -= removed;
	}
	CHECK(!lu.factor(s.n, s.colPtr.data(), s.rowIdx.data(), s.val.data()));
	CHECK(!lu.factored());
	std::vector<double> b(30, 1.0);
	CHECK(!lu.solve(b.data()));
	CHECK(lu.factor(a.n, a.colPtr.data(), a.rowIdx.data(), a.val.data()));
	CHECK(lu.solve(b.data()));

	// two equal columns: singular with no empty column or row
	Csc d;
	d.n = 3;
	d.colPtr = {0, 2, 4, 6};
	d.rowIdx = {0, 1, 0, 1, 1, 2};
	d.val = {1.0, 2.0, 1.0, 2.0, 1.0, 1.0};
	d.val[4] = 3.0;
	d.val[5] = 4.0;
	CHECK(lu.factor(d.n, d.colPtr.data(), d.rowIdx.data(), d.val.data()) || !lu.factored());
	// an empty row (no entry in row 2)
	Csc e;
	e.n = 3;
	e.colPtr = {0, 1, 2, 3};
	e.rowIdx = {0, 1, 1};
	e.val = {1.0, 1.0, 1.0};
	CHECK(!lu.factor(e.n, e.colPtr.data(), e.rowIdx.data(), e.val.data()));
}

int main() {
	testSolves();
	testRefactorAndSingular();
	return sastest::finish("sparse_lu_test");
}
