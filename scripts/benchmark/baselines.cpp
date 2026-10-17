// Compiled baseline integrators for the comparison with GenSAS (scripts/benchmark_solvers.py).
//
// RK45 (Dormand-Prince 5(4)) and DOP853 follow SciPy's implementation (same coefficients, error estimate and
// step size control), but run as compiled code. CVODE (BDF with a dense Newton iteration, from SUNDIALS) is
// the solver that many Modelica tools use for ODEs; it is built in when compiled with -DWITH_CVODE.
//
// Usage: baselines <method> <problem> <rtol> [n]       method: rk45, dop853, cvode
// Prints one line: steps rhs_evaluations seconds_per_run n final_state...
// The run is repeated until at least 0.3 s have passed and the time is the average of one run.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <string>
#include <vector>

#include "coefficients.h"

#ifdef WITH_CVODE
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#endif

typedef std::vector<double> Vec;
typedef std::function<void(double, const double *, double *)> Rhs;

struct Problem {
	int n = 0;
	double T = 1;
	Vec y0;
	Rhs f;
};

static Problem makeProblem(const std::string &name, int ringSize) {
	Problem p;
	if (name == "riccati") {
		p.n = 1;
		p.T = 0.99;
		p.y0 = {1.0};
		p.f = [](double, const double *y, double *d) { d[0] = y[0] * y[0]; };
	} else if (name == "runaway") {
		p.n = 1;
		p.T = 0.999;
		p.y0 = {0.0};
		p.f = [](double, const double *y, double *d) { d[0] = std::exp(y[0]); };
	} else if (name == "logistic") {
		p.n = 1;
		p.T = 40;
		p.y0 = {1e-4};
		p.f = [](double, const double *y, double *d) { d[0] = y[0] * (1 - y[0]); };
	} else if (name == "harmonic") {
		p.n = 2;
		p.T = 1000;
		p.y0 = {1.0, 0.0};
		p.f = [](double, const double *y, double *d) { d[0] = y[1]; d[1] = -y[0]; };
	} else if (name == "lorenz") {
		p.n = 3;
		p.T = 10;
		p.y0 = {1, 1, 1};
		p.f = [](double, const double *y, double *d) {
			d[0] = 10 * (y[1] - y[0]);
			d[1] = y[0] * (28 - y[2]) - y[1];
			d[2] = y[0] * y[1] - 8.0 / 3 * y[2];
		};
	} else if (name == "pendulum") {
		p.n = 2;
		p.T = 100;
		p.y0 = {3.0, 0.0};
		p.f = [](double, const double *y, double *d) { d[0] = y[1]; d[1] = -std::sin(y[0]); };
	} else if (name == "lotka") {
		p.n = 2;
		p.T = 100;
		p.y0 = {10, 10};
		p.f = [](double, const double *y, double *d) {
			d[0] = 1.1 * y[0] - 0.4 * y[0] * y[1];
			d[1] = 0.1 * y[0] * y[1] - 0.4 * y[1];
		};
	} else if (name == "kepler") {
		p.n = 4;
		p.T = 20 * 3.14159265358979323846;
		p.y0 = {1.9, 0, 0, 0.22941573387056177};
		p.f = [](double, const double *y, double *d) {
			const double r = std::sqrt(y[0] * y[0] + y[1] * y[1]);
			const double r3 = r * r * r;
			d[0] = y[2];
			d[1] = y[3];
			d[2] = -y[0] / r3;
			d[3] = -y[1] / r3;
		};
	} else if (name == "vdp1") {
		p.n = 2;
		p.T = 1000;
		p.y0 = {2, 0};
		p.f = [](double, const double *y, double *d) { d[0] = y[1]; d[1] = (1 - y[0] * y[0]) * y[1] - y[0]; };
	} else if (name == "vdp100") {
		p.n = 2;
		p.T = 200;
		p.y0 = {2, 0};
		p.f = [](double, const double *y, double *d) { d[0] = y[1]; d[1] = 100 * (1 - y[0] * y[0]) * y[1] - y[0]; };
	} else if (name == "robertson") {
		p.n = 2;
		p.T = 1000;
		p.y0 = {1, 0};
		p.f = [](double, const double *y, double *d) {
			const double y3 = 1 - y[0] - y[1];
			d[0] = -0.04 * y[0] + 1e4 * y[1] * y3;
			d[1] = 0.04 * y[0] - 1e4 * y[1] * y3 - 3e7 * y[1] * y[1];
		};
	} else if (name == "ring") {
		// a ring of coupled swing machines: states th[0..N-1], om[0..N-1]
		const int N = std::max(1, ringSize);
		std::vector<double> P(N);
		for (int i = 0; i < N; i++)
			P[i] = 0.5 * std::sin(1.7 * (i + 1));
		p.n = 2 * N;
		p.T = 20;
		p.y0 = Vec(2 * N, 0.0);
		p.f = [N, P](double, const double *y, double *d) {
			for (int i = 0; i < N; i++) {
				const double th = y[i];
				const double prev = y[(i + N - 1) % N];
				const double next = y[(i + 1) % N];
				d[i] = y[N + i];
				d[N + i] = P[i] - 0.2 * y[N + i] + (std::sin(prev - th) + std::sin(next - th));
			}
		};
	} else {
		std::fprintf(stderr, "unknown problem %s\n", name.c_str());
		std::exit(2);
	}
	return p;
}

static double rmsNorm(const Vec &x) {
	double s = 0;
	for (double v : x)
		s += v * v;
	return std::sqrt(s / x.size());
}

struct Stats {
	long steps = 0;
	long evals = 0;
};

// An explicit Runge-Kutta method with the step size control of scipy.integrate.RK45 / DOP853.
struct RkMethod {
	int stages;		 // without the FSAL stage
	int errorOrder;	 // order of the error estimate (4 and 7)
	const double *A; // stages x stages (row length `lda`)
	int lda;
	const double *B;
	const double *C;
	const double *E5;
	const double *E3; // null for RK45
};

static bool rkSolve(const RkMethod &m, const Problem &p, double rtol, double atol, Vec &y, Stats &st) {
	const int n = p.n;
	const int S = m.stages;
	std::vector<Vec> K(S + 1, Vec(n));
	Vec ytmp(n), ynew(n), scale(n), err5(n), err3(n);
	double t = 0;
	y = p.y0;
	auto rhs = [&](double tt, const Vec &yy, Vec &out) { st.evals++; p.f(tt, yy.data(), out.data()); };
	rhs(t, y, K[0]);
	// initial step (scipy's select_initial_step)
	{
		for (int i = 0; i < n; i++)
			scale[i] = atol + std::fabs(y[i]) * rtol;
		Vec a(n), b(n);
		for (int i = 0; i < n; i++) {
			a[i] = y[i] / scale[i];
			b[i] = K[0][i] / scale[i];
		}
		const double d0 = rmsNorm(a), d1 = rmsNorm(b);
		const double h0 = (d0 < 1e-5 || d1 < 1e-5) ? 1e-6 : 0.01 * d0 / d1;
		for (int i = 0; i < n; i++)
			ytmp[i] = y[i] + h0 * K[0][i];
		Vec f1(n);
		rhs(t + h0, ytmp, f1);
		Vec c(n);
		for (int i = 0; i < n; i++)
			c[i] = (f1[i] - K[0][i]) / scale[i];
		const double d2 = rmsNorm(c) / h0;
		double h1;
		if (d1 <= 1e-15 && d2 <= 1e-15)
			h1 = std::max(1e-6, h0 * 1e-3);
		else
			h1 = std::pow(0.01 / std::max(d1, d2), 1.0 / (m.errorOrder + 1));
		double habs = std::min({100 * h0, h1, p.T});
		const double exponent = -1.0 / (m.errorOrder + 1);
		bool rejected = false;
		while (t < p.T) {
			double hh = habs;
			double tnew = t + hh;
			if (tnew > p.T)
				tnew = p.T;
			hh = tnew - t;
			if (!(hh > 0.0) || tnew == t)
				return false; // the step size collapsed
			for (int s = 1; s < S; s++) {
				for (int i = 0; i < n; i++) {
					double acc = 0;
					for (int j = 0; j < s; j++)
						acc += K[j][i] * m.A[s * m.lda + j];
					ytmp[i] = y[i] + acc * hh;
				}
				rhs(t + m.C[s] * hh, ytmp, K[s]);
			}
			for (int i = 0; i < n; i++) {
				double acc = 0;
				for (int j = 0; j < S; j++)
					acc += K[j][i] * m.B[j];
				ynew[i] = y[i] + hh * acc;
			}
			rhs(tnew, ynew, K[S]);
			for (int i = 0; i < n; i++)
				scale[i] = atol + std::max(std::fabs(y[i]), std::fabs(ynew[i])) * rtol;
			double errNorm;
			if (m.E3 == nullptr) {
				double s2 = 0;
				for (int i = 0; i < n; i++) {
					double acc = 0;
					for (int j = 0; j <= S; j++)
						acc += K[j][i] * m.E5[j];
					const double e = hh * acc / scale[i];
					s2 += e * e;
				}
				errNorm = std::sqrt(s2 / n);
			} else {
				double n5 = 0, n3 = 0;
				for (int i = 0; i < n; i++) {
					double a5 = 0, a3 = 0;
					for (int j = 0; j <= S; j++) {
						a5 += K[j][i] * m.E5[j];
						a3 += K[j][i] * m.E3[j];
					}
					a5 /= scale[i];
					a3 /= scale[i];
					n5 += a5 * a5;
					n3 += a3 * a3;
				}
				errNorm = (n5 == 0 && n3 == 0) ? 0.0 : std::fabs(hh) * n5 / std::sqrt((n5 + 0.01 * n3) * n);
			}
			if (errNorm < 1) {
				double factor = errNorm == 0 ? 10.0 : std::min(10.0, 0.9 * std::pow(errNorm, exponent));
				if (rejected)
					factor = std::min(1.0, factor);
				habs *= factor;
				t = tnew;
				y = ynew;
				K[0] = K[S];
				st.steps++;
				rejected = false;
			} else {
				habs *= std::max(0.2, 0.9 * std::pow(errNorm, exponent));
				rejected = true;
			}
		}
	}
	return true;
}

#ifdef WITH_CVODE
struct CvCtx {
	const Problem *p;
	long evals = 0;
};
static int cvRhs(realtype t, N_Vector y, N_Vector yd, void *user) {
	CvCtx *c = static_cast<CvCtx *>(user);
	c->evals++;
	c->p->f(t, NV_DATA_S(y), NV_DATA_S(yd));
	return 0;
}
static bool cvodeSolve(const Problem &p, double rtol, double atol, Vec &y, Stats &st) {
	N_Vector yv = N_VNew_Serial(p.n);
	for (int i = 0; i < p.n; i++)
		NV_Ith_S(yv, i) = p.y0[i];
	void *mem = CVodeCreate(CV_BDF);
	CvCtx ctx;
	ctx.p = &p;
	CVodeInit(mem, cvRhs, 0.0, yv);
	CVodeSStolerances(mem, rtol, atol);
	CVodeSetUserData(mem, &ctx);
	CVodeSetMaxNumSteps(mem, 100000000);
	CVodeSetStopTime(mem, p.T);
	SUNMatrix A = SUNDenseMatrix(p.n, p.n);
	SUNLinearSolver LS = SUNLinSol_Dense(yv, A);
	CVodeSetLinearSolver(mem, LS, A);
	realtype t = 0;
	int flag = CVode(mem, p.T, yv, &t, CV_NORMAL);
	long steps = 0;
	CVodeGetNumSteps(mem, &steps);
	st.steps = steps;
	st.evals = ctx.evals;
	y.assign(NV_DATA_S(yv), NV_DATA_S(yv) + p.n);
	CVodeFree(&mem);
	SUNLinSolFree(LS);
	SUNMatDestroy(A);
	N_VDestroy(yv);
	return flag >= 0;
}
#endif

int main(int argc, char **argv) {
	if (argc < 4) {
		std::fprintf(stderr, "usage: baselines rk45|dop853|cvode problem rtol [n]\n");
		return 2;
	}
	const std::string method = argv[1];
	const double rtol = std::atof(argv[3]);
	const double atol = rtol * 1e-3;
	Problem p = makeProblem(argv[2], argc > 4 ? std::atoi(argv[4]) : 10);
	RkMethod rk45 = {6, 4, &coefficients::rk45_A[0][0], 5, coefficients::rk45_B, coefficients::rk45_C, coefficients::rk45_E, nullptr};
	RkMethod dop = {12, 7, &coefficients::dop853_A[0][0], 12, coefficients::dop853_B, coefficients::dop853_C, coefficients::dop853_E5, coefficients::dop853_E3};
	Vec y;
	Stats st;
	int runs = 0;
	const auto start = std::chrono::steady_clock::now();
	double elapsed = 0;
	do {
		st = Stats();
		if (method == "rk45" || method == "dop853") {
			if (!rkSolve(method == "rk45" ? rk45 : dop, p, rtol, atol, y, st)) {
				std::printf("failed\n");
				return 1;
			}
		}
#ifdef WITH_CVODE
		else if (method == "cvode") {
			if (!cvodeSolve(p, rtol, atol, y, st)) {
				std::printf("failed\n");
				return 1;
			}
		}
#endif
		else {
#ifndef WITH_CVODE
			if (method == "cvode") {
				std::fprintf(stderr, "cvode is not available: compile with -DWITH_CVODE and link SUNDIALS\n");
				return 2;
			}
#endif
			std::fprintf(stderr, "unknown method %s\n", method.c_str());
			return 2;
		}
		runs++;
		elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
	} while (elapsed < 0.3 && runs < 100000);
	std::printf("%ld %ld %.9g %d", st.steps, st.evals, elapsed / runs, p.n);
	for (double v : y)
		std::printf(" %.17g", v);
	std::printf("\n");
	return 0;
}
