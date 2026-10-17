# How the SAS solver compares with conventional integrators

This page says where the semi-analytical solver (SAS) of GenSAS is better than the integrators in common use, where it is not, and how to check it yourself. The benefit of a solver can be accuracy, speed, or both; the measurements below separate them, because the answer is different for each. They come from two scripts: [scripts/benchmark_solvers.py](../scripts/benchmark_solvers.py) (steps, error and run time against compiled integrators) and [scripts/compare_with_scipy.py](../scripts/compare_with_scipy.py) (steps and error against the methods of SciPy).

## In short

* **Accuracy for a given number of steps: better.** On smooth nonlinear problems GenSAS takes about 2 to 3 times fewer steps than the eighth order `DOP853` (the same number on the network below), 5 to 20 times fewer than a fifth order Runge-Kutta method, and 20 to 70 times fewer than BDF, for the same error. A smooth transition such as logistic growth is solved to 4e-14 in 8 steps.
* **Run time against the BDF solver of Modelica tools (CVODE): comparable, and better on some problems.** At errors of 1e-6 to 1e-8 on smooth problems GenSAS runs from 7 times slower (a ring of 10 machines) to 8 times faster (a blow-up) than CVODE, and it reaches errors that CVODE does not reach at any tolerance. On the stiff problems tested it is 10 to 20000 times slower.
* **Run time against an explicit Runge-Kutta method: much slower.** A compiled `DOP853` is about 50 to 200 times faster than GenSAS on every non-stiff problem here. For such a problem an explicit method is the right tool; the comparison matters for models where it cannot be used.
* **The highest accuracies are out of reach.** `errorTol` is limited to 1e-10 and the error of a long run grows beyond it: the Kepler orbit ends with an error of 2e-8, where `DOP853` gets 1e-10.
* **Not measured here: large sparse DAEs**, the case that GenSAS was built for. The ring of machines below is a first step in that direction.

The run times are those of the current main branch, which includes the speed-up described under [Run time](#run-time).

## What is different

A conventional integrator takes a step by combining a few evaluations of the right-hand side (Runge-Kutta), or by solving a nonlinear equation with Newton's method for the new state (BDF, Radau, the trapezoidal rule). SAS builds the solution of the step as a power series in time, one coefficient after another, each from a linear system with the same sparse matrix at every order, and turns the series into a Pade approximant. The result is:

* **A solution that is a function, not a list of points.** The approximant can be evaluated at any time in the step. Output samples, the check of equation residuals, and the search for event times all use this function; none of them needs extra steps or an interpolation of lower order than the solver.
* **No Newton iteration in time stepping.** Algebraic equations of a DAE are part of the series recurrence, so a step solves a sequence of linear systems with the same matrix, not a nonlinear iteration that may fail to converge.
* **A step size that follows the analytic structure of the solution.** The Pade approximant is useful beyond the radius of convergence of the series, and the step is chosen from the residual of the equations, not from a stability limit of the method.
* **A high and adjustable order** (default 15, up to 1000) instead of the fixed order of a Runge-Kutta pair.

The price is a more expensive step: one linear solve per order plus the series arithmetic of the nonlinear terms, and a model with many nonlinear functions has more series to carry.

## What was measured

Every problem is solved with GenSAS and with compiled C++ integrators whose right-hand side is compiled code, so that GenSAS is not compared with Python call overhead:

* `RK45`, the Dormand-Prince 5(4) pair, and `DOP853`, the eighth order method, both with the coefficients, error estimate and step size control of SciPy (the step counts agree with SciPy's exactly);
* `CVODE`, the BDF method with a dense Newton iteration from SUNDIALS, which is what many Modelica tools use for ODEs.

The methods of SciPy (`LSODA`, `BDF`, `Radau` as well) are used for step counts in the appendix. The reference is the exact solution where there is one (the logistic, blow-up, runaway, oscillator and Kepler problems), and otherwise a Taylor series solution at 30 digits computed with `mpmath` (kept in [scripts/benchmark/references.json](../scripts/benchmark/references.json), made by [scripts/benchmark/reference.py](../scripts/benchmark/reference.py)). Where that is too slow (the rings of 100 and 200 machines, Robertson) the reference is the SciPy solution with a tolerance of 1e-13, and its uncertainty is taken as the largest difference to two other double precision runs; for every problem where both exist this estimate was larger than the true error of that reference, by a factor of 2 to 10. The scripts print the uncertainty and mark with a dagger (†) any error that is within 10 times of it, because such an error is not resolved by the reference. An earlier version of the tables used the SciPy solution as the reference for all problems without a closed form. It was less accurate than its tolerance suggests (3e-10 for the pendulum, 6e-12 for Lotka-Volterra, 4e-13 for Lorenz), so entries of those problems below about 1e-9 in the tables that follow are not reliable, and the statement that `DOP853` reaches 1e-10 on the pendulum is wrong: against the accurate reference no method reaches 1e-10 there, because the pendulum starts near the unstable equilibrium and amplifies every error. The tables will be regenerated with the new references. The error is the largest error of the final state, scaled by `max(|reference|, 1)`.

Read the tables with these points in mind:

* **The tolerances are not the same quantity.** `errorTol` of GenSAS bounds the residual of the equations in a step; `rtol` of the others bounds an estimate of the local error. So the main table asks the question the other way round: **what does each method spend to reach a given error**, trying a ladder of tolerances from 1e-4 down to 1e-13 (1e-6 to 1e-10 for GenSAS, the range it accepts).
* **Times** are one run in milliseconds. For GenSAS they are the `Computation time` of the application, best of three, which includes building the equations from the model text. The baselines compute the final state only; GenSAS also writes a few samples. Compiled right-hand sides are an advantage that GenSAS, which reads the model as text, does not have.
* **Numbers are from one machine and one version** (GenSAS from the main branch of October 2026, SciPy 1.15, SUNDIALS 5.8). The pattern matters, not the third digit.

## Accuracy: examples that show the strengths

* **A smooth transition.** Logistic growth `x' = x (1 - x)` from 1e-4 to 1 is solved to an error of 4e-14 in 8 steps, at any `errorTol` from 1e-6 to 1e-10. `DOP853` needs 12 to 25 steps for the errors of 1e-4 to 1e-10, `RK45` 26 to 315, CVODE 84 to 574. The Pade approximant represents the sigmoid almost exactly.
* **A long oscillation.** The linear oscillator over 1000 s (160 periods) needs 732 steps for an error of 1e-6 where `DOP853` needs 1756, `RK45` 12523 and CVODE 55908.
* **A chaotic system.** Lorenz to t = 10: 190 steps for 1e-8 (`DOP853` 351, `RK45` 3195, CVODE 9837), and a final error of 9e-12 at `errorTol` 1e-10, which CVODE does not reach at all.
* **A solution that runs into a singularity.** `x' = x^2` up to t = 0.99, where the solution is 100 and infinite at t = 1: 2 to 4 steps for errors of 1e-4 to 1e-8, against 11 to 30 for `DOP853`, 32 to 122 for `RK45` and 238 to 1566 for CVODE. The thermal runaway `u' = exp(u)` is the same kind of problem, and there the advantage is small: 20 steps for GenSAS and 18 for `DOP853` at 1e-6.
* **A network.** In a ring of coupled swing machines the number of steps does not depend on the size: 24 to 25 steps for an error of 1e-6 and 33 to 34 for 1e-8 from 10 to 200 machines (20 to 400 states), about as many as `DOP853`, a fifth to a tenth of `RK45` and a twentieth to a fiftieth of CVODE.
* **Event times** are located on the continuous solution. For a pendulum that passes through zero six times in 100 s, GenSAS finds the crossings with an error of 8e-8 s in 193 steps; `DOP853` needs 386 steps for 3e-7 s, `RK45` 2682 steps for 4e-7 s, `LSODA` 2382 steps for 3e-6 s.
* **The output grid is free.** The Lorenz run takes the same 190 steps whether the result has 10, 1000 or 10000 samples; the time goes from 19 ms to 27 ms (before the speed-up below).

## Run time

Fewer steps do not yet mean less time. A step of GenSAS costs about 20 to 100 microseconds for a few states, where a step of `DOP853` costs well below one microsecond. At the same error:

* **Against CVODE** (the implicit method of the Modelica tools), the time of GenSAS is within a small factor, on either side. Lorenz: 4.4 ms against 2.8 ms for an error of 1e-6, 6.1 ms against 6.8 ms for 1e-8. Pendulum, 1e-6: 4.6 ms against 5.9 ms. Linear oscillator, 1e-6: 11 ms against 31 ms. The blow-up: 3 to 8 times faster. The Kepler orbit and the Van der Pol oscillator with mu = 1 take about twice as long, a ring of 10 machines five to seven times as long.
* **Against an explicit method**, GenSAS is 50 to 200 times slower on these non-stiff problems (Lorenz, 1e-8: 6.1 ms against 0.06 ms for `DOP853`).
* **On stiff problems GenSAS is not competitive**: the Van der Pol oscillator with mu = 100 takes 37 ms against 1 ms for CVODE, and the Robertson problem 3.6 s against 0.2 ms.
* **With the size of a network GenSAS gains on CVODE.** The ring of machines at an error of 1e-8: 158 ms for 100 machines and 351 ms for 200 (time grows like the size), against 92 ms and 495 ms for CVODE (like the size to the power 2.4). At 200 machines GenSAS is faster, and the lines cross near 150. The CVODE here uses a dense Jacobian, as a general purpose code does by default; with a sparse linear solver it would scale better, and that has not been measured. An explicit method is far ahead on this non-stiff model (1.9 ms for 200 machines), but a power grid is not non-stiff.

### Why a step is expensive, and what was done about it

A profile of a three-state model showed that most of the time of a step was not mathematics: 47% went into temporary matrices in the Levinson recursion that computes the Pade coefficients, and 27% into the set-up of SuperLU for a matrix of four unknowns, once per order. Replacing both (plain loops for the Pade routines, a dense LU, factored once per step, for algebraic systems of up to 64 unknowns) made the runs 1.4 to 3.3 times faster with the same results (Lorenz from 18 ms to 6 ms). What is left, in the order of the profile:

* the evaluation of the solution at the many trial times of the step size search (about 60% of the time, mostly allocations in the Pade evaluation),
* SuperLU, called once per order for algebraic systems of more than 64 unknowns, where a single factorization per step would do,
* the interpretation of the model: the series arithmetic walks the expression tree of the equations.

A further factor of several seems possible without changing the method; that is where the speed of GenSAS on small and medium models will come from.

## Where SAS is not the better choice

* **Very stiff problems with a fast transient.** The Robertson problem (a stiff chemical reaction written as a DAE) needs 80920 steps for GenSAS at 1e-6, against 192 to 287 for BDF methods, and GenSAS cannot start it at 1e-10. The cause has not been analyzed; the steps are short from the first one. Explicit methods are worse (`DOP853` needs 880000 steps), so GenSAS is in between, but not a good choice for this kind of model.
* **Stiffness between jumps is fine in steps, not in time.** The Van der Pol oscillator with mu = 100 needs 1165 steps at 1e-6, close to CVODE (1117 to 1533) and a fifth of `DOP853`, but each step costs more.
* **The tightest accuracies.** `errorTol` stops at 1e-10, and the error of a long run grows beyond it (the Kepler orbit ends at 2e-8, the pendulum at 5e-8, the oscillator at 3e-8). `DOP853` reaches 1e-10 on all of them.
* **Models whose solution is not smooth.** At each switching point the solver stops and restarts, and a model that switches thousands of times costs thousands of restarts. A Zeno model is reported as an error.
* **`piecewise = smooth` and `abs`.** The smoothed form of `abs` can leave the domain of `sqrt` while the series are evaluated (the Chua circuit stops at t = 0.12 in this mode). The exact mode, which uses events, works.

## What would show more

The problems here are small. The case that GenSAS was built for is **large sparse DAEs, in particular power systems**, where the absence of a Newton iteration in the time stepping matters most and the publications listed in the [README](../README.md) report the results. A comparable measurement in this repository needs a network model with algebraic equations (generators behind reactances, a network of buses) and a DAE solver to compare with (IDA, or a BDF code with a sparse linear solver); see [Roadmap](Roadmap.md#large-systems). The ring above is an ODE and shows the trend of the cost with the size only.

## Reproducing the numbers

Build the application (see [Installation](Installation.md)), install `numpy` and `scipy`, and run from the repository root:

```bash
python3 scripts/benchmark_solvers.py --app bazel-bin/app/app --sundials /usr --targets 1e-6 1e-8
python3 scripts/benchmark_solvers.py --app bazel-bin/app/app --problems lorenz ring --ring-sizes 10 100
python3 scripts/compare_with_scipy.py --app bazel-bin/app/app --problems lorenz kepler events
```

`benchmark_solvers.py` compiles `scripts/benchmark/baselines.cpp` with the C++ compiler (`CXX`). `--sundials` is the installation of SUNDIALS (the directory that has `include/` and `lib/`, or a prefix such as `/usr`) and adds CVODE; without it only the two Runge-Kutta methods are used. Without `--targets` the table shows fixed tolerances, with `--targets` the cost of reaching each error. Both scripts run GenSAS through `app -g`, count the lines `Stage start` of the verbose output for the steps, and read the final state from the CSV file.

## Results

### What each method spends to reach an error

Steps and milliseconds of the cheapest tolerance that reaches the error ("not reached": no tolerance of the method did).

#### Finite-time blow-up: x' = x^2, x(0) = 1, up to t = 0.99 (the solution is 1/(1-t), infinite at t = 1)

| error at most | GenSAS steps, ms | RK45 steps, ms | DOP853 steps, ms | CVODE (BDF) steps, ms |
|---|---:|---:|---:|---:|
| 0.0001 | 2, 0.097 | 32, 0.004 | 11, 0.002 | 238, 0.131 |
| 1e-06 | 4, 0.121 | 75, 0.006 | 22, 0.005 | 739, 0.373 |
| 1e-08 | 3, 0.099 | 122, 0.009 | 30, 0.006 | 1566, 0.771 |
| 1e-10 | not reached | 501, 0.036 | 68, 0.008 | not reached |

#### Logistic growth x' = x (1 - x), x(0) = 1e-4, t = 0 to 40 (a smooth transition from 1e-4 to 1)

| error at most | GenSAS steps, ms | RK45 steps, ms | DOP853 steps, ms | CVODE (BDF) steps, ms |
|---|---:|---:|---:|---:|
| 0.0001 | 8, 0.223 | 26, 0.002 | 12, 0.002 | 84, 0.043 |
| 1e-06 | 8, 0.242 | 56, 0.005 | 19, 0.003 | 84, 0.043 |
| 1e-08 | 8, 0.216 | 130, 0.010 | 19, 0.003 | 202, 0.114 |
| 1e-10 | 8, 0.226 | 315, 0.024 | 25, 0.004 | 574, 0.294 |

#### Lorenz system (chaotic), t = 0 to 10

| error at most | GenSAS steps, ms | RK45 steps, ms | DOP853 steps, ms | CVODE (BDF) steps, ms |
|---|---:|---:|---:|---:|
| 0.0001 | 134, 4.1 | 512, 0.044 | 116, 0.025 | 2116, 1.4 |
| 1e-06 | 134, 4.4 | 1275, 0.106 | 203, 0.042 | 4415, 2.8 |
| 1e-08 | 190, 6.1 | 3195, 0.275 | 351, 0.063 | 9837, 6.8 |
| 1e-10 | 229, 6.8 | 8019, 0.662 | 616, 0.102 | not reached |

#### Kepler orbit with eccentricity 0.9, ten periods (the exact end state is the start state)

| error at most | GenSAS steps, ms | RK45 steps, ms | DOP853 steps, ms | CVODE (BDF) steps, ms |
|---|---:|---:|---:|---:|
| 0.0001 | 151, 12.3 | 951, 0.129 | 473, 0.134 | 5827, 4.6 |
| 1e-06 | 259, 21.6 | 1482, 0.174 | 615, 0.179 | 18151, 13.8 |
| 1e-08 | not reached | 5938, 0.631 | 1044, 0.259 | 39951, 29.5 |
| 1e-10 | not reached | 14982, 1.6 | 1818, 0.383 | not reached |

#### Pendulum released at 3.0 rad (near the unstable position), t = 0 to 100

| error at most | GenSAS steps, ms | RK45 steps, ms | DOP853 steps, ms | CVODE (BDF) steps, ms |
|---|---:|---:|---:|---:|
| 0.0001 | 130, 4.2 | 1063, 0.125 | 230, 0.061 | 5050, 2.8 |
| 1e-06 | 154, 4.6 | 2682, 0.296 | 386, 0.092 | 10846, 5.9 |
| 1e-08 | not reached | 6770, 0.732 | 668, 0.145 | not reached |
| 1e-10 | not reached | not reached | 889, 0.198 | not reached |

#### Linear oscillator x'' = -x over 1000 s (about 160 periods)

| error at most | GenSAS steps, ms | RK45 steps, ms | DOP853 steps, ms | CVODE (BDF) steps, ms |
|---|---:|---:|---:|---:|
| 0.0001 | 532, 8.3 | 7784, 0.613 | 990, 0.118 | 22027, 12.4 |
| 1e-06 | 732, 11.2 | 12523, 0.945 | 1756, 0.212 | 55908, 30.9 |
| 1e-08 | not reached | 31960, 2.3 | 3126, 0.388 | 106135, 59.0 |
| 1e-10 | not reached | 81063, 5.8 | 5567, 0.673 | not reached |

#### Van der Pol oscillator, mu = 1 (a limit cycle, not stiff), t = 0 to 1000

| error at most | GenSAS steps, ms | RK45 steps, ms | DOP853 steps, ms | CVODE (BDF) steps, ms |
|---|---:|---:|---:|---:|
| 0.0001 | 2210, 68.1 | 13444, 1.3 | 2704, 0.503 | 28939, 17.9 |
| 1e-06 | 3153, 97.0 | 21019, 1.9 | 4507, 0.764 | 91216, 51.7 |
| 1e-08 | 4598, 144 | 52092, 4.1 | 7807, 1.4 | 284994, 158 |
| 1e-10 | not reached | 131072, 10.2 | 13361, 2.1 | not reached |

#### Van der Pol oscillator, mu = 100 (stiff between the jumps), t = 0 to 200

| error at most | GenSAS steps, ms | RK45 steps, ms | DOP853 steps, ms | CVODE (BDF) steps, ms |
|---|---:|---:|---:|---:|
| 0.0001 | 1165, 38.0 | 11639, 1.1 | 6023, 0.814 | 792, 0.475 |
| 1e-06 | 1165, 37.1 | 11892, 1.1 | 6087, 0.801 | 1533, 0.960 |
| 1e-08 | 1240, 40.7 | 13781, 1.2 | 6199, 0.869 | 4438, 2.6 |
| 1e-10 | 1504, 49.2 | 25982, 2.1 | 7511, 1.2 | 8667, 4.9 |

#### A ring of 10 coupled swing machines (20 states), t = 0 to 20

| error at most | GenSAS steps, ms | RK45 steps, ms | DOP853 steps, ms | CVODE (BDF) steps, ms |
|---|---:|---:|---:|---:|
| 1e-06 | 24, 11.6 | 152, 0.127 | 37, 0.066 | 571, 1.6 |
| 1e-08 | 33, 16.2 | 390, 0.304 | 63, 0.110 | 1200, 3.3 |

#### A ring of 40 coupled swing machines (80 states), t = 0 to 20

| error at most | GenSAS steps, ms | RK45 steps, ms | DOP853 steps, ms | CVODE (BDF) steps, ms |
|---|---:|---:|---:|---:|
| 1e-06 | 25, 56.8 | 137, 0.437 | 28, 0.219 | 485, 9.4 |
| 1e-08 | 34, 75.4 | 354, 1.1 | 57, 0.388 | 1085, 20.8 |

#### A ring of 100 coupled swing machines (200 states), t = 0 to 20

| error at most | GenSAS steps, ms | RK45 steps, ms | DOP853 steps, ms | CVODE (BDF) steps, ms |
|---|---:|---:|---:|---:|
| 1e-06 | 24, 117 | 129, 1.1 | 27, 0.535 | 478, 45.7 |
| 1e-08 | 33, 158 | 327, 2.4 | 43, 0.752 | 1015, 92.3 |

#### A ring of 200 coupled swing machines (400 states), t = 0 to 20

| error at most | GenSAS steps, ms | RK45 steps, ms | DOP853 steps, ms | CVODE (BDF) steps, ms |
|---|---:|---:|---:|---:|
| 1e-06 | 25, 263 | 128, 2.1 | 34, 1.3 | 478, 165 |
| 1e-08 | 34, 351 | 324, 4.7 | 54, 1.9 | 1532, 495 |

### Fixed tolerances

#### Finite-time blow-up: x' = x^2, x(0) = 1, up to t = 0.99 (the solution is 1/(1-t), infinite at t = 1)

| method | steps | error | time (ms) |
|---|---:|---:|---:|
| GenSAS, errorTol 1e-06 | 2 | 2.9e-05 | 0.098 |
| GenSAS, errorTol 1e-10 | 3 | 3.9e-09 | 0.099 |
| RK45, rtol 1e-06 | 32 | 2.8e-05 | 0.004 |
| RK45, rtol 1e-10 | 197 | 1.6e-09 | 0.014 |
| DOP853, rtol 1e-06 | 13 | 8.9e-06 | 0.003 |
| DOP853, rtol 1e-10 | 40 | 8.1e-10 | 0.008 |
| CVODE (BDF), rtol 1e-06 | 111 | 1.8e-03 | 0.068 |
| CVODE (BDF), rtol 1e-10 | 498 | 1.0e-06 | 0.255 |

#### Thermal runaway: u' = exp(u), u(0) = 0, up to t = 0.999 (the solution is -ln(1-t), infinite at t = 1)

| method | steps | error | time (ms) |
|---|---:|---:|---:|
| GenSAS, errorTol 1e-06 | 20 | 7.0e-06 | 0.477 |
| GenSAS, errorTol 1e-10 | 26 | 1.8e-09 | 0.511 |
| RK45, rtol 1e-06 | 27 | 5.2e-05 | 0.005 |
| RK45, rtol 1e-10 | 132 | 3.4e-09 | 0.015 |
| DOP853, rtol 1e-06 | 18 | 1.6e-06 | 0.006 |
| DOP853, rtol 1e-10 | 42 | 2.2e-10 | 0.014 |
| CVODE (BDF), rtol 1e-06 | 110 | 1.7e-03 | 0.073 |
| CVODE (BDF), rtol 1e-10 | 490 | 1.0e-06 | 0.260 |

#### Kepler orbit with eccentricity 0.9, ten periods (the exact end state is the start state)

| method | steps | error | time (ms) |
|---|---:|---:|---:|
| GenSAS, errorTol 1e-06 | 151 | 2.7e-05 | 12.4 |
| GenSAS, errorTol 1e-10 | 310 | 1.7e-08 | 24.9 |
| RK45, rtol 1e-06 | 607 | 6.7e-04 | 0.087 |
| RK45, rtol 1e-10 | 3732 | 5.9e-08 | 0.413 |
| DOP853, rtol 1e-06 | 290 | 1.1e-03 | 0.088 |
| DOP853, rtol 1e-10 | 806 | 6.0e-08 | 0.229 |
| CVODE (BDF), rtol 1e-06 | 2674 | 6.9e-03 | 2.3 |
| CVODE (BDF), rtol 1e-10 | 12613 | 2.2e-06 | 9.7 |

#### Logistic growth x' = x (1 - x), x(0) = 1e-4, t = 0 to 40 (a smooth transition from 1e-4 to 1)

| method | steps | error | time (ms) |
|---|---:|---:|---:|
| GenSAS, errorTol 1e-06 | 8 | 4.5e-14 | 0.227 |
| GenSAS, errorTol 1e-10 | 14 | 4.8e-14 | 0.328 |
| RK45, rtol 1e-06 | 56 | 1.6e-07 | 0.005 |
| RK45, rtol 1e-10 | 315 | 2.9e-12 | 0.024 |
| DOP853, rtol 1e-06 | 19 | 1.1e-09 | 0.003 |
| DOP853, rtol 1e-10 | 55 | 3.7e-14 | 0.007 |
| CVODE (BDF), rtol 1e-06 | 156 | 7.3e-07 | 0.084 |
| CVODE (BDF), rtol 1e-10 | 574 | 1.0e-11 | 0.294 |

#### Linear oscillator x'' = -x over 1000 s (about 160 periods)

| method | steps | error | time (ms) |
|---|---:|---:|---:|
| GenSAS, errorTol 1e-06 | 451 | 1.3e-04 | 6.9 |
| GenSAS, errorTol 1e-10 | 854 | 2.9e-08 | 12.9 |
| RK45, rtol 1e-06 | 4866 | 1.1e-04 | 0.422 |
| RK45, rtol 1e-10 | 31960 | 9.4e-09 | 2.3 |
| DOP853, rtol 1e-06 | 990 | 6.3e-05 | 0.121 |
| DOP853, rtol 1e-10 | 3126 | 5.5e-09 | 0.378 |
| CVODE (BDF), rtol 1e-06 | 10818 | 1.2e-03 | 6.1 |
| CVODE (BDF), rtol 1e-10 | 55908 | 2.5e-07 | 31.0 |

#### Van der Pol oscillator, mu = 1 (a limit cycle, not stiff), t = 0 to 1000

| method | steps | error | time (ms) |
|---|---:|---:|---:|
| GenSAS, errorTol 1e-06 | 2210 | 3.4e-06 | 67.9 |
| GenSAS, errorTol 1e-10 | 4598 | 8.1e-10 | 144 |
| RK45, rtol 1e-06 | 8668 | 1.4e-04 | 0.898 |
| RK45, rtol 1e-10 | 52092 | 7.4e-09 | 4.1 |
| DOP853, rtol 1e-06 | 3306 | 1.6e-04 | 0.558 |
| DOP853, rtol 1e-10 | 10209 | 6.2e-10 | 1.7 |
| CVODE (BDF), rtol 1e-06 | 28939 | 3.0e-05 | 17.8 |
| CVODE (BDF), rtol 1e-10 | 131928 | 1.6e-08 | 74.2 |

#### Lorenz system (chaotic), t = 0 to 10

| method | steps | error | time (ms) |
|---|---:|---:|---:|
| GenSAS, errorTol 1e-06 | 134 | 6.2e-08 | 4.2 |
| GenSAS, errorTol 1e-10 | 272 | 9.3e-12 | 8.3 |
| RK45, rtol 1e-06 | 327 | 2.4e-04 | 0.030 |
| RK45, rtol 1e-10 | 2018 | 3.0e-08 | 0.171 |
| DOP853, rtol 1e-06 | 116 | 8.4e-05 | 0.025 |
| DOP853, rtol 1e-10 | 351 | 9.8e-10 | 0.064 |
| CVODE (BDF), rtol 1e-06 | 991 | 1.2e-03 | 0.651 |
| CVODE (BDF), rtol 1e-10 | 4415 | 7.5e-07 | 2.8 |

#### Pendulum released at 3.0 rad (near the unstable position), t = 0 to 100

| method | steps | error | time (ms) |
|---|---:|---:|---:|
| GenSAS, errorTol 1e-06 | 88 | 1.8e-03 | 3.0 |
| GenSAS, errorTol 1e-10 | 186 | 5.3e-08 | 5.5 |
| RK45, rtol 1e-06 | 422 | 3.0e-03 | 0.054 |
| RK45, rtol 1e-10 | 2682 | 4.4e-07 | 0.306 |
| DOP853, rtol 1e-06 | 136 | 8.9e-04 | 0.039 |
| DOP853, rtol 1e-10 | 386 | 3.3e-07 | 0.094 |
| CVODE (BDF), rtol 1e-06 | 1224 | 1.3e-01 | 0.766 |
| CVODE (BDF), rtol 1e-10 | 5050 | 2.1e-05 | 2.9 |

#### Predator and prey (Lotka-Volterra), t = 0 to 100

| method | steps | error | time (ms) |
|---|---:|---:|---:|
| GenSAS, errorTol 1e-06 | 91 | 8.1e-04 | 2.8 |
| GenSAS, errorTol 1e-10 | 192 | 1.0e-08 | 5.7 |
| RK45, rtol 1e-06 | 453 | 3.4e-04 | 0.043 |
| RK45, rtol 1e-10 | 2801 | 3.9e-09 | 0.231 |
| DOP853, rtol 1e-06 | 159 | 1.5e-06 | 0.029 |
| DOP853, rtol 1e-10 | 491 | 5.9e-09 | 0.077 |
| CVODE (BDF), rtol 1e-06 | 1488 | 7.1e-03 | 0.883 |
| CVODE (BDF), rtol 1e-10 | 6831 | 3.1e-06 | 3.8 |

#### Van der Pol oscillator, mu = 100 (stiff between the jumps), t = 0 to 200

| method | steps | error | time (ms) |
|---|---:|---:|---:|
| GenSAS, errorTol 1e-06 | 1165 | 8.2e-08 | 37.7 |
| GenSAS, errorTol 1e-10 | 1504 | 1.2e-11 | 48.1 |
| RK45, rtol 1e-06 | 11892 | 2.7e-07 | 1.1 |
| RK45, rtol 1e-10 | 25982 | 2.6e-11 | 2.0 |
| DOP853, rtol 1e-06 | 6087 | 4.3e-08 | 0.802 |
| DOP853, rtol 1e-10 | 9534 | 4.8e-12 | 1.2 |
| CVODE (BDF), rtol 1e-06 | 1117 | 1.2e-05 | 0.714 |
| CVODE (BDF), rtol 1e-10 | 4438 | 2.2e-09 | 2.6 |

#### A ring of 10 coupled swing machines (20 states), t = 0 to 20

| method | steps | error | time (ms) |
|---|---:|---:|---:|
| GenSAS, errorTol 1e-06 | 24 | 1.5e-07 | 11.6 |
| GenSAS, errorTol 1e-10 | 48 | 2.8e-11 | 23.4 |
| RK45, rtol 1e-06 | 152 | 1.8e-07 | 0.126 |
| RK45, rtol 1e-10 | 1001 | 1.3e-11 | 0.745 |
| DOP853, rtol 1e-06 | 37 | 1.3e-07 | 0.065 |
| DOP853, rtol 1e-10 | 110 | 8.5e-12 | 0.186 |
| CVODE (BDF), rtol 1e-06 | 397 | 1.1e-06 | 1.1 |
| CVODE (BDF), rtol 1e-10 | 1771 | 4.4e-10 | 4.8 |

#### A ring of 40 coupled swing machines (80 states), t = 0 to 20

| method | steps | error | time (ms) |
|---|---:|---:|---:|
| GenSAS, errorTol 1e-06 | 25 | 2.3e-07 | 55.7 |
| GenSAS, errorTol 1e-10 | 47 | 4.3e-11 | 103 |
| RK45, rtol 1e-06 | 137 | 3.4e-07 | 0.432 |
| RK45, rtol 1e-10 | 912 | 2.6e-11 | 2.7 |
| DOP853, rtol 1e-06 | 34 | 1.1e-07 | 0.237 |
| DOP853, rtol 1e-10 | 100 | 9.6e-12 | 0.672 |
| CVODE (BDF), rtol 1e-06 | 347 | 2.0e-06 | 7.2 |
| CVODE (BDF), rtol 1e-10 | 1591 | 9.3e-10 | 29.4 |

#### A ring of 100 coupled swing machines (200 states), t = 0 to 20

| method | steps | error | time (ms) |
|---|---:|---:|---:|
| GenSAS, errorTol 1e-06 | 24 | 8.0e-07 | 116 |
| GenSAS, errorTol 1e-10 | 46 | 4.5e-11 | 226 |
| RK45, rtol 1e-06 | 129 | 3.3e-07 | 1.0 |
| RK45, rtol 1e-10 | 836 | 3.1e-11 | 6.0 |
| DOP853, rtol 1e-06 | 34 | 6.8e-08 | 0.630 |
| DOP853, rtol 1e-10 | 95 | 8.1e-12 | 1.6 |
| CVODE (BDF), rtol 1e-06 | 340 | 1.5e-06 | 33.7 |
| CVODE (BDF), rtol 1e-10 | 1426 | 1.3e-09 | 130 |

#### Robertson chemical kinetics (very stiff), written as a DAE with y3 = 1 - y1 - y2, t = 0 to 1000

| method | steps | error | time (ms) |
|---|---:|---:|---:|
| GenSAS, errorTol 1e-06 | 80920 | 2.0e-06 | 3508 |
| GenSAS, errorTol 1e-10 | 2 | failed: Error: the solver cannot advance beyond t = 0 within the error tolerance; the solution may be singular here | |
| RK45, rtol 1e-06 | 1707165 | 8.4e-08 | 169 |
| RK45, rtol 1e-10 | 1707617 | 4.8e-14 | 170 |
| DOP853, rtol 1e-06 | 883980 | 2.4e-08 | 126 |
| DOP853, rtol 1e-10 | 884072 | 1.8e-14 | 127 |
| CVODE (BDF), rtol 1e-06 | 287 | 4.2e-06 | 0.173 |
| CVODE (BDF), rtol 1e-10 | 1078 | 1.7e-10 | 0.719 |

### Event times and the methods of SciPy

#### Event times: the pendulum above passes th = 0 (downwards) 6 times in 100 s

| method | steps | largest error of an event time (s) |
|---|---:|---:|
| GenSAS, errorTol 1e-06 | 97 | 3.8e-03 |
| GenSAS, errorTol 1e-10 | 193 | 8.5e-08 |
| RK45, rtol 1e-06 | 422 | 2.8e-03 |
| RK45, rtol 1e-10 | 2682 | 4.0e-07 |
| DOP853, rtol 1e-06 | 136 | 9.3e-04 |
| DOP853, rtol 1e-10 | 386 | 3.1e-07 |
| LSODA, rtol 1e-06 | 1132 | 8.7e-03 |
| LSODA, rtol 1e-10 | 2382 | 3.3e-06 |
| Radau, rtol 1e-06 | 990 | 9.2e-05 |
| Radau, rtol 1e-10 | 9639 | 7.2e-10 |

### Appendix: steps of the implicit methods of SciPy

Steps and error for `LSODA`, `BDF` and `Radau` as well, from `compare_with_scipy.py` (steps only, no times).

#### Finite-time blow-up: x' = x^2, x(0) = 1, up to t = 0.99 (the solution is 1/(1-t), infinite at t = 1)

| method | steps | error |
|---|---:|---:|
| GenSAS, errorTol 1e-06 | 2 | 2.9e-05 |
| GenSAS, errorTol 1e-10 | 3 | 3.9e-09 |
| RK45, rtol 1e-06 | 32 | 2.8e-05 |
| RK45, rtol 1e-10 | 197 | 1.6e-09 |
| DOP853, rtol 1e-06 | 13 | 8.9e-06 |
| DOP853, rtol 1e-10 | 40 | 8.1e-10 |
| LSODA, rtol 1e-06 | 79 | 8.9e-04 |
| LSODA, rtol 1e-10 | 201 | 1.9e-07 |
| BDF, rtol 1e-06 | 120 | 1.2e-03 |
| BDF, rtol 1e-10 | 631 | 3.9e-07 |
| Radau, rtol 1e-06 | 100 | 7.9e-07 |
| Radau, rtol 1e-10 | 937 | 1.7e-11 |

#### Thermal runaway: u' = exp(u), u(0) = 0, up to t = 0.999 (the solution is -ln(1-t), infinite at t = 1)

| method | steps | error |
|---|---:|---:|
| GenSAS, errorTol 1e-06 | 20 | 7.0e-06 |
| GenSAS, errorTol 1e-10 | 26 | 1.8e-09 |
| RK45, rtol 1e-06 | 27 | 5.2e-05 |
| RK45, rtol 1e-10 | 132 | 3.4e-09 |
| DOP853, rtol 1e-06 | 18 | 1.6e-06 |
| DOP853, rtol 1e-10 | 42 | 2.2e-10 |
| LSODA, rtol 1e-06 | 67 | 6.6e-04 |
| LSODA, rtol 1e-10 | 182 | 1.6e-07 |
| BDF, rtol 1e-06 | 121 | 1.2e-03 |
| BDF, rtol 1e-10 | 614 | 2.8e-07 |
| Radau, rtol 1e-06 | 91 | 6.4e-07 |
| Radau, rtol 1e-10 | 845 | 2.0e-11 |

#### Kepler orbit with eccentricity 0.9, ten periods (the exact end state is the start state)

| method | steps | error |
|---|---:|---:|
| GenSAS, errorTol 1e-06 | 151 | 2.7e-05 |
| GenSAS, errorTol 1e-10 | 310 | 1.7e-08 |
| RK45, rtol 1e-06 | 607 | 6.7e-04 |
| RK45, rtol 1e-10 | 3732 | 5.9e-08 |
| DOP853, rtol 1e-06 | 290 | 1.1e-03 |
| DOP853, rtol 1e-10 | 806 | 6.0e-08 |
| LSODA, rtol 1e-06 | 2100 | 1.5e-03 |
| LSODA, rtol 1e-10 | 4565 | 8.9e-07 |
| BDF, rtol 1e-06 | 2758 | 9.4e-04 |
| BDF, rtol 1e-10 | 12298 | 1.2e-07 |
| Radau, rtol 1e-06 | 2179 | 1.5e-05 |
| Radau, rtol 1e-10 | 21505 | 1.4e-10 |

#### Lorenz system (chaotic), t = 0 to 10

| method | steps | error |
|---|---:|---:|
| GenSAS, errorTol 1e-06 | 134 | 6.2e-08 |
| GenSAS, errorTol 1e-10 | 272 | 9.3e-12 |
| RK45, rtol 1e-06 | 327 | 2.4e-04 |
| RK45, rtol 1e-10 | 2018 | 3.0e-08 |
| DOP853, rtol 1e-06 | 116 | 8.4e-05 |
| DOP853, rtol 1e-10 | 351 | 9.8e-10 |
| LSODA, rtol 1e-06 | 771 | 1.3e-04 |
| LSODA, rtol 1e-10 | 1601 | 2.9e-08 |
| BDF, rtol 1e-06 | 1014 | 1.2e-04 |
| BDF, rtol 1e-10 | 5018 | 1.9e-08 |
| Radau, rtol 1e-06 | 796 | 2.3e-06 |
| Radau, rtol 1e-10 | 7657 | 1.5e-11 |

#### Pendulum released at 3.0 rad (near the unstable position), t = 0 to 100

| method | steps | error |
|---|---:|---:|
| GenSAS, errorTol 1e-06 | 88 | 1.8e-03 |
| GenSAS, errorTol 1e-10 | 186 | 5.3e-08 |
| RK45, rtol 1e-06 | 422 | 3.0e-03 |
| RK45, rtol 1e-10 | 2682 | 4.4e-07 |
| DOP853, rtol 1e-06 | 136 | 8.9e-04 |
| DOP853, rtol 1e-10 | 386 | 3.3e-07 |
| LSODA, rtol 1e-06 | 1132 | 1.2e-02 |
| LSODA, rtol 1e-10 | 2382 | 3.4e-06 |
| BDF, rtol 1e-06 | 1222 | 8.9e-02 |
| BDF, rtol 1e-10 | 5429 | 1.8e-05 |
| Radau, rtol 1e-06 | 990 | 1.0e-04 |
| Radau, rtol 1e-10 | 9639 | 8.0e-10 |

#### Predator and prey (Lotka-Volterra), t = 0 to 100

| method | steps | error |
|---|---:|---:|
| GenSAS, errorTol 1e-06 | 91 | 8.1e-04 |
| GenSAS, errorTol 1e-10 | 192 | 1.0e-08 |
| RK45, rtol 1e-06 | 453 | 3.4e-04 |
| RK45, rtol 1e-10 | 2801 | 3.9e-09 |
| DOP853, rtol 1e-06 | 159 | 1.5e-06 |
| DOP853, rtol 1e-10 | 491 | 5.9e-09 |
| LSODA, rtol 1e-06 | 961 | 1.4e-03 |
| LSODA, rtol 1e-10 | 1916 | 9.2e-08 |
| BDF, rtol 1e-06 | 1626 | 3.3e-03 |
| BDF, rtol 1e-10 | 7598 | 1.4e-06 |
| Radau, rtol 1e-06 | 1416 | 4.7e-06 |
| Radau, rtol 1e-10 | 13857 | 9.6e-11 |

#### Van der Pol oscillator, mu = 100 (stiff between the jumps), t = 0 to 200

| method | steps | error |
|---|---:|---:|
| GenSAS, errorTol 1e-06 | 1165 | 8.2e-08 |
| GenSAS, errorTol 1e-10 | 1503 | 1.2e-11 |
| RK45, rtol 1e-06 | 11892 | 3.3e-07 |
| RK45, rtol 1e-10 | 25982 | 2.6e-11 |
| DOP853, rtol 1e-06 | 6087 | 4.3e-08 |
| DOP853, rtol 1e-10 | 9765 | 5.3e-12 |
| LSODA, rtol 1e-06 | 943 | 7.7e-06 |
| LSODA, rtol 1e-10 | 2684 | 9.2e-10 |
| BDF, rtol 1e-06 | 959 | 2.9e-06 |
| BDF, rtol 1e-10 | 4578 | 8.3e-10 |
| Radau, rtol 1e-06 | 688 | 3.2e-09 |
| Radau, rtol 1e-10 | 6808 | 3.3e-14 |

#### Robertson chemical kinetics (very stiff), written as a DAE with y3 = 1 - y1 - y2, t = 0 to 1000

| method | steps | error |
|---|---:|---:|
| GenSAS, errorTol 1e-06 | 80920 | 2.0e-06 |
| GenSAS, errorTol 1e-10 | 2 | failed: Error: the solver cannot advance beyond t = 0 within the error tolerance; the solution may be singular here |
| LSODA, rtol 1e-06 | 252 | 4.3e-09 |
| LSODA, rtol 1e-10 | 988 | 2.2e-10 |
| BDF, rtol 1e-06 | 192 | 8.2e-07 |
| BDF, rtol 1e-10 | 677 | 7.7e-10 |
| Radau, rtol 1e-06 | 84 | 1.1e-08 |
| Radau, rtol 1e-10 | 761 | 2.4e-12 |
