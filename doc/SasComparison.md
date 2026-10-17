# How the SAS solver compares with conventional integrators

This page says where the semi-analytical solver (SAS) of GenSAS is better than the integrators in common use, where it is not, and how to check it yourself. The numbers come from [scripts/compare_with_scipy.py](../scripts/compare_with_scipy.py), which runs GenSAS and the methods of SciPy on the same small problems.

## What is different

A conventional integrator takes a step by combining a few evaluations of the right-hand side (Runge-Kutta), or by solving a nonlinear equation with Newton's method for the new state (BDF, Radau, the trapezoidal rule). SAS builds the solution of the step as a power series in time, one coefficient after another, each from a linear system with the same sparse matrix at every order, and turns the series into a Pade approximant. The result is:

* **A solution that is a function, not a list of points.** The approximant can be evaluated at any time in the step. Output samples, the check of equation residuals, and the search for event times all use this function; none of them needs extra steps or an interpolation of lower order than the solver.
* **No Newton iteration in time stepping.** Algebraic equations of a DAE are part of the series recurrence, so a step solves a sequence of linear systems with the same matrix, not a nonlinear iteration that may fail to converge.
* **A step size that follows the analytic structure of the solution.** The Pade approximant is useful beyond the radius of convergence of the series, and the step is chosen from the residual of the equations, not from a stability limit of the method.
* **A high and adjustable order** (default 15, up to 1000) instead of the fixed order of a Runge-Kutta pair.

The price is a more expensive step: one linear solve per order plus the series arithmetic of the nonlinear terms, and a model with many nonlinear functions has more series to carry.

## What was measured

Every problem is solved with GenSAS at two error tolerances and with SciPy's `RK45`, `DOP853` (explicit), `LSODA`, `BDF` and `Radau` (implicit, the family that Modelica tools use through DASSL and CVODE) at two relative tolerances. The reference is the exact solution where there is one, and a solution with a tolerance of 1e-13 otherwise. The error is the largest error of the final state, scaled by `max(|reference|, 1)`.

Read the tables with these points in mind:

* **Steps, not seconds.** One step of GenSAS does more work than one step of a Runge-Kutta method, and a comparison of C++ and Python timings says nothing about the methods. Steps show how much the method needs to follow the solution; the cost of a step depends on the model. For the models here a GenSAS run takes a few milliseconds to a few tenths of a second.
* **The tolerances are not the same quantity.** `errorTol` of GenSAS bounds the residual of the equations in a step; `rtol` of SciPy bounds an estimate of the local error. Compare the errors that were reached, not the tolerances. In both cases the error at the end of a long run is larger than the tolerance.
* **Numbers are from one machine and one version** (GenSAS from the main branch of October 2026, SciPy 1.15). The pattern matters, not the third digit.

## Results in short

Where the solution is smooth and the demand on accuracy is high, GenSAS needs far fewer steps than the methods that are usually applied to such models:

* Lorenz system with an error tolerance of 1e-10: 272 steps and an error of 9e-12. `RK45` needs 2018 steps for 3e-8, `LSODA` 1601 for 3e-8, `BDF` 5018 for 2e-8, `Radau` 7657 for 2e-11, and the eighth order `DOP853` 351 for 1e-9.
* Kepler orbit with eccentricity 0.9 over ten periods, 1e-10: 310 steps, error 2e-8; `DOP853` 806 steps for 6e-8, `RK45` 3732 for 6e-8, `BDF` 12298 for 1e-7.
* The pendulum and the predator-prey model show the same: at a tolerance of 1e-10, GenSAS takes about 190 steps, `DOP853` 390 to 490, `RK45` 2700 to 2800, the implicit methods 1900 to 13900.
* The cost grows slowly with the demanded accuracy. Lowering the tolerance from 1e-6 to 1e-10 doubles the steps of GenSAS on the Lorenz problem; it multiplies the steps of `RK45` by six and those of `Radau` by ten.

Where the solution approaches a singularity, the series method adapts to it:

* `x' = x^2` up to t = 0.99, where the solution is 100 and infinite at t = 1: 2 to 3 steps. `DOP853` needs 13 to 40, `RK45` 32 to 197, `BDF` 120 to 631.
* `u' = exp(u)` is the same kind of problem, and here the advantage is mostly gone: 20 steps for GenSAS and 18 for `DOP853` at 1e-6, 26 and 42 at 1e-10.

Event times are located on the continuous solution. For the pendulum above, which passes through zero six times in 100 s, GenSAS finds the crossings with an error of 8e-8 s in 193 steps; `DOP853` needs 386 steps for 3e-7 s, `RK45` 2682 steps for 4e-7 s, `LSODA` 2382 steps for 3e-6 s.

The output grid is free. The Lorenz run above takes the same 190 steps whether the result has 10, 1000 or 10000 samples; the time goes from 19 ms to 27 ms (measured at an error tolerance of 1e-8).

## Where SAS is not the better choice

* **Very stiff problems with a fast transient.** The Robertson problem (a stiff chemical reaction written as a DAE) needs 80920 steps for GenSAS at 1e-6, against 84 to 252 for the implicit methods, and GenSAS cannot start it at 1e-10. The cause has not been analyzed; the steps are short from the first one. Explicit methods are much worse (`RK45`, which is not in the table, needs about 1.7 million steps), so GenSAS is in between, but not a good choice for this kind of model.
* **Stiffness between jumps is fine.** The Van der Pol oscillator with mu = 100 needs 1165 steps at 1e-6, close to `LSODA` and `BDF` (943 and 959) and a tenth of `RK45` (11892).
* **Models whose solution is not smooth.** At each switching point the solver stops and restarts, and a model that switches thousands of times costs thousands of restarts. A Zeno model is reported as an error.
* **Cost per step.** If the stable step of a simple explicit method is as long as the series step, the explicit method wins because its step is cheap.
* **`piecewise = smooth` and `abs`.** The smoothed form of `abs` can leave the domain of `sqrt` while the series are evaluated (the Chua circuit stops at t = 0.12 in this mode). The exact mode, which uses events, works.

## Examples that show the strengths

These problems are in the script, and the first two are good candidates for examples in the web client:

1. **High accuracy on a smooth nonlinear model:** the Lorenz system or the Kepler orbit at a tolerance of 1e-10. This is where the high order pays off.
2. **A solution that runs into a singularity:** `x' = x^2` or the thermal runaway `u' = exp(u)` ended just before the blow-up.
3. **Accurate event times of a nonlinear model:** the pendulum crossings above.
4. **Dense output:** any run where many samples are wanted.

What these problems do not show is the case that GenSAS was built for: **large sparse DAEs, in particular power systems**. There the absence of a Newton iteration in the time stepping matters most, and the publications listed in the [README](../README.md) report the results. A comparable measurement in this repository needs a scalable test model (a network of N machines with algebraic network equations) and a conventional DAE solver to compare with; see [Roadmap](Roadmap.md#large-systems).

## Reproducing the numbers

Build the application (see [Installation](Installation.md)), install `numpy` and `scipy`, and run from the repository root:

```bash
python3 scripts/compare_with_scipy.py --app bazel-bin/app/app
python3 scripts/compare_with_scipy.py --app bazel-bin/app/app --problems lorenz kepler events
python3 scripts/compare_with_scipy.py --app bazel-bin/app/app --tolerances 1e-4 1e-8
```

The script runs GenSAS through `app -g` with `-v`, counts the lines `Stage start` for the steps, and reads the final state from the CSV file.

## Full results

### Finite-time blow-up: x' = x^2, x(0) = 1, up to t = 0.99 (the solution is 1/(1-t), infinite at t = 1)

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

### Thermal runaway: u' = exp(u), u(0) = 0, up to t = 0.999 (the solution is -ln(1-t), infinite at t = 1)

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

### Kepler orbit with eccentricity 0.9, ten periods (the exact end state is the start state)

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

### Lorenz system (chaotic), t = 0 to 10

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

### Pendulum released at 3.0 rad (near the unstable position), t = 0 to 100

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

### Predator and prey (Lotka-Volterra), t = 0 to 100

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

### Van der Pol oscillator, mu = 100 (stiff between the jumps), t = 0 to 200

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

### Robertson chemical kinetics (very stiff), written as a DAE with y3 = 1 - y1 - y2, t = 0 to 1000

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

### Event times: the pendulum above passes th = 0 (downwards) 6 times in 100 s

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
