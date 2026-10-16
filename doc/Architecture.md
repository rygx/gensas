# Architecture of the Modelica path (ModelicaSAS)

This page describes how a `.mo` model becomes a simulation result, which parts are hybrid-aware, and
how to test the solver quickly. Power flow (`-p`) is a separate path (`pf/`, `io/`) and is not covered.

## Pipeline

```
.mo text
  -> SasModelParser (sas/ModelParser.cpp, ExprParser.cpp, TokenScanner.cpp)   parse into SasModel / AST
  -> checkModel     (sas/ModelChecker.cpp)                                    what can the solver run?
  -> plain model:   SasComputationModel::solve                                one DAE, Pade segments
     hybrid model:  HybridSimulator::run (sas/HybridSimulator.cpp)            events between segments
  -> SasSolutionSet (sas/SasComputation.cpp)                                  buildResult -> JSON / CSV / MAT writers (sas/SasResult*.cpp)
```

`app/main.cpp` (`-g` mode) chooses the path: a model with discrete variables, if/when statements or
`time` (`isHybridModel`) goes through the simulator, everything else straight to the solver.

### Parsing

* The model parser accepts a flat hybrid subset: `Real`/`Integer`/`Boolean`/`discrete` declarations
  with `start`/`unit`/`fixed` modifications, parameters with constant expressions, `der()` and algebraic
  equations, if-expressions, if/when equations (with `elseif`/`elsewhen`), `reinit`, `pre`, `sample`,
  and initial equations of the form `x = const` / `der(x) = const`.
* `if` and `when` statements are AST nodes with `nodeKind` `NK_IfEquation` / `NK_WhenEquation`; their
  layout is documented next to the enum in `sas/SasModel.h`. Statement lists are linked through
  `AstNode::next`.
* Identifiers are injected into the model's id table by the scanner; attribute names inside a
  modification are scanned with `noInject`.

### Checking

`checkModel` takes a `SolverFeatures` set. With every flag false it describes the plain solver
(continuous `Real`, parameters, `der`, algebraic equations). The application enables the hybrid
flags and `timeEventsOnly` + `stateEvents`, which additionally forbid constructs the simulator
cannot handle: relational operators in continuous equations, `edge`/`change`
outside when-conditions (and of anything but one discrete variable), `when` inside `if`, assignments to discrete variables outside `when`.

### Piecewise functions

`rewritePiecewiseFunctions` (run by the application before `isHybridModel`/`checkModel`) rewrites
`abs`, `min`, `max` and `sign` in continuous equations into if-expressions
(`abs(e)` becomes `if e >= 0 then e else -e`). The conditions are ordinary dynamic guards, so the
switching points are located like any other state event and every mode is smooth for the Pade solver.
Conditions and when-bodies are not rewritten; they are evaluated directly.
With `PW_SMOOTH` (`--piecewise smooth`) the calls become smooth `sqrt(e*e + eps*eps)` forms instead, with
no events.

### Elementary functions

`sin`, `cos`, `exp`, `sqrt`, `log`, `atan`, `asin` and `acos` are solved natively: each is a temporary variable
with an algebraic equation `f(x) = y` and its own series recurrence (`log`, `atan`, `asin` and `acos` use
`y' = x' g(x)^alpha`). `tan`, `sinh`, `cosh`, `tanh` and `log10` are rewritten before solving, together with the
piecewise functions, into `sin/cos`, `exp` and `log` expressions. `noEvent(e)` in a continuous equation is
replaced by `e`: events for the relations inside are still generated, which the Modelica specification allows.
`evalMathFunction` (`MathFunctions.h`) evaluates all of them numerically for parameters, start values,
conditions and event expressions.

### Arrays and for-equations

Arrays never reach the solver. `parseDeclarator` declares one scalar variable per element, named `x[1]` or
`A[2,3]` (the base name is recorded in `SasModelParser::arrayDims`). Each equation section statement is then
*elaborated* (`elabChain` in `ModelParser.cpp`): loop variables become integer constants, `x[i+1]` is
evaluated with `evalConst` and replaced by the identifier of that element, `size`/`sum`/`prod` of an array are
expanded, and an `NK_ForEquation` node is unrolled into its statements (one `AstTree` each). The same
elaboration resolves subscripts in parameter values, start values and bindings. The iterator and element
counts are limited (1e5 elements, 1e6 unrolled iterations).

### assert and terminate

A top-level `assert(cond, msg [, level])` is a guard of kind `G_ASSERT` whose truth value is the violation
(`not cond`); the simulator treats it like a when-condition, so crossings are found by the same sampling and
bisection, and a violation either fails the run (`error`) or is added to `HybridSimulator::warnings()`. The
call is dropped by `specializeModel`. `assert` and `terminate` inside a `when` body run in `fireEvents`;
`terminate` sets `HybridEvent::terminate`, the simulator stops after that event, and the solution set (and the
result) records the time and message. Both are rejected inside if-equations by the checker.

### Hybrid simulation

`HybridSimulator::run` repeats:

1. Find the next **time event** (`nextEventTime`: `time >= e`, `sample(start, interval)`).
2. Build a fresh `SasComputationModel` from the parsed model and **specialize** it
   (`specializeModel`): `when` statements are dropped, if-equations and if-expressions are resolved
   for the current discrete state, discrete variables are replaced by their values. Conditions that
   depend on continuous variables or `time` ("dynamic guards") take their truth value from the
   simulator's tracked state.
3. Solve until the event time. The `onSolution` callback scans each new solution segment: it samples
   the dynamic guards (`guardSamples`, finer right after an event) and locates a flip by bisection on
   the Pade solution; the solve stops there. Between two samples it also follows the margin (left minus
   right side) of every relation: a parabola through the interval tells whether the margin may cross
   zero twice, and the extremum is then searched and used as an extra sample, so pulses shorter than
   the sample spacing are found.
   `edge(b)`/`change(d)` are true only while the event instant is handled (`HybridEnv::atEvent`),
   compared with the discrete values at the start of the instant (`eventPre`).
4. At the event, `fireEvents` runs the `when` bodies with the values before the event (`pre`),
   applies discrete assignments and returns `reinit`s; conditions made true by the event fire in
   further rounds (up to 50). Events that pile up in finite time stop the run with a Zeno error.
5. Continue from the new state. Continuous values are carried over by name, so a mode may reorder
   the equations or turn a state into an algebraic variable; the set of variables must stay the same.

Guards are numbered by `collectGuards`; the id is stored in `AstNode::index` of the if/when node so it
survives copying the AST into the computation model.

### Solver notes (`sas/SasComputation.cpp`)

* `solve()` computes Pade segments of length up to `segment`; `SasSolutionSet` keeps all of them with
  their absolute `[absStart, absEnd]`. `startTime` and the `onSolution` callback exist for the hybrid
  simulator.
* A `CheSolutionPade` caches its coefficients on the first evaluation. Anything that changes the series
  afterwards must reset `ready` (the AE calibration in `solveSegment` does).
* Nonlinear operations are decomposed into temporary variables (`modelFlatten`), so `y = k*x` stays
  linear but `x*y` or `sin(x)` add algebraic equations.
* `copyAstTree` copies `next`-linked nodes with their subtrees; call arguments and statement lists
  rely on that.
* Segments shorter than `alphaTol` are not solved.

### Library entry point

`sas/Simulation.h` wraps the whole pipeline for one Modelica run: `ModelicaJob` (model file or text, solver options, piecewise mode, output step and files; filled by `parseModelicaJob` from a job file, see [Job format](JobFormat.md)), `Simulation::run(progress)` (parse, check, solve, sample; the optional progress callback receives the time reached and the end time after each solution segment and returns true to cancel) and `Simulation::result()`. The command line is a thin layer over it: `app run` (`runJobFile` in `app/main.cpp`) reads a job file, applies `--set` overrides and calls `parseModelicaJob`; `app -g` builds the same `ModelicaJob` from flags; both then go through `runModelica`, and the C API (`capi/gensas.h`, see [C API](CApi.md)) wraps the same class.

### Re-entrancy

There is no process-wide mutable state in the Modelica path: each `SasModelParser` owns a `ModelPool`
(the top-level models it creates), string literals live in the model that owns the AST pointing at
them, and every simulation works on its own `SasComputationModel`. Independent runs (parse, check,
solve) may therefore execute on different threads at the same time; `concurrent_run_test` does this.
A single model, parser or simulator is not thread safe and must stay on one thread.

## Testing

* Unit tests live in `sas/test` and run in CI with `bazel test //...`.
* To iterate without Bazel, use `scripts/local_solver_test.sh`. It needs a C++ compiler, Armadillo,
  SuperLU, jsoncpp and BLAS/LAPACK (on Debian/Ubuntu:
  `sudo apt install g++ libarmadillo-dev libsuperlu-dev libjsoncpp-dev libblas-dev liblapack-dev`)
  A full run takes seconds after the first compile.
* Coverage: `bash scripts/coverage.sh` (needs `pip install gcovr coverage`) runs the same tests and the
  Python tests with `--coverage` and writes line and branch coverage of `sas/`, `io/`, `util/` and
  `capi/` plus line coverage of `python/gensas` to `/tmp/gensas-coverage` (`summary.md`, HTML reports).
  The Coverage workflow does this on every pull request, shows the numbers in the job summary, comments
  them on the pull request and keeps the HTML reports as an artifact. Power flow and the PSAT reader are
  built and tested through Bazel only, so the script does not see them.
* `scripts/run_integration_test.sh` runs the command line on a power-flow case and on the bouncing
  ball example (used by the Integration workflow on pushes to `main`).
* Formatting: `bash scripts/format.sh` (clang-format 22.1.2, as in the lint check).
