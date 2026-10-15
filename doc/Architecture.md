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

## Testing

* Unit tests live in `sas/test` and run in CI with `bazel test //...`.
* To iterate without Bazel, use `scripts/local_solver_test.sh`. It needs a C++ compiler, Armadillo,
  SuperLU, jsoncpp and BLAS/LAPACK (on Debian/Ubuntu:
  `sudo apt install g++ libarmadillo-dev libsuperlu-dev libjsoncpp-dev libblas-dev liblapack-dev`)
  A full run takes seconds after the first compile.
* `scripts/run_integration_test.sh` runs the command line on a power-flow case and on the bouncing
  ball example (used by the Integration workflow on pushes to `main`).
* Formatting: `bash scripts/format.sh` (clang-format 22.1.2, as in the lint check).
