# Roadmap and open topics

This page collects what is not built yet and the thinking about it, so that the decisions are visible and can be changed. Nothing here is promised. What exists is described in the other pages: the [Modelica subset](Modelica.md), the [server](Server.md) and the [web client](Web.md).

## Large systems

The target of GenSAS is large sparse systems, power grids first of all. The engine already solves sparse DAEs, and the language has arrays, `for` equations, components and connectors. What is missing is mostly around the solver: how a model of thousands of components is written, simulated and looked at.

### Where it breaks first

* **Writing the model.** A grid of a thousand buses cannot be typed. The model has to come from data: array and matrix parameters filled from a table or a data file, or a model generated from a [case file](CaseFormat.md) (which already holds the power-flow data). Text of megabytes also raises the limits on the size of a request (`maxSourceBytes`) and makes a text editor slow.
* **Result size.** Ten thousand variables sampled a thousand times are 80 MB. The JSON result that the client plots does not scale to that. The result request should name the variables it wants (or a pattern such as `gen[*].omega`), the server should thin the samples (keeping the minimum and maximum of each bucket, so that spikes survive), and the transfer should be a binary column format that is loaded per variable.
* **The solver.** Nobody knows yet where the time goes for 10 000 equations: assembling the Jacobian, the sparse factorization, the series arithmetic, or the number of steps. This has to be measured before it is designed around (see below).
* **Initial conditions.** A grid study starts from a power flow. The power-flow code in `pf/` could set the start values of the dynamic model; today the start values are typed in.
* **Hosting.** A public server needs a limit on the time and memory of a run, and a way to say "this model is too big" before the run starts (a model summary with the counts of states and equations in the check result).

### What the client needs

* A **tree of the variables** by component name (`a.b.x` below `a.b` below `a`), with search and filters, and lists that draw only the visible rows. A flat list of checkboxes is fine for ten variables and useless for ten thousand.
* A **parameter panel grouped by component**, with search and a way to change many values at once (`gen[*].H = 5`).
* **Plots of the selection only**, never of everything; the choice of the variables is part of the result request.
* **Progress and partial results.** Time to go, and curves that grow while the run goes on. The engine can hand over each solution segment as it is finished; the server and the client do not use that yet.

### How to proceed

1. **Measure.** Write a scalable test model (N machines of the classical swing model, coupled through a network with algebraic equations, N from 10 to 1000 and more). Time it, find the largest case that is comfortable, and profile the steps. [How the SAS solver compares with conventional integrators](SasComparison.md) says why such a measurement is needed against a conventional DAE solver: the claim of the method is about this kind of model.
2. **Results.** Variable selection, decimation and a binary format in the result API; selective loading in the client.
3. **Data in models.** Array parameters from data, and a way to build a model from a case file.
4. **Client.** The variable tree, the grouped parameters, the model summary.
5. **Power-flow start values.**

## Run time

[The comparison with conventional integrators](SasComparison.md) shows that GenSAS needs fewer steps for a given error, but that a step is expensive: 20 to 100 microseconds for a few states, against well below one microsecond for a compiled Runge-Kutta step. Against the BDF solver of Modelica tools the run time is comparable, against explicit methods on non-stiff models it is 50 to 200 times worse. Making the step cheaper is the main way to turn the advantage in steps into an advantage in time, and it needs no change of the method. In the order of the profile:

* **The trial evaluations of the step size search.** The solution is evaluated at many trial times; each evaluation of the Pade approximant allocates temporaries. (The derivative of the solution, needed for the residual, is no longer rebuilt at every trial; that made the runs 2 to 2.5 times faster.) Evaluating on plain arrays, and caching what does not depend on the trial time, should remove most of the 60% of the time that go there.
* **One factorization per step for larger systems (done).** Algebraic systems of up to 64 unknowns are factored densely once per step; larger ones are factored once per step with SuperLU (`sas/SparseLU`) and each order is a triangular solve. The rings and networks of machines with 120 to 840 algebraic unknowns became 2.1 to 2.5 times faster. For the networks of the next section the remaining costs are the series arithmetic and the Pade approximants of all rows.
* **The interpretation of the model.** The series arithmetic walks the expression tree of the equations for every coefficient. Compiling the equations into a flat list of operations, or generating code, would cut the cost of the nonlinear terms.
* **Larger models: measure the scaling.** The ring of machines in the comparison grows like the size to the power 1.2 in time; the number of steps does not depend on the size. A network with algebraic equations, compared with an implicit DAE solver with a sparse linear solver, is the measurement that is still missing.
* **Accuracy at the top.** `errorTol` stops at 1e-10 and long runs end with errors of 1e-8 (the Kepler orbit), where `DOP853` reaches 1e-10. Whether a smaller tolerance and a higher order are enough, or the accumulation of the error needs attention, has not been looked at.

## Several runs and parameter sweeps

Typical uses: how does the result depend on one or two parameters; scenarios (the clearing time of a fault); random samples with a fixed seed (Monte Carlo); sensitivity; fitting a parameter to data; the pictures that need many runs, such as a bifurcation diagram.

Each run is independent, so the problem is bookkeeping and presentation, not the solver. The engine is re-entrant, and the server already runs jobs in separate worker processes with limits.

### Server

A **batch** resource next to the runs:

* A request holds the base run and the sweep: for each parameter a list of values or a range (linear or logarithmic, with a count), combined as a product (a grid) or element by element, or random samples from a distribution with a seed.
* The runs go through the existing pool of workers with the same limits. A batch has a progress of n of N runs, can be cancelled as a whole, and a failed run does not stop the others; its state and message are kept.
* Keeping all curves of 500 runs is not affordable, so a batch keeps what was asked for: **summary values** per run (final value, minimum, maximum, time of the maximum, time when a threshold is crossed), and the full curves only for a chosen few or on request. This needs the selective results of the section above.
* A public server needs an upper limit on the size of a batch and on the share of the workers one batch may take.

### Client

* A parameter can be marked "vary": a range with a count and a scale, or a list. A second parameter makes a grid.
* The result of a batch is shown as **curves of all runs on one plot**, colored by the value of the parameter, with a color bar; as a **summary value against the parameter**; and for two parameters as a **heat map** of a summary value. Runs that failed, or ended in a Zeno point, are marked in the picture.
* The first step, which is also useful alone, is to **compare runs**: keep the result of a finished run and draw a new one over it.

### Scripts

`gensas.sweep(job, {"k": [...]}, processes=n)` in the [Python interface](Python.md) for work that is better done in a notebook. Sweeps are the same list of jobs with different `model.parameters`, so this is a thin layer.

### Order

Compare runs, then selective results (shared with large systems), then one-parameter sweeps with overlay and summary plots, then grids and heat maps, then random samples and the Python function.

## Known problems and small fixes

* **Smoothed `abs` can leave the domain of `sqrt`.** With `--piecewise smooth` the Chua circuit stops at t = 0.118 with a message about `sqrt`: the series of `sqrt(e*e + eps*eps)` is evaluated outside its domain. The exact mode (events) works.
* **Very stiff problems are slow.** The Robertson problem needs about 100 000 steps (see [the comparison](SasComparison.md)). The cause is [the radius of convergence of the series](SasComparison.md#why-stiff-problems-cost-many-steps): the step is about 40 over the fast eigenvalue whenever the fast dynamics are nonlinear. Ideas, none tried: a higher default order for stiff models, an exponential integrating factor for the fast linear part, and measuring whether the fast modes of a power system are of this kind.
* **A Zeno point ends the run with an error.** The message explains it, but the run could keep its results up to the accumulation point, or treat the accumulation point as the end of a run that ended normally. Not decided.
* **A model with an unbound `input` cannot be run on its own.** The message is correct, but the web client could mark such models as components in the model selector, or point at the declaration.
* **Algebraic variables need a start value when the equations cannot be evaluated at 0.** A division by an algebraic variable (`r` in a Kepler model) fails at t = 0 with the message that the equations cannot be evaluated; the fix is `Real r(start = 1)`. The message could say that.
* **A state without a start value is reported in plain models only.** Hybrid models do not report it.
* **The solver has trouble with terms that cancel by a factor of zero** (`der(v) = -f(v)*0 - 1`).
* **Small gaps in the language:** `while` loops with constant bounds, poles of even order (`1/u^2`) are not detected, the time of a domain error inside a hybrid run is only known to the step, and a failed run may not keep the results of its steps.
* **Several files in one project** are not supported by the client (a single text only).
