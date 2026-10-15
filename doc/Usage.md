# Usage

### PowerSAS
Current version only provides AC power flow. Assume the executable is `bazel-bin/app/app` (`bazel-bin/app/app.exe` on Windows, see [Installation](Installation.md)), then the power flow computation can be called as follows.

```bash
bazel-bin/app/app -p \
    -f/--file <input-file-name> \
    -o/--output <output-file-name> \
    [-l/--level <sas-order>] \
    [-s/--segment <segment-length>] \
    [-a/--alphatol <alpha-tolerance>] \
    [-d/--difftol <error-tolerance>] \
    [-r/--repeat <count>]
```

Explanations:
* The first argument `-p` (mandatory, and must be the first argument) means GenSAS runs under PowerSAS mode. For the rest of the arguments, the order does not matter.
* `-f/--file <input-file-name>` (mandatory) specifies the input data: a GenSAS case file (see [Case format](CaseFormat.md)), made from a PSAT .mat file with `python3 scripts/psat_to_case.py`. PSAT .mat files (for example `resources/psat_mat/d_014_syn_ind_zip_export.mat`) are not read directly.
* `-o/--output <output-file-name>` (mandatory) specifies the output curve data. Currently the output is a .mat file containing the power flow solution as a vector. The name is taken relative to the current directory, so do not pass an absolute path.
* `-l/--level <sas-order>` (optional) specifies the order of SAS. If not specified, the order of SAS is 15.
* `-s/--segment <segment-length>` (optional) specfies the length of a segment of SAS computation. If not specified, the the segment length is 1.0.
* `-a/--alphatol <alpha-tolerance>` (optional) specifies the tolerance of embedding variable in SAS. If not specified, the tolerance is set as 1e-4.
* `-d/--difftol <error-tolerance>` (optional) specifies the error tolerance of the equations. If not specified, the error tolerance is set as 1e-6.
* `-r/--repeat <count>` (optional) repeats the computation `<count>` times, e.g. for timing. If not specified, the computation runs once.

Example:
Try running power flow of the modified synthetic eastern-interconnection (EI) 70,000-bus system in the project root directory:
```bash
python3 scripts/psat_to_case.py resources/psat_mat/d_70k_070.mat d_70k_070.case
bazel-bin/app/app -p -f d_70k_070.case -s 0.5 -l 28 -d 1e-5 -o res.mat
```

### ModelicaSAS
ModelicaSAS simulates a single flat Modelica .mo model. Besides continuous `Real` variables, parameters, `der()` and algebraic equations, it supports a hybrid subset driven by time events: `Integer`/`Boolean`/`discrete` variables, `if` expressions and `if` equations whose conditions depend only on parameters and discrete variables, and `when time >= e` / `when sample(start, interval)` / `when initial()` (once, at the start time) equations with discrete assignments (`d = expr;`) and `reinit(x, expr);`. State events are supported too: `when` conditions and `if` conditions may depend on continuous variables or `time` (for example `when h < 0 then reinit(v, -e*pre(v)); end when;`). They are detected along the solution and located by bisection; `if` equations may change which equations are active (a variable can be a state in one mode and algebraic in another) as long as every mode involves the same set of variables; `abs`, `min`, `max` and `sign` may be used in continuous equations: they are turned into if-expressions and their switching points are state events, so the solution is exact (but a variable that hovers around a switching point, such as `sign(v)` friction at `v = 0`, makes the events pile up and stops with an error; `--piecewise smooth` replaces the functions by smooth approximations instead, see below); `edge(b)` (Boolean) and `change(d)` can be used in `when` conditions on a single discrete variable, and a simulation whose events accumulate in finite time (Zeno behavior) stops with an error. See `resources/mofile/test_hybrid_time_events.mo` and `resources/mofile/test_bouncing_ball.mo`. The elementary functions `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh`, `atan2`, `exp`, `log`, `log10` and `sqrt` can be used in equations, parameters and conditions, together with `^`; `atan2` in an equation is rewritten into `atan` forms with if-branches, so the switches at `x = 0` and at the branch cut (`y = 0`, `x < 0`) are state events; the solver does not check function domains, so a model that leaves one (`log` of a negative number, for example) gives wrong results or a failed solve. `noEvent(e)` in a continuous equation is accepted and treated like `e`. `terminate("message")` in a `when`-equation ends the run successfully at that event: the result stops there and records the message. `assert(condition, "message")` in the equation section is checked along the solution (a state event when the condition turns false) and in `when`-bodies at the event; a violation stops the run with an error naming the time, or with `AssertionLevel.warning` as a third argument it is reported once as a warning and the run goes on. `terminate` and `assert` are not supported inside `if`-equations. A binding `Real y = expr;` (also `Real y(start = 1) = expr;`) is the equation `y = expr`; it may refer to variables declared later. Bindings on `discrete`, `Integer` and `Boolean` variables are not supported (assign them in a `when`-equation). Several variables can be declared together (`Real x(start = 1), z;`), `input` (which needs a binding) and `output` prefixes are accepted, `//` and `/* */` comments are skipped, and `annotation(...)` clauses are ignored. Arrays are supported as sets of scalar variables: `Real x[3];`, `Real[2,3] A;`, `parameter Real k[3] = {1, 2, 3};`, `Real x[n](each start = 1)` or `start = {1, 2, 3}`. Subscripts are resolved while the model is read, so they must be constants: numbers, parameters declared earlier or loop variables (`x[i+1]`). Elements are named `x[2]` or `A[1,3]` in the results. `for i in 1:n loop ... end for;` (also `first:step:last` and several iterators `for i in 1:n, j in 1:m loop`) is unrolled in `equation` and `initial equation` sections and inside `if`/`when` branches, and `size(x, k)`, `sum(x)` and `prod(x)` are available; whole-array expressions (`x = y`, `x + y`) are not, write a loop instead. See `resources/mofile/test_array_cascade.mo`. `extends`, `algorithm` sections, functions and nested classes are not supported and give an error. Algebraic equations may be implicit (`y*y = x`); they are solved with Newton's method from the start values, so give the unknowns start values near the solution (a singular Jacobian is reported as an error). The simulation can be called as follows:

```bash
bazel-bin/app/app -g \
    -m/--mode file/string \
    -i/--input <input> \
    [-o/--output <mat-file-name>] \
    [-j/--json <json-output-file-name>] \
    [-c/--csv <csv-output-file-name>] \
    [--piecewise events|smooth] [--smooth-eps <eps>] \
    [-l/--level <sas-order>] \
    [-s/--segment <segment-length>] \
    [-a/--atol <alpha-tolerance>] \
    [-e/--etol <error-tolerance>] \
    [-t/--time <max-time>] \
    [-k/--step <output-step>] \
    [-v/--verbose]
```

Explanations:
* The first argument `-g` (mandatory, and must be the first argument) means GenSAS runs under ModelicaSAS mode. For the rest of the arguments, the order does not matter.
* `-m/--mode file/string` (mandatory) specifies the format of input. If `-m file`, then expect the input file name after `-i`; and if `-m string`, then expect the Modelica model content after `-i`.
* `-i/--input <input>` (mandatory) specifies the input data. The `<input>` depends on the mode specified after `-m`.
* `-o/--output <mat-file-name>` (optional) writes the result as a MATLAB level 5 .mat file (no HDF5 needed to read it). It holds the sample times `t`, the continuous variables `s` (one row per variable, in the order of the char matrix `names`) and, for hybrid models, the discrete variables `d` (names in `dnames`) and the event times `event_t`.
* `-j/--json <json-output-file-name>` (optional) specifies the name of the .json file containing output curve. For hybrid models the file also has a `discrete` object (the discrete variables, sampled like the curves; the two samples at an event time hold the values before and after the event) and an `events` array (`time`, `assignments` and `reinits` of each event).
* `-c/--csv <csv-output-file-name>` (optional) writes a CSV file with one column per variable (`time`, continuous variables, then discrete variables). For hybrid models with events, the events go to a second file named like the CSV file with `.events.csv` instead of `.csv`.

* `--piecewise events|smooth` (optional, default `events`) selects how `abs`, `min`, `max` and `sign` in continuous equations are treated. `events` is exact: the switching points are located as state events. `smooth` uses `abs(e) = sqrt(e^2+eps^2)`, `max(a,b) = (a+b+sqrt((a-b)^2+eps^2))/2` (`min` with a minus sign) and `sign(e) = e/sqrt(e^2+eps^2)`: there are no events, a model that would chatter at a switching point settles near it, and the result is off by an amount of the order of `eps` near the switching points (a very small `eps` makes the solver take small steps there).
* `--smooth-eps <eps>` (optional, default 1e-6) is the smoothing parameter.

All result files share one schema, version 1: see [Result format](ResultFormat.md).
* `-l/--level <sas-order>` (optional) specifies the order of SAS. If not specified, the order of SAS is 15.
* `-t/--time <max-time>` (optional) specifies the maximum time of simulated process. If not specified, the time is set as 10.0.
* `-s/--segment <segment-length>` (optional) specfies the length of a segment of SAS computation. If not specified, the the segment length is 1.0.
* `-a/--atol <alpha-tolerance>` (optional) specifies the tolerance of embedding variable in SAS. If not specified, the tolerance is set as 1e-3.
* `-e/--etol <error-tolerance>` (optional) specifies the error tolerance of the equations. If not specified, the error tolerance is set as 1e-5.
* `-k/--step <output-step>` (optional) specifies the time step of the output curves. If not specified, the time step is set as 0.01.
* `-v/--verbose` (optional) if used, will print intermediate result in SAS computation.

Example:
Try running simulation of the model in `resources/mofile/test_solve_ode.mo` in the project root directory.

```bash
 bazel-bin/app/app -g \
    -m file -i resources/mofile/test_solve_ode.mo \
    -o resources/mofile/test_solve_ode.mat \
    -t 15
```

#### Job files
The same simulation can be described in a JSON job file and started with `bazel-bin/app/app run job.json [--set key=value] [--validate]`; paths in the file are relative to the job file. See [Job format](JobFormat.md). The flags above remain available. A missing or non-numeric flag value is an error (exit code 1), and an unknown flag is rejected. The last sample of every result is at the end time.

### Paths on different operating systems
Run the commands from the repository root, using relative paths; they behave identically on Linux, macOS and Windows (MSYS2 shell, with `app.exe`). The native Windows executable cannot open MSYS-style absolute paths like `/d/work/...`; convert them with `cygpath -m <path>`.
