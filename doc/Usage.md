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
* `-f/--file <input-file-name>` (mandatory) specifies the input data. Currently the input file needs to be .mat file containing PSAT data structure. See `resources/psat_mat/d_014_syn_ind_zip_export.mat` for example.
* `-o/--output <output-file-name>` (mandatory) specifies the output curve data. Currently the output is a .mat file containing the power flow solution as a vector. The name is taken relative to the current directory, so do not pass an absolute path.
* `-l/--level <sas-order>` (optional) specifies the order of SAS. If not specified, the order of SAS is 15.
* `-s/--segment <segment-length>` (optional) specfies the length of a segment of SAS computation. If not specified, the the segment length is 1.0.
* `-a/--alphatol <alpha-tolerance>` (optional) specifies the tolerance of embedding variable in SAS. If not specified, the tolerance is set as 1e-4.
* `-d/--difftol <error-tolerance>` (optional) specifies the error tolerance of the equations. If not specified, the error tolerance is set as 1e-6.
* `-r/--repeat <count>` (optional) repeats the computation `<count>` times, e.g. for timing. If not specified, the computation runs once.

Example:
Try running power flow of the modified synthetic eastern-interconnection (EI) 70,000-bus system in the project root directory:
```bash
bazel-bin/app/app -p -f resources/psat_mat/d_70k_070.mat -s 0.5 -l 28 -d 1e-5 -o res.mat
```

### ModelicaSAS
ModelicaSAS simulates a single flat Modelica .mo model. Besides continuous `Real` variables, parameters, `der()` and algebraic equations, it supports a hybrid subset driven by time events: `Integer`/`Boolean`/`discrete` variables, `if` expressions and `if` equations whose conditions depend only on parameters and discrete variables, and `when time >= e` / `when sample(start, interval)` equations with discrete assignments (`d = expr;`) and `reinit(x, expr);`. State events are supported too: `when` conditions and `if` conditions may depend on continuous variables or `time` (for example `when h < 0 then reinit(v, -e*pre(v)); end when;`). They are detected along the solution and located by bisection; `if` equations may change which equations are active (a variable can be a state in one mode and algebraic in another) as long as every mode involves the same set of variables; `edge()`/`change()` are not supported yet, and a simulation whose events accumulate in finite time (Zeno behavior) stops with an error. See `resources/mofile/test_hybrid_time_events.mo` and `resources/mofile/test_bouncing_ball.mo`. The simulation can be called as follows:

```bash
bazel-bin/app/app -g \
    -m/--mode file/string \
    -i/--input <input> \
    -o/--output <output-file-name> \
    [-j/--json <json-output-file-name>] \
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
* `-o/--output <output-file-name>` (mandatory) specifies the output curve data file name. Currently the output is a .mat file containing the output as a matrix.
* `-j/--json <json-output-file-name>` (optional) specifies the name of the .json file containing output curve. For hybrid models the file also has a `discrete` object (the discrete variables, sampled like the curves; the two samples at an event time hold the values before and after the event) and an `events` array (`time`, `assignments` and `reinits` of each event). 
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

### Paths on different operating systems
Run the commands from the repository root, using relative paths; they behave identically on Linux, macOS and Windows (MSYS2 shell, with `app.exe`). The native Windows executable cannot open MSYS-style absolute paths like `/d/work/...`; convert them with `cygpath -m <path>`.
