# GenSAS

**GenSAS** is a C++ based generalized simulation tool based on semi-analytical solutions (SAS) technology. It provides two modules:

1. **ModelicaSAS** is a SAS-based simulation tool that supports the Modelica modeling language. It simulates a single flat Modelica model: ODEs and implicit algebraic equations, `algorithm` sections, user-defined functions, arrays (also as whole-array equations) and `for` loops, the usual elementary functions, and hybrid models with discrete variables, `if` and `when` equations, time and state events, `reinit`, `terminate` and `assert` (see [Supported Modelica](doc/Modelica.md)).

2. **PowerSAS** is a robust, efficient and scalable power grid analysis framework based on semi-analytical solutions (SAS) technology. This is the version for C++ users. It currently provides AC power flow functionality, reading PSAT-style cases (converted once to the compact [GenSAS case format](doc/CaseFormat.md)) and writing JSON, CSV or MAT results.

Backed by the SAS approach, the GenSAS tool provides much better convergence than the tools using traditional Newton-type equation solvers. Moreover, due to the analytical nature, GenSAS provides model-adaptive high-accuracy approximation, which brings significantly extended effective range for power grid analysis. PowerSAS has been used to solve large-scale system cases with 200,000+ buses.

### History and status
GenSAS was created at Argonne National Laboratory, where Rui Yao was the principal investigator of the project before 2022. Since then Rui Yao has continued the development independently in this repository, which started as a fork of the original [ANL-CEEESA/gensas](https://github.com/ANL-CEEESA/gensas) repository and keeps its history. The original work, its funding and its publications are acknowledged below.

Recent work in this repository includes Linux, macOS and Windows builds checked in CI, hybrid Modelica simulation with event handling, a self-contained input and output stack (no MATLAB, HDF5 or zlib libraries needed), a solver that can run several independent simulations in one process, JSON job files (`app run job.json`), a C library interface with a Python wrapper (`ctypes`) for use from other programs, and solver checks that stop a run with a clear error instead of returning wrong results (function domains, stalls). Planned next: Modelica packages and `extends`.

### Build and run
GenSAS builds with Bazel on **Linux, macOS and Windows** (MSYS2 MinGW-w64). See [Installation](doc/Installation.md) for prerequisites and per-OS instructions, and [Usage](doc/Usage.md) for command-line options and [job files](doc/JobFormat.md). [Supported Modelica](doc/Modelica.md) lists the language subset. The [C API](doc/CApi.md) and the [Python interface](doc/Python.md) run simulations from other programs, and the [server](doc/Server.md) offers them over HTTP ([API](api/openapi.yaml)); the [web client](doc/Web.md) is a browser front end for it. A [Docker image](doc/Deployment.md) packages both for a home network or Google Cloud.

### Companion Project
PowerSAS.m -- SAS-based power grid analysis toolbox for Matlab/GNU Octave users. 

### Acknowledgement
This work is supported by the Laboratory Directed Research and Development (LDRD) program of Argonne National Laboratory, provided by the U.S. Department of Energy Office of Science under Contract No. DE-AC02-06CH11357, and the U.S. Department of Energy Office of Electricity, Advanced Grid Modeling program under Grant DE-OE0000875. This acknowledgement refers to the work done at Argonne National Laboratory that this project builds on.

### Publications
* Rui Yao, Kai Sun, Feng Qiu, “Vectorized Efficient Computation of Padé Approximation for Semi-Analytical Simulation of Large-Scale Power Systems,” IEEE Transactions on Power Systems, 34 (5), 3957-3959, 2019.
* Rui Yao, Yang Liu, Kai Sun, Feng Qiu, Jianhui Wang,"Efficient and Robust Dynamic Simulation of Power Systems with Holomorphic Embedding", IEEE Transactions on Power Systems, 35 (2), 938 - 949, 2020.
* Rui Yao, Feng Qiu, "Novel AC Distribution Factor for Efficient Outage Analysis", IEEE Transactions on Power Systems, 35 (6), 4960-4963, 2020.
* Xin Xu, Rui Yao, Kai Sun, Feng Qiu, "A Semi-Analytical Solution Approach for Solving Constant-Coefficient First-Order Partial Differential Equations", IEEE Control Systems Letters, 6, 704-709, 2021.
* Rui Yao, Feng Qiu, Kai Sun, “Contingency Analysis Based on Partitioned and Parallel Holomorphic Embedding”, IEEE Transactions on Power Systems, in press.

### License
This software is under 3-clause BSD license. Please refer to LICENSE.md for details.