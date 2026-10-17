#!/usr/bin/env python3
"""Speed and accuracy of GenSAS against compiled integrators.

compare_with_scipy.py counts steps; this script also measures run time, and does so fairly: the baselines
(scripts/benchmark/baselines.cpp) are compiled C++ with a compiled right-hand side, so GenSAS is not compared
with Python call overhead. They are

  * RK45  - Dormand-Prince 5(4), the same method and step control as SciPy's RK45,
  * DOP853 - the eighth order method of SciPy's DOP853,
  * CVODE - BDF with a dense Newton iteration from SUNDIALS (what many Modelica tools use for ODEs);
            only if SUNDIALS is given with --sundials.

For every problem the table shows steps, the error of the final state against a reference, and the run time
(one run, in milliseconds; for GenSAS the "Computation time" of the application, best of three, which includes
building the equations but not reading the file). The baselines only compute the final state, GenSAS also
writes a few samples of the result.

Usage (from the repository root, application built):

    python3 scripts/benchmark_solvers.py --app bazel-bin/app/app
    python3 scripts/benchmark_solvers.py --app bazel-bin/app/app --sundials /usr   # with CVODE
    python3 scripts/benchmark_solvers.py --app bazel-bin/app/app --problems lorenz ring --ring-sizes 10 40

Needs numpy, scipy (for the reference solutions) and a C++ compiler. --sundials takes the directory that has
include/ and lib/ of SUNDIALS 5 or later (the libraries are linked statically if the .a files exist).
"""

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_with_scipy as cmp  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def find_sundials(root):
    """The include directory and the library directory of a SUNDIALS installation below `root`."""
    include = libdir = None
    for dirpath, dirs, files in os.walk(root):
        if include is None and "cvode" in dirs and os.path.exists(os.path.join(dirpath, "cvode", "cvode.h")):
            include = dirpath
        if libdir is None and ("libsundials_cvode.a" in files or "libsundials_cvode.so" in files):
            libdir = dirpath
    if include is None or libdir is None:
        sys.exit("no SUNDIALS (cvode/cvode.h and libsundials_cvode) below " + root)
    return include, libdir


def build_baselines(outdir, sundials):
    exe = os.path.join(outdir, "baselines")
    cmd = [os.environ.get("CXX", "c++"), "-std=c++17", "-O2", os.path.join(HERE, "benchmark", "baselines.cpp"), "-I" + os.path.join(HERE, "benchmark"), "-o", exe]
    if sundials:
        include, libdir = find_sundials(sundials)
        cmd += ["-DWITH_CVODE", "-I" + include]
        names = ["cvode", "nvecserial", "sunlinsoldense", "sunmatrixdense"]
        static = [os.path.join(libdir, "libsundials_%s.a" % n) for n in names]
        if all(os.path.exists(f) for f in static):
            cmd += static
        else:
            cmd += ["-L" + libdir, "-Wl,-rpath," + libdir] + ["-lsundials_" + n for n in names]
        cmd += ["-lm"]
    subprocess.run(cmd, check=True)
    return exe


def run_baseline(exe, method, key, rtol, ring_n):
    p = subprocess.run([exe, method, key, repr(rtol), str(ring_n)], capture_output=True, text=True)
    out = p.stdout.split()
    if p.returncode != 0 or not out or out[0] == "failed":
        return None
    steps, evals, seconds, n = int(out[0]), int(out[1]), float(out[2]), int(out[3])
    return steps, evals, seconds, np.array([float(v) for v in out[4 : 4 + n]])


def run_gensas(app, model, T, tol):
    """Returns (steps, final values by column name, best computation time in seconds) or an error text."""
    best, steps, final, err = None, 0, {}, ""
    with tempfile.TemporaryDirectory() as d:
        mo, out = os.path.join(d, "m.mo"), os.path.join(d, "out.csv")
        with open(mo, "w", encoding="utf-8") as f:
            f.write(model)
        for rep in range(3):
            cmd = [app, "-g", "-m", "file", "-i", mo, "-t", str(T), "-s", "10", "-e", str(tol), "-k", str(min(1.0, T)), "-c", out]
            if rep == 0:
                cmd.append("-v")
            p = subprocess.run(cmd, capture_output=True, text=True)
            text = p.stdout + p.stderr
            if rep == 0:
                steps = p.stdout.count("Stage start")
            m = re.search(r"Computation time: ([\d.e+-]+)", text)
            if p.returncode != 0 or not m or not os.path.exists(out):
                lines = text.strip().splitlines()
                return steps, {}, None, lines[-1] if lines else "failed"
            t = float(m.group(1))
            best = t if best is None else min(best, t)
        with open(out, newline="", encoding="utf-8") as f:
            rows = list(csv.reader(f))
        final = {name: float(v) for name, v in zip(rows[0], rows[-1])}
    return steps, final, best, err


def fmt_time(seconds):
    ms = seconds * 1e3
    return "%.3f" % ms if ms < 1 else ("%.1f" % ms if ms < 100 else "%.0f" % ms)


def benchmark(app, exe, key, problem, tolerances, ring_n, methods):
    ref = problem.reference()

    def err(v):
        return cmp.scaled_error(list(v), ref)

    print("\n### %s\n" % problem.title)
    print("| method | steps | error | time (ms) |")
    print("|---|---:|---:|---:|")
    for tol in tolerances:
        steps, final, t, bad = run_gensas(app, problem.model, problem.T, tol)
        if bad:
            print("| GenSAS, errorTol %g | %d | failed: %s | |" % (tol, steps, bad))
        else:
            print("| GenSAS, errorTol %g | %d | %.1e | %s |" % (tol, steps, err([final[n] for n in problem.names]), fmt_time(t)))
    for m in methods:
        for tol in tolerances:
            r = run_baseline(exe, m, key, tol, ring_n)
            label = {"rk45": "RK45", "dop853": "DOP853", "cvode": "CVODE (BDF)"}[m]
            if r is None:
                print("| %s, rtol %g | | failed | |" % (label, tol))
            else:
                print("| %s, rtol %g | %d | %.1e | %s |" % (label, tol, r[0], err(r[3]), fmt_time(r[2])))


LADDER = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13]


def time_to_accuracy(app, exe, key, problem, targets, ring_n, methods):
    """For each target error, the cheapest tolerance of each method that reaches it (a ladder of tolerances is tried)."""
    ref = problem.reference()

    def err(v):
        return cmp.scaled_error(list(v), ref)

    rows = {}
    for target in targets:
        row = {}
        for name in ["gensas"] + methods:
            found = None
            for tol in LADDER:
                if name == "gensas":
                    if tol < 1e-10 or tol > 1e-6:  # the range GenSAS accepts for errorTol
                        continue
                    steps, final, t, bad = run_gensas(app, problem.model, problem.T, tol)
                    if bad:
                        continue
                    e, st, sec = err([final[n] for n in problem.names]), steps, t
                else:
                    r = run_baseline(exe, name, key, tol, ring_n)
                    if r is None:
                        continue
                    e, st, sec = err(r[3]), r[0], r[2]
                if e <= target:
                    found = (st, sec)
                    break
            row[name] = found
        rows[target] = row
    return rows


def print_accuracy_table(title, results, methods):
    names = ["gensas"] + methods
    labels = {"gensas": "GenSAS", "rk45": "RK45", "dop853": "DOP853", "cvode": "CVODE (BDF)"}
    print("\n### %s\n" % title)
    print("| error at most | " + " | ".join(labels[n] + " steps, ms" for n in names) + " |")
    print("|---|" + "---:|" * len(names))
    for target, row in results.items():
        cells = []
        for n in names:
            f = row[n]
            cells.append("not reached" if f is None else "%d, %s" % (f[0], fmt_time(f[1])))
        print("| %g | %s |" % (target, " | ".join(cells)))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--app", required=True, help="the GenSAS executable (bazel-bin/app/app)")
    ap.add_argument("--sundials", help="SUNDIALS installation (the directory with include/ and lib/); enables CVODE")
    keys = [k for k in cmp.PROBLEMS if k not in ("ring",)] + ["ring"]
    ap.add_argument("--problems", nargs="*", default=[k for k in keys if k != "robertson"] + ["robertson"], help="default: all of " + ", ".join(keys))
    ap.add_argument("--tolerances", nargs="*", type=float, default=[1e-6, 1e-10])
    ap.add_argument("--ring-sizes", nargs="*", type=int, default=[20], help="machines in the ring problem")
    ap.add_argument("--targets", nargs="*", type=float, help="instead of fixed tolerances, report the steps and time each method needs to reach these errors")
    args = ap.parse_args()
    if not os.path.exists(args.app) and os.path.exists(args.app + ".exe"):
        args.app += ".exe"
    if not os.path.exists(args.app):
        sys.exit("cannot find the application: " + args.app)
    with tempfile.TemporaryDirectory() as d:
        exe = build_baselines(d, args.sundials)
        methods = ["rk45", "dop853"] + (["cvode"] if args.sundials else [])
        for key in args.problems:
            if key == "ring":
                for n in args.ring_sizes:
                    prob = cmp.ring_problem(n)
                    if args.targets:
                        print_accuracy_table(prob.title, time_to_accuracy(args.app, exe, "ring", prob, args.targets, n, methods), methods)
                    else:
                        benchmark(args.app, exe, "ring", prob, args.tolerances, n, methods)
            elif key in cmp.PROBLEMS:
                prob = cmp.PROBLEMS[key]
                if args.targets:
                    print_accuracy_table(prob.title, time_to_accuracy(args.app, exe, key, prob, args.targets, 0, methods), methods)
                else:
                    benchmark(args.app, exe, key, prob, args.tolerances, 0, methods)
            else:
                sys.exit("unknown problem " + key)


if __name__ == "__main__":
    main()
