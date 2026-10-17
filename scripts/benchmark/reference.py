#!/usr/bin/env python3
"""High-precision reference solutions for the benchmark problems.

The benchmark scripts compare the final state of every solver with a reference. A reference computed with
a double precision integrator is only as good as that integrator: SciPy's DOP853 with rtol = 1e-13 ends the
pendulum problem with an error of 3e-10 (see doc/SasComparison.md), so errors below about 1e-9 would be
measured against noise. This module computes the references for the problems that have no closed form with
the Taylor series method of mpmath at 30 digits, an independent method in arithmetic that is far more
precise than the solvers that are compared. The results are kept in references.json (with the method and
an uncertainty), so that the benchmark does not have to compute them again.

    python3 scripts/benchmark/reference.py list
    python3 scripts/benchmark/reference.py compute lorenz pendulum lotka    # slow: minutes each
    python3 scripts/benchmark/reference.py compute ring --ring-sizes 10     # key ring10

Needs numpy and mpmath.
"""

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "references.json")
DIGITS = 30
# The Taylor series method of mpmath is run at 30 digits; the stated uncertainty is far larger than the
# method's own error and only has to be below anything a double precision solver can resolve.
MP_UNCERTAINTY = 1e-20


def load():
    if os.path.exists(CACHE):
        with open(CACHE, encoding="utf-8") as f:
            return json.load(f)
    return {}


def lookup(key):
    """(final state as floats, uncertainty, source text) of a cached reference, or None."""
    entry = load().get(key)
    if entry is None:
        return None
    return np.array([float(v) for v in entry["value"]]), entry["uncertainty"], entry["source"]


def problem_function(key, ring_n=0):
    """The right-hand side (for mpmath), the initial state and the end time of a problem.

    Constants that the models write as decimals are the double precision numbers that GenSAS reads, so that
    the reference solves exactly the problem the solvers solve."""
    from mpmath import mpf, sin

    if key == "lorenz":
        sigma, rho, beta = mpf(10), mpf(28), mpf(8.0 / 3)  # 8.0 / 3 in double precision, as the model reads it
        return (lambda t, y: [sigma * (y[1] - y[0]), y[0] * (rho - y[2]) - y[1], y[0] * y[1] - beta * y[2]]), [1, 1, 1], 10
    if key == "pendulum":
        return (lambda t, y: [y[1], -sin(y[0])]), [3, 0], 100
    if key == "lotka":
        a, b, c, d = mpf(1.1), mpf(0.4), mpf(0.1), mpf(0.4)
        return (lambda t, y: [a * y[0] - b * y[0] * y[1], c * y[0] * y[1] - d * y[1]]), [10, 10], 100
    if key == "vdp1":
        return (lambda t, y: [y[1], (1 - y[0] ** 2) * y[1] - y[0]]), [2, 0], 1000
    if key == "vdp100":
        mu = mpf(100)
        return (lambda t, y: [y[1], mu * (1 - y[0] ** 2) * y[1] - y[0]]), [2, 0], 200
    if key == "ring":
        n = ring_n
        P = [mpf(float(0.5 * np.sin(1.7 * (i + 1)))) for i in range(n)]
        damping = mpf(0.2)

        def f(t, y):
            th, om = y[:n], y[n:]
            return list(om) + [P[i] - damping * om[i] + sin(th[i - 1] - th[i]) + sin(th[(i + 1) % n] - th[i]) for i in range(n)]

        return f, [0] * (2 * n), 20
    raise KeyError(key)


def cache_key(key, ring_n=0):
    return "ring%d" % ring_n if key == "ring" else key


def compute(key, ring_n=0):
    from mpmath import mp, mpf, odefun, nstr

    mp.dps = DIGITS
    f, y0, T = problem_function(key, ring_n)
    sol = odefun(f, 0, [mpf(v) for v in y0])
    return [nstr(v, 25) for v in sol(T)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("action", choices=["list", "compute"])
    ap.add_argument("problems", nargs="*")
    ap.add_argument("--ring-sizes", nargs="*", type=int, default=[10])
    args = ap.parse_args()
    cache = load()
    if args.action == "list":
        for k, e in cache.items():
            print("%-8s uncertainty %.0e  %s  %s" % (k, e["uncertainty"], e["source"], " ".join(e["value"][:2]) + " ..."))
        return
    for key in args.problems:
        for n in args.ring_sizes if key == "ring" else [0]:
            name = cache_key(key, n)
            print("computing", name, flush=True)
            try:
                cache[name] = {
                    "value": compute(key, n),
                    "uncertainty": MP_UNCERTAINTY,
                    "source": "mpmath odefun (Taylor series), %d digits" % DIGITS,
                }
            except KeyError:
                sys.exit("no high-precision definition for " + key)
            with open(CACHE, "w", encoding="utf-8") as f:
                json.dump(cache, f, indent=1)
                f.write("\n")


if __name__ == "__main__":
    main()
