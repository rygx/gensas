#!/usr/bin/env python3
"""Why GenSAS needs many steps on some stiff problems (see doc/SasComparison.md, "Stiff problems").

Three experiments, each printing a markdown table:

  steps    the steps of GenSAS on the Robertson problem against the fast eigenvalue of the Jacobian:
           the step length is a constant times 1/|lambda| (needs --app)
  radius   the radius of convergence of the Taylor series of the exact solution about a point of the Robertson
           trajectory, for states moved off the slow manifold by delta: it is about ln(1/delta)/|lambda|
  models   steps of GenSAS on small stiff models, linear and nonlinear in the fast variable, started on and
           off the slow manifold, for growing stiffness L (needs --app)
  first    the first step from a state at a distance delta from the slow manifold of the nonlinear model
           (needs --app)

    python3 scripts/benchmark/stiffness_diagnosis.py --app bazel-bin/app/app
    python3 scripts/benchmark/stiffness_diagnosis.py --app bazel-bin/app/app radius

Needs numpy, scipy and mpmath.
"""

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile

import numpy as np
from scipy.integrate import solve_ivp

try:
    from mpmath import mp, mpf
except ImportError:
    sys.exit("this script needs mpmath (pip install mpmath)")

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import compare_with_scipy as cmp  # noqa: E402


def robertson_rhs(t, y):
    y1, y2 = y
    y3 = 1 - y1 - y2
    return [-0.04 * y1 + 1e4 * y2 * y3, 0.04 * y1 - 1e4 * y2 * y3 - 3e7 * y2 * y2]


def fast_eigenvalue(y1, y2):
    """The eigenvalue of the Jacobian of the Robertson problem (with y3 = 1 - y1 - y2) of the largest magnitude."""
    y3 = 1 - y1 - y2
    J = np.array([[-0.04 - 1e4 * y2, 1e4 * y3 - 1e4 * y2], [0.04 + 1e4 * y2, -1e4 * y3 + 1e4 * y2 - 6e7 * y2]])
    e = np.linalg.eigvals(J)
    return abs(e[np.argmax(abs(e))])


def run_gensas(app, model, T, tol, extra=()):
    """(list of (start, end) of the steps, final values by name, error text)."""
    with tempfile.TemporaryDirectory() as d:
        mo, out = os.path.join(d, "m.mo"), os.path.join(d, "o.csv")
        with open(mo, "w", encoding="utf-8") as f:
            f.write(model)
        cmd = [app, "-g", "-m", "file", "-i", mo, "-t", str(T), "-s", "10", "-e", str(tol), "-k", "1", "-v", "-c", out, *extra]
        p = subprocess.run(cmd, capture_output=True, text=True)
        steps = [(float(a), float(b)) for a, b in re.findall(r"Stage start ([\d.e+-]+) end ([\d.e+-]+), length", p.stdout)]
        if p.returncode != 0 or not os.path.exists(out):
            text = (p.stdout + p.stderr).strip().splitlines()
            return steps, {}, text[-1] if text else "failed"
        with open(out, newline="", encoding="utf-8") as f:
            rows = list(csv.reader(f))
        if len(rows) < 2:
            return steps, {}, "the output file has no data rows"
        return steps, {n: float(v) for n, v in zip(rows[0], rows[-1])}, ""


def experiment_steps(app, T=100.0, tol=1e-6):
    steps, _, err = run_gensas(app, cmp.PROBLEMS["robertson"].model, T, tol)
    if err:
        sys.exit(err)
    ref = solve_ivp(robertson_rhs, (0, T), [1, 0], method="Radau", rtol=1e-11, atol=1e-14, dense_output=True)
    print("Robertson, t = 0 to %g, errorTol %g: %d steps.\n" % (T, tol, len(steps)))
    print("| t | step length h | fast eigenvalue | h · abs(lambda) |")
    print("|---:|---:|---:|---:|")
    pick = [s for s in steps if s[0] > 0.01]
    for a, b in pick[:: max(1, len(pick) // 12)]:
        lam = fast_eigenvalue(*ref.sol((a + b) / 2))
        print("| %.3g | %.2e | %.0f | %.1f |" % (a, b - a, lam, (b - a) * lam))


def taylor_y2(y1, y2, terms=120):
    """Taylor coefficients of y2 for the Robertson problem about a state, in high precision."""
    mp.dps = 60
    a1, a2 = [mpf(y1)], [mpf(y2)]
    for k in range(terms):
        y3 = [(1 if i == 0 else 0) - a1[i] - a2[i] for i in range(k + 1)]
        p23 = sum(a2[i] * y3[k - i] for i in range(k + 1))
        p22 = sum(a2[i] * a2[k - i] for i in range(k + 1))
        a1.append((-mpf("0.04") * a1[k] + 10**4 * p23) / (k + 1))
        a2.append((mpf("0.04") * a1[k] - 10**4 * p23 - 3 * 10**7 * p22) / (k + 1))
    return a2


def experiment_radius(t=5.0):
    r = solve_ivp(robertson_rhs, (0, t), [1, 0], method="Radau", rtol=1e-13, atol=1e-17)
    y1, y2 = r.y[:, -1]
    lam = fast_eigenvalue(y1, y2)
    print("Robertson at t = %g: y1 = %.9g, y2 = %.6e, fast eigenvalue %.0f. Radius of convergence rho of the Taylor series of y2 (root test on the" % (t, y1, y2, lam))
    print("coefficients 60 to 120), for y2 moved by delta, relative to y2, from the state computed by Radau.\n")
    print("| delta / y2 | rho | rho · abs(lambda) | ln(1 / (delta / y2)) |")
    print("|---:|---:|---:|---:|")
    for rel in (0, 1e-3, 1e-5, 1e-7, 1e-9, 1e-11):
        a = taylor_y2(y1, y2 * (1 + rel))
        roots = [abs(a[k]) ** (mpf(1) / k) for k in range(60, 120) if a[k] != 0]
        if not roots:
            sys.exit("all Taylor coefficients 60 to 120 are zero")
        rho = 1.0 / float(max(roots))
        print("| %s | %.3e | %.1f | %s |" % ("0" if rel == 0 else "%.0e" % rel, rho, rho * lam, "" if rel == 0 else "%.1f" % np.log(1 / rel)))


def experiment_models(app, T=10.0, tol=1e-8):
    start = {"on": "2", "off": "2.1"}
    lin = "der(y) = %.1f * (y - s) + c;"
    nl = "der(y) = %.1f * (y * y - s * s) + c;"
    template = "model P\n  Real s(start = 2);\n  Real c(start = 1);\n  Real y(start = %s);\nequation\n  der(s) = c;\n  der(c) = -(s - 2);\n  %s\nend P;\n"
    Ls = (-1e1, -1e2, -1e3, -1e4, -1e5)
    print("y relaxes to the slow solution s with rate about |L| (errorTol %g, alphaTol 1e-7, t = 0 to %g). Steps of GenSAS:\n" % (tol, T))
    print("| L | linear, on the manifold | nonlinear, on the manifold | linear, started off | nonlinear, started off |")
    print("|---:|---:|---:|---:|---:|")
    for L in Ls:
        cells = []
        for where, eq in (("on", lin), ("on", nl), ("off", lin), ("off", nl)):
            steps, _, err = run_gensas(app, template % (start[where], eq % L), T, tol, ("-a", "1e-7"))
            cells.append("failed" if err else str(len(steps)))
        print("| %g | %s |" % (L, " | ".join(cells)))


def experiment_first(app, L=-1e4):
    """The first step from a state at a distance delta from the slow manifold of the nonlinear model."""
    template = "model P\n  Real s(start = 2);\n  Real c(start = 1);\n  Real y(start = %.17g);\nequation\n  der(s) = c;\n  der(c) = -(s - 2);\n  der(y) = %.1f * (y * y - s * s) + c;\nend P;\n"
    lam = 4 * abs(L)  # 2 s |L| with s = 2
    print("Nonlinear model, L = %g (fast rate %g). The first step of GenSAS from y(0) = s(0) + delta (errorTol 1e-8, alphaTol 1e-9):\n" % (L, lam))
    print("| delta | ln(1 + 2s/delta) | first step · abs(lambda) | ratio |")
    print("|---:|---:|---:|---:|")
    for delta in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8):
        steps, _, err = run_gensas(app, template % (2 + delta, L), 0.5, 1e-8, ("-a", "1e-9"))
        if err or not steps:
            print("| %g | failed |" % delta)
            continue
        h = (steps[0][1] - steps[0][0]) * lam
        pole = np.log(1 + 4 / delta)
        print("| %g | %.1f | %.1f | %.2f |" % (delta, pole, h, h / pole))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--app", help="the GenSAS executable")
    ap.add_argument("experiments", nargs="*", choices=["steps", "radius", "models", "first"], default=["steps", "radius", "models", "first"])
    args = ap.parse_args()
    if args.app and not os.path.exists(args.app) and os.path.exists(args.app + ".exe"):
        args.app += ".exe"
    for name in args.experiments:
        print("\n### %s\n" % name)
        if name == "radius":
            experiment_radius()
            continue
        if not args.app or not os.path.exists(args.app):
            sys.exit("the experiment %s needs --app" % name)
        {"steps": experiment_steps, "models": experiment_models, "first": experiment_first}[name](args.app)


if __name__ == "__main__":
    main()
