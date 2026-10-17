#!/usr/bin/env python3
"""Compares the GenSAS solver with the integrators of SciPy on a few small problems.

For every problem the script prints the number of solver steps and the largest error of the final state
(scaled by max(|reference|, 1)) for GenSAS at several error tolerances and for SciPy methods at several
relative tolerances. It counts steps, not seconds: a step of the series method is more expensive than a
step of a Runge-Kutta method, and Python and C++ timings are not comparable. See doc/SasComparison.md for
the results and what to read from them.

Usage (from the repository root, with the application built):

    python3 scripts/compare_with_scipy.py --app bazel-bin/app/app
    python3 scripts/compare_with_scipy.py --app bazel-bin/app/app --problems lorenz robertson

Needs numpy and scipy. GenSAS is run through its command line (`app -g`), with the verbose output used to
count the solver steps ("Stage start" lines) and the CSV output to read the final state.
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile

import numpy as np
from scipy.integrate import solve_ivp

ALL_METHODS = ("RK45", "DOP853", "LSODA", "BDF", "Radau")


class Problem:
    def __init__(self, title, model, f, y0, T, names, exact=None, ref_method="DOP853", methods=ALL_METHODS):
        self.title, self.model, self.f, self.y0, self.T, self.names = title, model, f, y0, T, names
        self.exact, self.ref_method, self.methods = exact, ref_method, methods

    def reference(self):
        if self.exact is not None:
            return np.array(self.exact, dtype=float)
        r = solve_ivp(self.f, (0, self.T), self.y0, method=self.ref_method, rtol=1e-13, atol=1e-16)
        return r.y[:, -1]


PROBLEMS = {
    "riccati": Problem(
        "Finite-time blow-up: x' = x^2, x(0) = 1, up to t = 0.99 (the solution is 1/(1-t), infinite at t = 1)",
        "model R\n  Real x(start = 1);\nequation\n  der(x) = x * x;\nend R;\n",
        lambda t, y: [y[0] ** 2], [1.0], 0.99, ["x"], exact=[1 / (1 - 0.99)],
    ),
    "runaway": Problem(
        "Thermal runaway: u' = exp(u), u(0) = 0, up to t = 0.999 (the solution is -ln(1-t), infinite at t = 1)",
        "model T\n  Real u(start = 0);\nequation\n  der(u) = exp(u);\nend T;\n",
        lambda t, y: [np.exp(y[0])], [0.0], 0.999, ["u"], exact=[-np.log(1 - 0.999)],
    ),
    "kepler": Problem(
        "Kepler orbit with eccentricity 0.9, ten periods (the exact end state is the start state)",
        "model K\n  Real x(start = 1.9);\n  Real y(start = 0);\n  Real vx(start = 0);\n"
        "  Real vy(start = 0.22941573387056177);\n  Real r(start = 1.9);\nequation\n  r = sqrt(x * x + y * y);\n"
        "  der(x) = vx;\n  der(y) = vy;\n  der(vx) = -x / (r * r * r);\n  der(vy) = -y / (r * r * r);\nend K;\n",
        lambda t, s: [s[2], s[3], -s[0] / (s[0] ** 2 + s[1] ** 2) ** 1.5, -s[1] / (s[0] ** 2 + s[1] ** 2) ** 1.5],
        [1.9, 0.0, 0.0, 0.22941573387056177], 20 * np.pi, ["x", "y", "vx", "vy"],
        exact=[1.9, 0.0, 0.0, 0.22941573387056177],
    ),
    "lorenz": Problem(
        "Lorenz system (chaotic), t = 0 to 10",
        "model Z\n  Real x(start = 1);\n  Real y(start = 1);\n  Real z(start = 1);\nequation\n"
        "  der(x) = 10 * (y - x);\n  der(y) = x * (28 - z) - y;\n  der(z) = x * y - 8.0 / 3 * z;\nend Z;\n",
        lambda t, y: [10 * (y[1] - y[0]), y[0] * (28 - y[2]) - y[1], y[0] * y[1] - 8.0 / 3 * y[2]],
        [1.0, 1.0, 1.0], 10.0, ["x", "y", "z"],
    ),
    "pendulum": Problem(
        "Pendulum released at 3.0 rad (near the unstable position), t = 0 to 100",
        "model P\n  Real th(start = 3.0);\n  Real om(start = 0);\nequation\n  der(th) = om;\n  der(om) = -sin(th);\nend P;\n",
        lambda t, y: [y[1], -np.sin(y[0])], [3.0, 0.0], 100.0, ["th", "om"],
    ),
    "lotka": Problem(
        "Predator and prey (Lotka-Volterra), t = 0 to 100",
        "model L\n  Real a(start = 10);\n  Real b(start = 10);\nequation\n  der(a) = 1.1 * a - 0.4 * a * b;\n"
        "  der(b) = 0.1 * a * b - 0.4 * b;\nend L;\n",
        lambda t, y: [1.1 * y[0] - 0.4 * y[0] * y[1], 0.1 * y[0] * y[1] - 0.4 * y[1]], [10.0, 10.0], 100.0, ["a", "b"],
    ),
    "vdp": Problem(
        "Van der Pol oscillator, mu = 100 (stiff between the jumps), t = 0 to 200",
        "model V\n  Real x(start = 2);\n  Real y(start = 0);\nequation\n  der(x) = y;\n"
        "  der(y) = 100 * (1 - x * x) * y - x;\nend V;\n",
        lambda t, y: [y[1], 100 * (1 - y[0] ** 2) * y[1] - y[0]], [2.0, 0.0], 200.0, ["x", "y"],
    ),
    "robertson": Problem(
        "Robertson chemical kinetics (very stiff), written as a DAE with y3 = 1 - y1 - y2, t = 0 to 1000",
        "model Rb\n  Real y1(start = 1);\n  Real y2(start = 0);\n  Real y3;\nequation\n"
        "  der(y1) = -0.04 * y1 + 1e4 * y2 * y3;\n  der(y2) = 0.04 * y1 - 1e4 * y2 * y3 - 3e7 * y2 * y2;\n"
        "  y1 + y2 + y3 = 1;\nend Rb;\n",
        lambda t, y: [-0.04 * y[0] + 1e4 * y[1] * (1 - y[0] - y[1]), 0.04 * y[0] - 1e4 * y[1] * (1 - y[0] - y[1]) - 3e7 * y[1] ** 2],
        [1.0, 0.0], 1000.0, ["y1", "y2"], ref_method="Radau", methods=("LSODA", "BDF", "Radau"),
    ),
}

EVENT_MODEL = (
    "model P\n  Real th(start = 3.0);\n  Real om(start = 0);\n  discrete Real n(start = 0);\nequation\n"
    "  der(th) = om;\n  der(om) = -sin(th);\n  when th < 0 then\n    n = pre(n) + 1;\n  end when;\nend P;\n"
)


def run_gensas(app, model, T, tol, extra=()):
    """Returns (steps, final values by name, event times, error text)."""
    with tempfile.TemporaryDirectory() as d:
        mo, out = os.path.join(d, "m.mo"), os.path.join(d, "out.csv")
        with open(mo, "w") as f:
            f.write(model)
        cmd = [app, "-g", "-m", "file", "-i", mo, "-t", str(T), "-s", "10", "-e", str(tol), "-k", str(min(1.0, T)), "-v", "-c", out, *extra]
        p = subprocess.run(cmd, capture_output=True, text=True)
        steps = p.stdout.count("Stage start")
        if not os.path.exists(out):
            return steps, {}, [], (p.stdout + p.stderr).strip().splitlines()[-1] if (p.stdout + p.stderr).strip() else "no output"
        with open(out, newline="") as f:
            rows = list(csv.reader(f))
        final = {name: float(v) for name, v in zip(rows[0], rows[-1])}
        events = []
        ev = out[:-4] + ".events.csv"
        if os.path.exists(ev):
            with open(ev, newline="") as f:
                events = [float(r[0]) for r in list(csv.reader(f))[1:]]
        return steps, final, events, ""


def scaled_error(values, ref):
    return max(abs(v - r) / max(abs(r), 1.0) for v, r in zip(values, ref))


def compare_problem(app, key, tolerances):
    p = PROBLEMS[key]
    ref = p.reference()
    print("\n### %s\n" % p.title)
    print("| method | steps | error |")
    print("|---|---:|---:|")
    for tol in tolerances:
        steps, final, _, err = run_gensas(app, p.model, p.T, tol)
        if err:
            print("| GenSAS, errorTol %g | %d | failed: %s |" % (tol, steps, err))
        else:
            print("| GenSAS, errorTol %g | %d | %.1e |" % (tol, steps, scaled_error([final[n] for n in p.names], ref)))
    for m in p.methods:
        for tol in tolerances:
            s = solve_ivp(p.f, (0, p.T), p.y0, method=m, rtol=tol, atol=tol * 1e-3)
            print("| %s, rtol %g | %d | %.1e |" % (m, tol, len(s.t) - 1, scaled_error(s.y[:, -1], ref)))


def compare_events(app, tolerances):
    T = 100.0
    f = lambda t, y: [y[1], -np.sin(y[0])]
    cross = lambda t, y: y[0]
    cross.direction = -1
    ref = solve_ivp(f, (0, T), [3.0, 0.0], method="DOP853", rtol=1e-13, atol=1e-16, events=cross).t_events[0]
    print("\n### Event times: the pendulum above passes th = 0 (downwards) %d times in 100 s\n" % len(ref))
    print("| method | steps | largest error of an event time (s) |")
    print("|---|---:|---:|")
    for tol in tolerances:
        steps, _, events, err = run_gensas(app, EVENT_MODEL, T, tol)
        bad = err or len(events) != len(ref)
        print("| GenSAS, errorTol %g | %d | %s |" % (tol, steps, "failed" if bad else "%.1e" % max(abs(a - b) for a, b in zip(events, ref))))
    for m in ("RK45", "DOP853", "LSODA", "Radau"):
        for tol in tolerances:
            s = solve_ivp(f, (0, T), [3.0, 0.0], method=m, rtol=tol, atol=tol * 1e-3, events=cross)
            t = s.t_events[0]
            bad = len(t) != len(ref)
            print("| %s, rtol %g | %d | %s |" % (m, tol, len(s.t) - 1, "failed" if bad else "%.1e" % max(abs(a - b) for a, b in zip(t, ref))))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--app", required=True, help="the GenSAS executable (bazel-bin/app/app)")
    ap.add_argument("--problems", nargs="*", default=[*PROBLEMS, "events"], help="problems to run (default: all): %s, events" % ", ".join(PROBLEMS))
    ap.add_argument("--tolerances", nargs="*", type=float, default=[1e-6, 1e-10])
    args = ap.parse_args()
    if not os.path.exists(args.app):
        sys.exit("cannot find the application: " + args.app)
    for key in args.problems:
        if key == "events":
            compare_events(args.app, args.tolerances)
        elif key in PROBLEMS:
            compare_problem(args.app, key, args.tolerances)
        else:
            sys.exit("unknown problem %s" % key)


if __name__ == "__main__":
    main()
