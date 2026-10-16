"""Reads the result files of a run and serves them as JSON or raw float64 columns."""

import json
import math
import os

import numpy as np


class ResultError(ValueError):
    """The request names a variable the run does not have."""


class Result:
    def __init__(self, directory):
        with open(os.path.join(directory, "meta.json"), encoding="utf-8") as f:
            meta = json.load(f)
        data = np.load(os.path.join(directory, "result.npz"))
        self.time = data["time"]
        self.names = meta["variables"]
        self.discrete_names = meta["discrete"]
        self.variables = {n: data["v%d" % i] for i, n in enumerate(self.names)}
        self.discrete = {n: data["d%d" % i] for i, n in enumerate(self.discrete_names)}
        self.events = meta["events"]
        self.termination = meta["termination"]

    def select(self, names):
        """The (continuous, discrete) variables to return: all of them, or those in `names`."""
        if not names:
            return dict(self.variables), dict(self.discrete)
        unknown = [n for n in names if n not in self.variables and n not in self.discrete]
        if unknown:
            raise ResultError("unknown variable(s): " + ", ".join(unknown))
        return {n: self.variables[n] for n in self.names if n in names}, {n: self.discrete[n] for n in self.discrete_names if n in names}


def reduce_indices(time, columns, max_points):
    """Indices of the samples to keep when there are more than `max_points`: the first and last sample, both
    samples of every event time, an even spread, and the extremes of the first few columns in each time bucket.
    The result has about `max_points` entries (a little more if there are many events)."""
    n = len(time)
    if max_points is None or n <= max_points:
        return np.arange(n)
    keep = {0, n - 1}
    same = np.nonzero(time[1:] == time[:-1])[0]  # two samples at one time: an event
    for i in same:
        keep.update((int(i), int(i) + 1))
    spread = max(2, max_points // 2)
    keep.update(int(i) for i in np.linspace(0, n - 1, spread).round())
    used = columns[:4]
    if used:
        buckets = max(1, (max_points - spread) // (2 * len(used)))
        edges = np.linspace(0, n, buckets + 1).astype(int)
        for lo, hi in zip(edges[:-1], edges[1:]):
            if hi <= lo:
                continue
            for col in used:
                part = col[lo:hi]
                if np.isfinite(part).any():
                    keep.add(lo + int(np.nanargmin(part)))
                    keep.add(lo + int(np.nanargmax(part)))
    return np.array(sorted(keep))


def _list(a):
    return [x if math.isfinite(x) else None for x in a.tolist()]  # JSON has no NaN or infinity


def _columns(result, names, max_points):
    variables, discrete = result.select(names)
    cols = list(variables.values()) + list(discrete.values())
    idx = reduce_indices(result.time, cols, max_points)
    return variables, discrete, idx


def to_json(result, names=None, max_points=None):
    variables, discrete, idx = _columns(result, names, max_points)
    return {
        "time": _list(result.time[idx]),
        "variables": {n: _list(a[idx]) for n, a in variables.items()},
        "discrete": {n: _list(a[idx]) for n, a in discrete.items()},
        "events": result.events,
        "terminated": result.termination,
    }


def to_f64(result, names=None, max_points=None):
    """Little-endian float64 columns: time, variables, discrete variables. Returns (bytes, samples)."""
    variables, discrete, idx = _columns(result, names, max_points)
    cols = [result.time[idx]] + [a[idx] for a in variables.values()] + [a[idx] for a in discrete.values()]
    return np.ascontiguousarray(np.stack(cols).astype("<f8")).tobytes(), len(idx)
