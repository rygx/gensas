#!/usr/bin/env python3
"""Turns the coverage data written by scripts/coverage.sh into a Markdown summary.

Usage: coverage_summary.py DIR   (DIR holds cpp.json from gcovr and, optionally, python.json)
"""
import json
import os
import sys


def pct(covered, total):
    return 100.0 * covered / total if total else 100.0


def cpp_table(path):
    with open(path) as f:
        data = json.load(f)
    rows = []
    total = [0, 0, 0, 0]
    for entry in sorted(data["files"], key=lambda e: e["file"]):
        lines = [ln for ln in entry["lines"] if not ln.get("gcovr/noncode")]
        covered = sum(1 for ln in lines if ln["count"] > 0)
        branches = [b for ln in lines for b in (ln.get("branches") or [])]
        branch_covered = sum(1 for b in branches if b["count"] > 0)
        rows.append((entry["file"], covered, len(lines), branch_covered, len(branches)))
        for k, v in enumerate((covered, len(lines), branch_covered, len(branches))):
            total[k] += v
    return rows, total


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: coverage_summary.py DIR   (DIR is the output directory of scripts/coverage.sh)")
    out = sys.argv[1]
    print("## Test coverage\n")
    rows, total = cpp_table(os.path.join(out, "cpp.json"))
    print("C++ library (`sas/`, `io/`, `util/`, `capi/`), measured with the C++ unit tests and the Python tests.")
    print("Power flow and the PSAT reader are tested through Bazel only and are not included.\n")
    print("| | Covered | Total | Percent |")
    print("|---|---:|---:|---:|")
    print("| Lines | %d | %d | %.1f%% |" % (total[0], total[1], pct(total[0], total[1])))
    print("| Branches | %d | %d | %.1f%% |" % (total[2], total[3], pct(total[2], total[3])))
    py = os.path.join(out, "python.json")
    if os.path.exists(py):
        with open(py) as f:
            t = json.load(f)["totals"]
        print("| Python wrapper lines | %d | %d | %.1f%% |" % (t["covered_lines"], t["num_statements"], t["percent_covered"]))
    print("\n<details><summary>C++ per file</summary>\n")
    print("| File | Lines | Line % | Branches | Branch % |")
    print("|---|---:|---:|---:|---:|")
    for name, c, n, bc, bn in sorted(rows, key=lambda r: pct(r[1], r[2])):
        print("| `%s` | %d/%d | %.1f%% | %d/%d | %.1f%% |" % (name, c, n, pct(c, n), bc, bn, pct(bc, bn)))
    print("\n</details>")


main()
