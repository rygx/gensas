#!/usr/bin/env python3
"""Turns the coverage data written by scripts/coverage.sh into a Markdown summary.

Usage: coverage_summary.py DIR   (DIR holds cpp.json from gcovr and, optionally, python.json)
"""
import json
import os
import sys


# Source files that only power flow and the PSAT reader use. They are tested through Bazel only, so the Bazel-free
# coverage run cannot reach them. The summary reports C++ coverage with and without them, so that a coverage goal can
# follow the code that is being worked on. A new power-flow-only file has to be added here; when power flow is tested by
# this run, empty the list.
POWER_FLOW_FILES = (
    "io/CheIoUtil.",
    "io/CheDataFormat.",
    "io/MatPsatDataRW.",
    "util/CheCompUtil.",
    "util/CheState.",
)


def is_power_flow(name):
    return name.startswith(POWER_FLOW_FILES)


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


def sums(rows):
    return [sum(r[k] for r in rows) for k in (1, 2, 3, 4)]


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: coverage_summary.py DIR   (DIR is the output directory of scripts/coverage.sh)")
    out = sys.argv[1]
    print("## Test coverage\n")
    rows, total = cpp_table(os.path.join(out, "cpp.json"))
    print("C++ library (`sas/`, `io/`, `util/`, `capi/`), measured with the C++ unit tests and the Python and server tests.")
    print("Power flow and the PSAT reader are tested through Bazel only; their files are counted in the totals, and left out of the figures \"without power flow\".\n")
    print("| | Covered | Total | Percent |")
    print("|---|---:|---:|---:|")
    print("| C++ lines | %d | %d | %.1f%% |" % (total[0], total[1], pct(total[0], total[1])))
    print("| C++ branches | %d | %d | %.1f%% |" % (total[2], total[3], pct(total[2], total[3])))
    kept = sums([r for r in rows if not is_power_flow(r[0])])
    if kept != total:
        print("| **C++ lines, without power flow** | %d | %d | **%.1f%%** |" % (kept[0], kept[1], pct(kept[0], kept[1])))
        print("| C++ branches, without power flow | %d | %d | %.1f%% |" % (kept[2], kept[3], pct(kept[2], kept[3])))
    py = os.path.join(out, "python.json")
    if os.path.exists(py):
        with open(py) as f:
            t = json.load(f)["totals"]
        print("| Python wrapper lines | %d | %d | %.1f%% |" % (t["covered_lines"], t["num_statements"], t["percent_covered"]))
    srv = os.path.join(out, "server.json")
    if os.path.exists(srv):
        with open(srv) as f:
            t = json.load(f)["totals"]
        print("| Server lines (Python, without the worker processes) | %d | %d | %.1f%% |" % (t["covered_lines"], t["num_statements"], t["percent_covered"]))
    web = os.path.join(out, "web", "coverage-summary.json")
    if os.path.exists(web):
        with open(web) as f:
            t = json.load(f)["total"]
        print("| Web client lines (unit tests; the browser tests are not measured) | %d | %d | %.1f%% |" % (t["lines"]["covered"], t["lines"]["total"], t["lines"]["pct"]))
        print("| Web client branches | %d | %d | %.1f%% |" % (t["branches"]["covered"], t["branches"]["total"], t["branches"]["pct"]))
    print("\n<details><summary>C++ per file</summary>\n")
    print("| File | Lines | Line % | Branches | Branch % |")
    print("|---|---:|---:|---:|---:|")
    for name, c, n, bc, bn in sorted(rows, key=lambda r: pct(r[1], r[2])):
        tag = " *(power flow)*" if is_power_flow(name) else ""
        print("| `%s`%s | %d/%d | %.1f%% | %d/%d | %.1f%% |" % (name, tag, c, n, pct(c, n), bc, bn, pct(bc, bn)))
    print("\n</details>")


main()
