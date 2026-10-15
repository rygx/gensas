#!/usr/bin/env python3
"""Converts a PSAT .mat file (MATLAB v5, v7 or v7.3) into a GenSAS case file.

Usage:
    python3 scripts/psat_to_case.py [--json] INPUT.mat [OUTPUT]

The default output is INPUT with the suffix .case (binary: fast, for large cases);
--json writes the readable JSON encoding instead (suffix .case.json). Both encodings are
described in doc/CaseFormat.md.

Requires numpy and scipy; v7.3 files (HDF5) also need h5py.
"""
import argparse
import json
import struct
import sys

import numpy as np

# Keep in sync with kPsatMatrixNames in io/CaseFile.cpp.
PSAT_MATRICES = ["bus", "sw", "pv", "pq", "shunt", "line", "zip", "syn", "ind", "tg", "exc"]


def load_matrices(path):
    """Returns {name: 2-D float64 array in MATLAB shape (rows, cols)}."""
    out = {}
    try:
        from scipy.io import loadmat

        data = loadmat(path)
        for name in PSAT_MATRICES:
            value = data.get(name)
            if isinstance(value, np.ndarray) and value.ndim == 2 and value.dtype.kind in "fiub":
                out[name] = np.asarray(value, dtype="<f8")
    except NotImplementedError:
        # MATLAB v7.3 files are HDF5 files.
        import h5py

        with h5py.File(path, "r") as h5:
            for name in PSAT_MATRICES:
                if name in h5 and isinstance(h5[name], h5py.Dataset) and h5[name].ndim == 2:
                    # MATLAB stores column-major; HDF5 shows the dimensions reversed.
                    out[name] = np.asarray(h5[name][()], dtype="<f8").T
    return out


def write_binary(path, matrices):
    names = sorted(matrices)
    entries, offset = [], 0
    for name in names:
        rows, cols = matrices[name].shape
        entries.append({"name": name, "rows": rows, "cols": cols, "offset": offset})
        offset += rows * cols * 8
    header = json.dumps(
        {"gensas": {"format": "case", "version": 1}, "matrices": entries}, separators=(",", ":")
    ).encode()
    used = 16 + len(header)
    with open(path, "wb") as f:
        f.write(b"GSCASE01")
        f.write(struct.pack("<Q", len(header)))
        f.write(header)
        f.write(b"\0" * (-used % 8))
        for name in names:
            f.write(np.asfortranarray(matrices[name], dtype="<f8").tobytes(order="F"))


def write_json(path, matrices):
    doc = {
        "gensas": {"format": "case", "version": 1},
        "matrices": {name: m.tolist() for name, m in sorted(matrices.items())},
    }
    with open(path, "w") as f:
        json.dump(doc, f, separators=(",", ":"))
        f.write("\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input")
    ap.add_argument("output", nargs="?")
    ap.add_argument("--json", action="store_true", help="write the JSON encoding")
    args = ap.parse_args()

    matrices = load_matrices(args.input)
    if "bus" not in matrices:
        sys.exit(f"{args.input}: no 'bus' matrix; is this a PSAT data file?")
    stem = args.input[:-4] if args.input.endswith(".mat") else args.input
    output = args.output or stem + (".case.json" if args.json else ".case")
    (write_json if args.json else write_binary)(output, matrices)
    print(f"{output}: " + ", ".join(f"{n} {m.shape[0]}x{m.shape[1]}" for n, m in sorted(matrices.items())))


if __name__ == "__main__":
    main()
