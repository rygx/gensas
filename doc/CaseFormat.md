# GenSAS case file (version 1)

The power-flow mode (`-p`) reads its data from a case file: the PSAT data matrices `bus`, `sw`, `pv`, `pq`, `shunt`, `line`, `zip`, `syn`, `ind`, `tg` and `exc` (only `bus` is mandatory). Each matrix is a real double matrix with the PSAT column layout.

Convert a PSAT `.mat` file (MATLAB v5, v7 or v7.3) with

```
python3 scripts/psat_to_case.py d_70k.mat            # writes d_70k.case
python3 scripts/psat_to_case.py --json d_014.mat     # writes d_014.case.json
```

The script needs numpy and scipy, and h5py for v7.3 files. Loading a case file takes a few milliseconds even for the 210k-bus case (about 70 MB); the same data takes several times longer from `.mat`.

## Binary encoding (`.case`)

| Bytes | Content |
|-------|---------|
| 0-7 | `GSCASE01` |
| 8-15 | length H of the header, unsigned 64-bit little endian |
| 16 to 16+H | JSON header (below) |
| up to the next multiple of 8 | zero padding; the data start here |
| data | the matrices as little-endian float64, column by column |

Header:

```json
{"gensas": {"format": "case", "version": 1},
 "matrices": [{"name": "bus", "rows": 14, "cols": 6, "offset": 0}, ...]}
```

`offset` is in bytes from the start of the data.

## JSON encoding (`.case.json`, for small cases)

```json
{"gensas": {"format": "case", "version": 1},
 "matrices": {"bus": [[1, 138, 1.06, 0, 1, 1], ...], "pq": [[...]]}}
```

Each matrix is an array of rows. Both encodings are recognized automatically by `readCaseFile` (`io/CaseFile.h`), which reports a damaged or truncated file or an unknown version with an error message.

## Notes

The solver reads no `.mat` files: it has no matio, HDF5 or zlib dependency. MATLAB, Octave and SciPy read the `.mat` result files it writes (level 5 format).
