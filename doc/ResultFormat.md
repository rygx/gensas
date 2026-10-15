# Result format (version 1)

ModelicaSAS (`-g`) writes the simulation result sampled every `--step` seconds. All formats carry the same data; the format version is 1 and changes only when an existing field changes meaning.

The samples start at the start time, are `step` apart, and the last one is at the end time (closer to the previous one if the end time is not a multiple of `step`).

A hybrid model produces two samples at each event time: the first shows the values before the event, the second the values after it. Times are non-decreasing.

## JSON (`-j`)

One object:

| Key | Content |
|-----|---------|
| `gensas` | `{"format": "result", "version": 1}` |
| `time` | array of sample times |
| `<variable>` | one array per continuous variable, same length as `time`; `null` marks a value the solver did not provide (JSON has no NaN) |
| `discrete` | hybrid models only: object with one array per discrete variable, sampled like `time` |
| `events` | hybrid models with events only: array of `{"time", "assignments": {name: value}, "reinits": {name: value}}` |

A model variable must not be called `gensas`, `time`, `discrete` or `events`.

## CSV (`-c`)

Header `time,<continuous variables>,<discrete variables>`, then one row per sample, values with 17 significant digits. Fields containing a comma, quote or line break are quoted. Events are written to `<name>.events.csv` with the columns `time,kind,name,value`, where `kind` is `assign` or `reinit`.

## MAT (`-o`)

MATLAB level 5 file, readable by MATLAB, Octave, SciPy (`scipy.io.loadmat`) and matio; the writer is `io/Mat5Writer`, with no third-party code.

| Variable | Content |
|----------|---------|
| `gensas_version` | 1x1 format version |
| `t` | 1xN sample times |
| `s` | nVar x N continuous variables |
| `names` | nVar x L char matrix, one variable name per row, blank padded |
| `d`, `dnames` | hybrid models only: nDisc x N discrete variables and their names |
| `event_t` | hybrid models with events only: 1xE event times |

## In code

`buildResult(solutionSet, interval)` in `sas/SasResult.h` returns a `SasResult`; `makeResultWriter("json" | "csv")` and `makeMatResultWriter()` (`sas/SasResultMat.h`) write it. New formats implement `ResultWriter`.
