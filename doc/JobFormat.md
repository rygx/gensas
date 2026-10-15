# Job file format (version 1)

A job file describes one ModelicaSAS simulation as JSON. Run it with

```bash
bazel-bin/app/app run job.json [--set key=value ...] [--validate]
```

`--set solver.endTime=5` overrides one value (the value is JSON, so numbers and quoted strings both work; a value that is not valid JSON is taken as a string). `--validate` checks the job, prints the resolved options and exits without running. The exit code is 0 on success and 1 on any error.

```json
{
  "gensas": {"format": "job", "version": 1},
  "kind": "modelica",
  "model": {"file": "resources/mofile/test_bouncing_ball.mo"},
  "solver": {"endTime": 1.2, "segment": 1.0, "level": 15, "alphaTol": 1e-3, "errorTol": 1e-5},
  "output": {"step": 0.01, "json": "out/ball.json", "csv": "out/ball.csv", "mat": "out/ball.mat"}
}
```

| Key | Meaning | Default |
|---|---|---|
| `gensas.format`, `gensas.version` | must be `"job"` and `1` | required |
| `kind` | must be `"modelica"` (power flow is not available as a job yet) | required |
| `model.file` or `model.source` | a Modelica file, or the model text itself; exactly one of them | required |
| `model.piecewise` | `"events"` (exact) or `"smooth"` for `abs`, `min`, `max`, `sign` | `"events"` |
| `model.smoothEps` | smoothing parameter, positive | 1e-6 |
| `solver.endTime` | simulated time | 10 |
| `solver.segment` | length of one solver segment | 1.0 |
| `solver.level` | order of the semi-analytical solution, integer 1 to 1000 | 15 |
| `solver.alphaTol` | tolerance of the embedding variable | 1e-3 |
| `solver.errorTol` | error tolerance of the equations | 1e-5 |
| `output.step` | target spacing of the samples in the result; the last sample is always at the end time | 0.01 |
| `output.json`, `output.csv`, `output.mat` | result files to write (see [Result format](ResultFormat.md)); each is optional | none |

Rules:

* Relative paths (`model.file` and the output files) are resolved against the directory of the job file, not the current directory. Absolute paths are used as they are.
* Unknown keys, wrong types and unsupported versions or kinds are errors, and the message names the key (`unknown key solver.endTim`, `solver.endTime: expected a number`).
* Numeric values outside the supported ranges are not errors: they are moved to the nearest limit and a warning is printed (end time at most 1000 s, segment between 1e-4 and 10 s, error tolerance between 1e-10 and 1e-2, output step between 1e-3 and 1 s, alpha tolerance between 1e-9 and 1 times the segment).
* The flags of `app -g` build the same job, so both entry points behave the same way. In `-g` mode a missing or non-numeric flag value is an error.

A job file can make the program write files at the paths it names, so treat job files from untrusted sources like scripts.

The same structure is what the library takes (`ModelicaJob` and `Simulation` in `sas/Simulation.h`): fill it from a job file with `parseModelicaJob` or set the fields directly, call `Simulation::run` with an optional progress callback that can cancel the run, and read the sampled `SasResult`.
