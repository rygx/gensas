# Server

`server/` is an HTTP server that checks Modelica models, runs simulations and serves the results. It is the back end of the web client and can also be used from scripts. It runs on your own machine or on a remote one; the client and the server are separate programs that only talk through the API in [api/openapi.yaml](../api/openapi.yaml).

## Running it

```bash
pip install -r server/requirements.txt
bash scripts/build_shared_lib.sh /tmp/gensas-lib            # or: bazel build //capi:libgensas.so
export GENSAS_LIBRARY=/tmp/gensas-lib/libgensas.so
cd server && python -m gensas_server --port 8765
curl http://127.0.0.1:8765/version
```

By default it listens on `127.0.0.1` only. To serve other machines use `--host 0.0.0.0` **and a token** (`--token SECRET` or `GENSAS_TOKEN`): every request except `GET /version` then needs `Authorization: Bearer SECRET`. Put the server behind a TLS proxy when it is reachable from the internet.

## How it works

* Every check and every run is a separate worker process (started with `spawn`). A crash of the solver, or a run that does not end, therefore never takes the server down: the server notices that the process is gone and marks the run `failed`, or kills it when it exceeds the time limit.
* `POST /runs` returns at once (`202`, state `queued`). At most `--max-concurrent-runs` run at a time, up to `--max-queued-runs` wait, and more are refused with `429`.
* Progress is available by polling `GET /runs/{id}` or as server-sent events from `GET /runs/{id}/events`. Browsers cannot set headers on `EventSource`, so a client that needs a token reads the stream with `fetch`.
* Results are stored in a temporary directory (`--work-dir`) and removed after `GENSAS_RETENTION` seconds or when more than `GENSAS_MAX_KEPT_RUNS` finished runs exist. `GET /runs/{id}/result` gives JSON (optionally only some variables and at most about `maxPoints` samples, which keep the first and last sample, event times and the extremes), CSV, MAT, or raw float64 columns (`format=f64`) for large results.
* The server never reads or writes files named by a client: models are sent as text, and the output is always taken from memory or from the run's own directory.

## Settings

| Flag | Variable | Default | Meaning |
|---|---|---|---|
| `--host`, `--port` | `GENSAS_HOST`, `GENSAS_PORT` | `127.0.0.1`, `8765` | where to listen |
| `--token` | `GENSAS_TOKEN` | none | bearer token required for all requests but `/version` |
| `--cors-origin` (repeatable) | `GENSAS_CORS_ORIGINS` (comma separated) | local origins when bound to loopback, else none | browser origins allowed to call the API; `*` allows all |
| `--static-dir` | `GENSAS_STATIC_DIR` | none | serve a built web client at `/` |
| `--max-concurrent-runs` | `GENSAS_MAX_CONCURRENT_RUNS` | 2 | simultaneous worker processes |
| `--max-queued-runs` | `GENSAS_MAX_QUEUED_RUNS` | 16 | waiting runs |
| `--run-timeout` | `GENSAS_RUN_TIMEOUT` | 120 | wall-clock seconds per run |
| | `GENSAS_CHECK_TIMEOUT` | 10 | seconds per check |
| `--max-source-bytes` | `GENSAS_MAX_SOURCE_BYTES` | 200000 | largest model text |
| | `GENSAS_MAX_KEPT_RUNS`, `GENSAS_RETENTION` | 50, 3600 | how many finished runs are kept and for how long |
| `--work-dir` | `GENSAS_WORK_DIR` | a temporary directory | where results are stored |

The engine itself also limits what a job can ask for (end time at most 1000 s and so on, see the [job format](JobFormat.md)); values outside the limits are moved to the nearest limit and reported as `warnings` of the run.

## Tests

`python -m pytest server/tests` (after `pip install -r server/requirements-dev.txt` and building the shared library). The contract test checks that the server implements exactly the paths of `api/openapi.yaml` and that the responses in the other tests match its schemas. The Server workflow runs them on every push.
