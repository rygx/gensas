# GenSAS server

HTTP API over the GenSAS Python interface. Start here: [doc/Server.md](../doc/Server.md). The contract is
[api/openapi.yaml](../api/openapi.yaml).

```bash
pip install -r server/requirements.txt
bazel build //capi:libgensas.so          # or: bash scripts/build_shared_lib.sh /tmp/gensas-lib
export GENSAS_LIBRARY=$PWD/bazel-bin/capi/libgensas.so
python -m gensas_server --port 8764      # run from the server/ directory
```

Tests: `pip install -r server/requirements-dev.txt` and `python -m pytest server/tests`.
