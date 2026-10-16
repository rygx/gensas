import json
import os
import time

import numpy as np
import pytest

from conftest import BALL, DECAY, make_client, start, wait_for
from contract import validate
from gensas_server.examples import EXAMPLES


def test_version(client):
    r = client.get("/version")
    assert r.status_code == 200
    validate("Version", r.json())
    assert r.json()["api"] == 1 and r.json()["engine"]["apiVersion"] >= 2


def test_examples_have_the_documented_shape_and_run(client):
    r = client.get("/examples")
    assert r.status_code == 200
    for item in r.json():
        validate("Example", item)
    assert [e["id"] for e in r.json()] == [e["id"] for e in EXAMPLES]


@pytest.mark.parametrize("example", EXAMPLES, ids=[e["id"] for e in EXAMPLES])
def test_every_example_runs(client, example):
    checked = client.post("/check", json={"source": example["source"]}).json()
    assert checked["ok"], checked
    run = start(client, example["source"], example["endTime"])
    done = wait_for(client, run["id"])
    assert done["state"] in ("completed", "terminated"), done


def test_check_good_model_returns_info(client):
    r = client.post("/check", json={"source": DECAY, "parameters": {"k": 3}})
    assert r.status_code == 200
    body = r.json()
    validate("CheckResponse", body)
    assert body["ok"] and body["diagnostics"] == []
    assert body["info"]["name"] == "D"
    assert body["info"]["parameters"][0]["value"] == 3.0
    assert [(v["name"], v["kind"]) for v in body["info"]["variables"]] == [("x", "state")]


def test_check_reports_errors_with_lines(client):
    r = client.post("/check", json={"source": "model M\n  Real x(start = 1);\nequation\n  der(x) = -y;\nend M;\n"})
    body = r.json()
    validate("CheckResponse", body)
    assert not body["ok"]
    assert body["diagnostics"][0]["line"] == 4 and "undeclared identifier 'y'" in body["diagnostics"][0]["message"]


def test_check_bad_override_is_a_diagnostic(client):
    body = client.post("/check", json={"source": DECAY, "parameters": {"nope": 1}}).json()
    validate("CheckResponse", body)
    assert not body["ok"] and "parameters.nope" in body["diagnostics"][0]["message"]


def test_invalid_requests(client):
    for path, body in [("/check", {}), ("/check", {"source": 5}), ("/check", {"source": DECAY, "bogus": 1}),
                       ("/runs", {"model": {"source": DECAY}, "solver": {"endTime": "long"}}), ("/runs", {"model": {"source": DECAY, "piecewise": "x"}})]:
        r = client.post(path, json=body)
        assert r.status_code == 400, (path, body, r.text)
        validate("Error", r.json())
        assert r.json()["error"]["code"] == "invalid_request"


def test_too_large_source():
    with make_client(max_source_bytes=100) as c:
        for path, body in [("/check", {"source": "x" * 101}), ("/runs", {"model": {"source": "x" * 101}})]:
            r = c.post(path, json=body)
            assert r.status_code == 413
            validate("Error", r.json())


def test_run_lifecycle_and_json_result(client):
    run = start(client, DECAY, 2.0, parameters={"k": 2})
    validate("Run", run)
    assert run["state"] in ("queued", "running")
    done = wait_for(client, run["id"])
    validate("Run", done)
    assert done["state"] == "completed", done
    assert done["summary"]["variables"] == ["x"] and done["summary"]["endTime"] == pytest.approx(2.0)
    assert done["progress"]["time"] == pytest.approx(2.0)
    assert any(r["id"] == run["id"] for r in client.get("/runs").json())
    result = client.get("/runs/%s/result" % run["id"]).json()
    validate("ResultJson", result)
    assert len(result["time"]) == done["summary"]["samples"] == len(result["variables"]["x"])
    assert result["variables"]["x"][-1] == pytest.approx(np.exp(-4.0), abs=1e-4)
    assert result["events"] == [] and result["terminated"] is None


def test_hybrid_result_has_events_and_discrete_samples(client):
    run = start(client, BALL, 1.2)
    done = wait_for(client, run["id"])
    assert done["state"] == "completed", done
    assert done["summary"]["events"] >= 1
    result = client.get("/runs/%s/result" % run["id"]).json()
    validate("ResultJson", result)
    assert result["events"][0]["reinits"]


def test_result_formats(client):
    run = wait_for(client, start(client, DECAY, 1.0)["id"])
    rid = run["id"]
    csv = client.get("/runs/%s/result?format=csv" % rid)
    assert csv.status_code == 200 and csv.text.splitlines()[0] == "time,x"
    mat = client.get("/runs/%s/result?format=mat" % rid)
    assert mat.status_code == 200 and mat.content[:6] == b"MATLAB"
    raw = client.get("/runs/%s/result?format=f64" % rid)
    n = int(raw.headers["x-gensas-samples"])
    assert n == run["summary"]["samples"]
    matrix = np.frombuffer(raw.content, dtype="<f8").reshape(2, n)
    assert matrix[0][-1] == pytest.approx(1.0) and matrix[1][-1] == pytest.approx(np.exp(-1.0), abs=1e-4)
    assert client.get("/runs/%s/result?format=xml" % rid).status_code == 400


def test_variable_selection_and_downsampling(client):
    source = "model P\n  Real a(start = 1);\n  Real b(start = 2);\nequation\n  der(a) = -a;\n  der(b) = -b;\nend P;\n"
    run = wait_for(client, start(client, source, 5.0, )["id"])
    # outputStep is 0.01 by default: about 500 samples
    assert run["summary"]["samples"] > 200
    only_b = client.get("/runs/%s/result?variables=b" % run["id"]).json()
    assert list(only_b["variables"]) == ["b"]
    small = client.get("/runs/%s/result?maxPoints=50" % run["id"]).json()
    assert 2 <= len(small["time"]) <= 80
    assert small["time"][0] == 0.0 and small["time"][-1] == pytest.approx(5.0)
    assert client.get("/runs/%s/result?variables=nope" % run["id"]).status_code == 400
    assert client.get("/runs/%s/result?maxPoints=1" % run["id"]).status_code == 400


def test_downsampling_keeps_the_peak(client):
    source = "model S\n  Real x(start = 0);\n  Real v(start = 1);\nequation\n  der(x) = v;\n  der(v) = -x;\nend S;\n"
    run = wait_for(client, start(client, source, 20.0)["id"])
    full = client.get("/runs/%s/result" % run["id"]).json()
    small = client.get("/runs/%s/result?maxPoints=40&variables=x" % run["id"]).json()
    assert max(small["variables"]["x"]) == pytest.approx(max(full["variables"]["x"]), abs=1e-9)


def test_failed_run_reports_diagnostics(client):
    run = wait_for(client, start(client, "model M\n  Real x(start = 1);\nequation\n  der(x) = -y;\nend M;\n")["id"])
    validate("Run", run)
    assert run["state"] == "failed" and run["error"]
    assert run["diagnostics"][0]["line"] == 4
    r = client.get("/runs/%s/result" % run["id"])
    assert r.status_code == 409
    validate("Error", r.json())


def test_terminate_gives_state_terminated(client):
    source = "model T\n  Real x(start = 0);\nequation\n  der(x) = 1;\n  when x > 0.5 then\n    terminate(\"done\");\n  end when;\nend T;\n"
    run = wait_for(client, start(client, source, 3.0)["id"])
    assert run["state"] == "terminated", run
    assert run["termination"]["message"] == "done"
    result = client.get("/runs/%s/result" % run["id"]).json()
    validate("ResultJson", result)
    assert result["terminated"]["message"] == "done"


def test_warnings_for_clamped_options(client):
    run = wait_for(client, start(client, DECAY, 5000.0)["id"], timeout=120)
    validate("Run", run)
    assert any("endTime" in w or "simulation length" in w for w in run["warnings"])


def test_not_found_and_conflicts(client):
    for path in ("/runs/nope", "/runs/nope/result", "/runs/nope/events"):
        r = client.get(path)
        assert r.status_code == 404
        validate("Error", r.json())
    assert client.delete("/runs/nope").status_code == 404
    assert client.get("/nothing-here").status_code == 404


def test_delete_finished_run(client):
    run = wait_for(client, start(client)["id"])
    r = client.delete("/runs/" + run["id"])
    assert r.status_code == 204
    assert client.get("/runs/" + run["id"]).status_code == 404


def test_cancel_running_run():
    with make_client(max_concurrent_runs=1) as c:
        first = start(c)
        r = c.delete("/runs/" + first["id"])
        # either stopped (200) or it had just finished (204)
        assert r.status_code in (200, 204)
        if r.status_code == 200:
            validate("Run", r.json())
            done = wait_for(c, first["id"])
            assert done["state"] == "cancelled"


def test_queue_limit_gives_429():
    with make_client(max_concurrent_runs=1, max_queued_runs=0) as c:
        start(c)
        r = c.post("/runs", json={"model": {"source": DECAY}})
        assert r.status_code == 429
        validate("Error", r.json())


def test_run_timeout():
    with make_client(run_timeout=0.01) as c:
        done = wait_for(c, start(c)["id"])
        assert done["state"] == "failed" and "time limit" in done["error"]


def test_worker_crash_is_survived(monkeypatch):
    monkeypatch.setenv("GENSAS_SERVER_TEST_CRASH", "1")
    with make_client() as c:
        done = wait_for(c, start(c)["id"])
        assert done["state"] == "failed" and "crashed" in done["error"]
        monkeypatch.delenv("GENSAS_SERVER_TEST_CRASH")
        ok = wait_for(c, start(c)["id"])  # the server still works
        assert ok["state"] == "completed"


def test_event_stream(client):
    run = start(client, DECAY, 3.0)
    events = []
    with client.stream("GET", "/runs/%s/events" % run["id"]) as r:
        assert r.status_code == 200 and r.headers["content-type"].startswith("text/event-stream")
        name = None
        for line in r.iter_lines():
            if line.startswith("event: "):
                name = line[7:]
            elif line.startswith("data: "):
                events.append((name, json.loads(line[6:])))
    kinds = [k for k, _ in events]
    assert kinds[0] == "state" and kinds[-1] == "done"
    assert events[-1][1]["state"] == "completed"
    validate("Run", events[-1][1])


def test_token_protects_everything_but_version():
    with make_client(token="s3cret") as c:
        assert c.get("/version").json()["authRequired"] is True
        for method, path in [("get", "/examples"), ("get", "/runs"), ("post", "/check")]:
            r = getattr(c, method)(path) if method == "get" else c.post(path, json={"source": DECAY})
            assert r.status_code == 401, path
            validate("Error", r.json())
        assert c.get("/examples", headers={"Authorization": "Bearer wrong"}).status_code == 401
        assert c.get("/examples", headers={"Authorization": "Bearer s3cret"}).status_code == 200
        assert c.post("/check", json={"source": DECAY}, headers={"Authorization": "Bearer s3cret"}).json()["ok"]


def test_cors_defaults_to_local_origins_only():
    with make_client() as c:  # bound to the loopback address
        ok = c.get("/version", headers={"Origin": "http://localhost:5173"})
        assert ok.headers.get("access-control-allow-origin") == "http://localhost:5173"
        other = c.get("/version", headers={"Origin": "http://evil.example"})
        assert "access-control-allow-origin" not in other.headers
    with make_client(host="0.0.0.0", cors_origins=["https://gensas.example"]) as c:
        assert c.get("/version", headers={"Origin": "https://gensas.example"}).headers["access-control-allow-origin"] == "https://gensas.example"
        assert "access-control-allow-origin" not in c.get("/version", headers={"Origin": "http://localhost:5173"}).headers


def test_finished_runs_are_collected():
    with make_client(max_kept_runs=2) as c:
        ids = [wait_for(c, start(c)["id"])["id"] for _ in range(4)]
        start(c)  # creating a run collects the old ones
        known = {r["id"] for r in c.get("/runs").json()}
        assert ids[0] not in known and ids[-1] in known
