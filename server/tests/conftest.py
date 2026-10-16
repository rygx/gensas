import os
import sys
import time

import pytest
from fastapi.testclient import TestClient

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

from gensas_server.app import create_app  # noqa: E402
from gensas_server.config import Config  # noqa: E402

DECAY = "model D\n  parameter Real k = 1;\n  Real x(start = 1);\nequation\n  der(x) = -k * x;\nend D;\n"
BALL = (
    "model B\n  parameter Real g = 9.81;\n  Real h(start = 1);\n  Real v(start = 0);\nequation\n  der(h) = v;\n  der(v) = -g;\n"
    "  when h < 0 then\n    reinit(v, -0.5 * pre(v));\n  end when;\nend B;\n"
)


def make_client(**overrides):
    cfg = Config(**overrides)
    return TestClient(create_app(cfg))


@pytest.fixture(scope="module")
def client():
    with make_client() as c:
        yield c


def wait_for(client, run_id, states=("completed", "terminated", "failed", "cancelled"), timeout=60.0):
    deadline = time.time() + timeout
    while time.time() < deadline:
        run = client.get("/runs/" + run_id).json()
        if run["state"] in states:
            return run
        time.sleep(0.05)
    raise AssertionError("run %s did not finish: %s" % (run_id, run))


def start(client, source=DECAY, end_time=2.0, **model):
    body = {"model": dict(source=source, **model), "solver": {"endTime": end_time}}
    r = client.post("/runs", json=body)
    assert r.status_code == 202, r.text
    return r.json()
