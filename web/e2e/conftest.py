import os
import socket
import subprocess
import sys
import time
import urllib.request

import pytest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
DIST = os.path.join(ROOT, "web", "dist")


def free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def start_server(*args):
    port = free_port()
    proc = subprocess.Popen(
        [sys.executable, "-m", "gensas_server", "--port", str(port), "--static-dir", DIST, *args],
        cwd=os.path.join(ROOT, "server"),
    )
    url = "http://127.0.0.1:%d" % port
    for _ in range(100):
        try:
            urllib.request.urlopen(url + "/version", timeout=1)
            return proc, url
        except OSError:
            time.sleep(0.1)
    proc.kill()
    raise RuntimeError("the server did not start")


@pytest.fixture(scope="session")
def server():
    assert os.path.isdir(DIST), "run `npm run build` in web/ first"
    proc, url = start_server()
    yield url
    proc.terminate()
    proc.wait(10)


@pytest.fixture(scope="session")
def secured_server():
    proc, url = start_server("--token", "s3cret")
    yield url
    proc.terminate()
    proc.wait(10)


@pytest.fixture(scope="session")
def browser():
    from playwright.sync_api import sync_playwright

    with sync_playwright() as p:
        b = p.chromium.launch()
        yield b
        b.close()


@pytest.fixture
def page(browser):
    context = browser.new_context(viewport={"width": 1400, "height": 900})
    page = context.new_page()
    page.errors = []
    page.on("pageerror", lambda e: page.errors.append(str(e)))
    yield page
    assert page.errors == []
    context.close()
