"""Drives the built client in a real browser against a real server."""

import json
import urllib.request

import pytest

TIMEOUT = 30000


def open_app(page, url):
    page.goto(url + "/")
    page.wait_for_selector(".conn.ok", timeout=TIMEOUT)
    page.wait_for_selector(".status-line .ok-text", timeout=TIMEOUT)


def api(url, path):
    return json.load(urllib.request.urlopen(url + path))


def replace_source(page, text):
    page.click(".cm-content")
    page.keyboard.press("Control+A")
    page.keyboard.insert_text(text)


def test_example_runs_and_plots(page, server):
    open_app(page, server)
    assert "Decay" in page.inner_text(".status-line")
    page.get_by_role("button", name="Run").click()
    page.wait_for_selector(".results .plot canvas", timeout=TIMEOUT)
    assert "completed" in page.inner_text(".run-status")
    assert page.locator("ul.legend li").count() == 1


def test_changed_parameter_is_used(page, server):
    open_app(page, server)
    page.get_by_label("k", exact=True).fill("2")
    page.get_by_label("end time").fill("5")
    page.get_by_role("button", name="Run").click()
    page.wait_for_selector(".results .plot canvas", timeout=TIMEOUT)
    runs = api(server, "/runs")
    result = api(server, "/runs/%s/result" % runs[0]["id"])
    assert result["variables"]["x"][-1] == pytest.approx(pytest.importorskip("math").exp(-10.0), abs=1e-4)
    page.get_by_text("reset", exact=True).first.click()
    assert page.get_by_label("k", exact=True).input_value() == "1"


def test_errors_show_with_line_and_block_running(page, server):
    open_app(page, server)
    replace_source(page, "model M\n  Real x(start = 1);\nequation\n  der(x) = -y;\nend M;\n")
    page.wait_for_selector("ul.diagnostics li.error", timeout=TIMEOUT)
    assert "line 4" in page.inner_text("ul.diagnostics")
    assert "undeclared identifier 'y'" in page.inner_text("ul.diagnostics")
    assert page.get_by_role("button", name="Run").is_disabled()
    assert page.locator(".cm-lintRange-error").count() >= 1


def test_events_and_stacked_plots(page, server):
    open_app(page, server)
    page.select_option("select[aria-label=example]", "bouncing_ball")
    page.wait_for_selector(".status-line .ok-text:has-text('with events')", timeout=TIMEOUT)
    page.get_by_role("button", name="Run").click()
    page.wait_for_selector(".results .plot canvas", timeout=TIMEOUT)
    assert page.locator("details summary", has_text="events").count() == 1
    page.get_by_label("one plot per variable").check()
    page.wait_for_selector(".plot.stacked .plot-cell canvas", timeout=TIMEOUT)
    assert page.locator(".plot-cell").count() == 2


def test_downloads(page, server):
    open_app(page, server)
    page.get_by_role("button", name="Run").click()
    page.wait_for_selector(".results .plot canvas", timeout=TIMEOUT)
    with page.expect_download() as d:
        page.get_by_role("button", name="CSV").click()
    assert d.value.suggested_filename.endswith(".csv")


def test_token_flow(page, secured_server):
    page.goto(secured_server + "/")
    page.wait_for_selector(".banner.bad", timeout=TIMEOUT)
    assert "token" in page.inner_text(".banner.bad")
    page.get_by_role("button", name="Server").click()
    page.get_by_label("Access token", exact=False).fill("s3cret")
    page.get_by_role("button", name="Connect").click()
    page.wait_for_selector(".conn.ok", timeout=TIMEOUT)
    page.wait_for_selector(".status-line .ok-text", timeout=TIMEOUT)


def test_unreachable_server_is_explained(page, server):
    page.goto(server + "/")
    page.wait_for_selector(".conn.ok", timeout=TIMEOUT)
    page.get_by_role("button", name="Server").click()
    page.get_by_label("Server address").fill("http://127.0.0.1:1")
    page.get_by_role("button", name="Connect").click()
    page.wait_for_selector(".banner.bad", timeout=TIMEOUT)
    assert "cannot reach the server" in page.inner_text(".banner.bad")
