"""The HTTP application. The contract is api/openapi.yaml; tests/test_contract.py compares them."""

import asyncio
import hmac
import json
from contextlib import asynccontextmanager
from typing import Dict, Literal, Optional, Union

from fastapi import Depends, FastAPI, Query, Request
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse, Response, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, ConfigDict, Field
from starlette.concurrency import run_in_threadpool
from starlette.exceptions import HTTPException

from . import __version__
from .config import Config
from .examples import EXAMPLES
from .results import Result, ResultError, to_f64, to_json
from .runs import ACTIVE, HAS_RESULT, Busy, RunManager

API_VERSION = 1
Numbers = Dict[str, Union[bool, float]]


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


class ModelSpec(_Strict):
    source: str
    className: Optional[str] = None
    piecewise: Optional[Literal["events", "smooth"]] = None
    smoothEps: Optional[float] = Field(default=None, gt=0)
    parameters: Optional[Numbers] = None
    start: Optional[Numbers] = None


class SolverSpec(_Strict):
    endTime: Optional[float] = None
    segment: Optional[float] = None
    level: Optional[int] = None
    alphaTol: Optional[float] = None
    errorTol: Optional[float] = None


class OutputSpec(_Strict):
    step: Optional[float] = None


class RunRequest(_Strict):
    kind: Literal["modelica"] = "modelica"
    model: ModelSpec
    solver: Optional[SolverSpec] = None
    output: Optional[OutputSpec] = None


def _drop_none(d):
    return {k: v for k, v in d.items() if v is not None}


def build_job(model, solver=None, output=None):
    """The engine's job dictionary (doc/JobFormat.md) for a request."""
    job = {
        "gensas": {"format": "job", "version": 1},
        "kind": "modelica",
        "model": _drop_none(model.model_dump()),
    }
    if solver is not None and _drop_none(solver.model_dump()):
        job["solver"] = _drop_none(solver.model_dump())
    if output is not None and _drop_none(output.model_dump()):
        job["output"] = _drop_none(output.model_dump())
    return job


def _error(status, code, message):
    return JSONResponse(status_code=status, content={"error": {"code": code, "message": message}})


class ApiError(Exception):
    def __init__(self, status, code, message):
        self.status, self.code, self.message = status, code, message


def create_app(config: Optional[Config] = None) -> FastAPI:
    config = config or Config()
    manager = RunManager(config)

    @asynccontextmanager
    async def lifespan(app):
        yield
        manager.close()  # stops the worker processes and removes the stored results

    app = FastAPI(title="GenSAS server", version=__version__, docs_url=None, redoc_url=None, openapi_url=None, lifespan=lifespan)
    app.state.manager = manager
    app.state.config = config

    # ---- errors ---------------------------------------------------------------------------

    @app.exception_handler(ApiError)
    async def _api_error(request, exc):
        return _error(exc.status, exc.code, exc.message)

    @app.exception_handler(RequestValidationError)
    async def _invalid(request, exc):
        first = exc.errors()[0] if exc.errors() else {}
        where = ".".join(str(p) for p in first.get("loc", ()) if p != "body")
        return _error(400, "invalid_request", (where + ": " if where else "") + str(first.get("msg", "invalid request")))

    @app.exception_handler(HTTPException)
    async def _http(request, exc):
        codes = {401: "unauthorized", 404: "not_found", 405: "method_not_allowed"}
        return _error(exc.status_code, codes.get(exc.status_code, "error"), str(exc.detail))

    # ---- security and CORS ----------------------------------------------------------------

    def require_token(request: Request):
        if config.token is None:
            return
        header = request.headers.get("authorization", "")
        supplied = header[7:] if header.lower().startswith("bearer ") else ""
        if not hmac.compare_digest(supplied.encode(), config.token.encode()):
            raise ApiError(401, "unauthorized", "a valid bearer token is required")

    auth = [Depends(require_token)]
    if config.cors_origins or config.loopback:
        app.add_middleware(
            CORSMiddleware,
            allow_origins=config.cors_origins,
            allow_origin_regex=None if config.cors_origins else r"^https?://(localhost|127\.0\.0\.1)(:\d+)?$",
            allow_methods=["GET", "POST", "DELETE"],
            allow_headers=["Authorization", "Content-Type"],
            expose_headers=["X-Gensas-Samples", "Location"],
        )

    def check_size(model: ModelSpec):
        if len(model.source.encode("utf-8")) > config.max_source_bytes:
            raise ApiError(413, "too_large", "the model text is larger than %d bytes" % config.max_source_bytes)

    def find(run_id):
        run = manager.get(run_id)
        if run is None:
            raise ApiError(404, "not_found", "no run with this id (it may have been removed)")
        return run

    # ---- endpoints ------------------------------------------------------------------------

    @app.get("/version")
    def version():
        import gensas

        return {
            "server": __version__,
            "api": API_VERSION,
            "engine": {"version": gensas.version(), "apiVersion": gensas.api_version()},
            "authRequired": config.token is not None,
            "limits": {
                "maxSourceBytes": config.max_source_bytes,
                "maxConcurrentRuns": config.max_concurrent_runs,
                "maxQueuedRuns": config.max_queued_runs,
                "runTimeoutSeconds": config.run_timeout,
                "checkTimeoutSeconds": config.check_timeout,
                "retentionSeconds": config.retention,
            },
        }

    @app.get("/examples", dependencies=auth)
    def examples():
        return EXAMPLES

    @app.post("/check", dependencies=auth)
    async def check(spec: ModelSpec):
        check_size(spec)
        try:
            return await run_in_threadpool(manager.check, build_job(spec))
        except Busy as err:
            raise ApiError(503, "busy", str(err))

    @app.post("/runs", status_code=202, dependencies=auth)
    def create_run(req: RunRequest, response: Response):
        check_size(req.model)
        try:
            run = manager.create(build_job(req.model, req.solver, req.output))
        except Busy as err:
            raise ApiError(429, "too_many_runs", str(err))
        response.headers["Location"] = "/runs/" + run.id
        return manager.snapshot(run)

    @app.get("/runs", dependencies=auth)
    def list_runs():
        return manager.list()

    @app.get("/runs/{run_id}", dependencies=auth)
    def get_run(run_id: str):
        return manager.snapshot(find(run_id))

    @app.delete("/runs/{run_id}", dependencies=auth)
    def delete_run(run_id: str):
        run = find(run_id)
        if manager.cancel(run):
            return manager.snapshot(run)
        manager.delete(run)
        return Response(status_code=204)

    @app.get("/runs/{run_id}/events", dependencies=auth)
    async def run_events(run_id: str):
        run = find(run_id)

        async def stream():
            last_state, last_progress, idle = None, None, 0.0
            while True:
                snap = manager.snapshot(run)
                if snap["state"] != last_state:
                    last_state = snap["state"]
                    yield "event: state\ndata: %s\n\n" % json.dumps(snap)
                progress = (snap["progress"]["time"], snap["progress"]["end"])
                if progress != last_progress and snap["state"] == "running":
                    last_progress = progress
                    yield "event: progress\ndata: %s\n\n" % json.dumps(snap["progress"])
                    idle = 0.0
                if snap["state"] not in ACTIVE:
                    yield "event: done\ndata: %s\n\n" % json.dumps(manager.snapshot(run))
                    return
                await asyncio.sleep(0.1)
                idle += 0.1
                if idle >= 15.0:
                    idle = 0.0
                    yield ": keep-alive\n\n"

        return StreamingResponse(stream(), media_type="text/event-stream", headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"})

    @app.get("/runs/{run_id}/result", dependencies=auth)
    def run_result(
        run_id: str,
        format: Literal["json", "csv", "mat", "f64"] = "json",
        variables: Optional[str] = None,
        maxPoints: Optional[int] = Query(default=None, ge=2),
    ):
        run = find(run_id)
        if run.state not in HAS_RESULT:
            raise ApiError(409, "no_result", "the run has no result (state: %s)" % run.state)
        if format == "csv":
            return FileResponse(run.directory + "/result.csv", media_type="text/csv", filename="%s.csv" % run.id)
        if format == "mat":
            return FileResponse(run.directory + "/result.mat", media_type="application/octet-stream", filename="%s.mat" % run.id)
        names = [n for n in (variables or "").split(",") if n] or None
        try:
            result = Result(run.directory)
            if format == "f64":
                body, samples = to_f64(result, names, maxPoints)
                return Response(body, media_type="application/octet-stream", headers={"X-Gensas-Samples": str(samples)})
            return Response(json.dumps(to_json(result, names, maxPoints)), media_type="application/json")
        except ResultError as err:
            raise ApiError(400, "invalid_request", str(err))

    if config.static_dir:
        app.mount("/", StaticFiles(directory=config.static_dir, html=True), name="client")
    return app
