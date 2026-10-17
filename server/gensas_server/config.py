"""Server settings. Every setting has an environment variable (GENSAS_*) and a command-line flag."""

import os
from dataclasses import dataclass, field, fields
from typing import List, Optional


@dataclass
class Config:
    host: str = "127.0.0.1"
    port: int = 8764
    token: Optional[str] = None  # if set, requests need "Authorization: Bearer <token>"
    cors_origins: List[str] = field(default_factory=list)  # allowed browser origins; "*" allows all
    static_dir: Optional[str] = None  # a built web client to serve at /
    max_concurrent_runs: int = 2  # worker processes running at the same time
    max_queued_runs: int = 16  # runs waiting for a worker; more are refused with 429
    max_concurrent_checks: int = 2
    run_timeout: float = 120.0  # wall-clock seconds before a run is stopped
    check_timeout: float = 10.0
    max_source_bytes: int = 200_000
    max_kept_runs: int = 50  # finished runs kept (oldest are removed first)
    retention: float = 3600.0  # seconds a finished run is kept
    work_dir: Optional[str] = None  # where run results are stored (default: a temporary directory)
    rate_limit: int = 0  # POST /check and POST /runs allowed per minute for each client address; 0 = no limit
    trust_proxy: bool = False  # take the client address from X-Forwarded-For (only behind a proxy you control)

    @property
    def loopback(self) -> bool:
        return self.host in ("127.0.0.1", "localhost", "::1")


_ENV = {
    "host": ("GENSAS_HOST", str),
    "port": ("GENSAS_PORT", int),
    "token": ("GENSAS_TOKEN", str),
    "static_dir": ("GENSAS_STATIC_DIR", str),
    "max_concurrent_runs": ("GENSAS_MAX_CONCURRENT_RUNS", int),
    "max_queued_runs": ("GENSAS_MAX_QUEUED_RUNS", int),
    "max_concurrent_checks": ("GENSAS_MAX_CONCURRENT_CHECKS", int),
    "run_timeout": ("GENSAS_RUN_TIMEOUT", float),
    "check_timeout": ("GENSAS_CHECK_TIMEOUT", float),
    "max_source_bytes": ("GENSAS_MAX_SOURCE_BYTES", int),
    "max_kept_runs": ("GENSAS_MAX_KEPT_RUNS", int),
    "retention": ("GENSAS_RETENTION", float),
    "work_dir": ("GENSAS_WORK_DIR", str),
    "rate_limit": ("GENSAS_RATE_LIMIT", int),
}


def from_env(environ=None) -> Config:
    env = os.environ if environ is None else environ
    cfg = Config()
    for f in fields(Config):
        if f.name in _ENV:
            name, kind = _ENV[f.name]
            if env.get(name):
                setattr(cfg, f.name, kind(env[name]))
    if not env.get("GENSAS_PORT") and env.get("PORT"):  # Cloud Run and similar platforms name the port this way
        cfg.port = int(env["PORT"])
    if env.get("GENSAS_TRUST_PROXY", "").lower() in ("1", "true", "yes"):
        cfg.trust_proxy = True
    if env.get("GENSAS_CORS_ORIGINS"):
        cfg.cors_origins = [o.strip() for o in env["GENSAS_CORS_ORIGINS"].split(",") if o.strip()]
    return cfg
