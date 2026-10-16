"""Starts the server:  python -m gensas_server --port 8765"""

import argparse

import uvicorn

from .app import create_app
from .config import from_env


def main(argv=None):
    cfg = from_env()
    p = argparse.ArgumentParser(prog="gensas_server", description="GenSAS HTTP server (see doc/Server.md)")
    p.add_argument("--host", default=cfg.host, help="address to listen on (default %(default)s; 0.0.0.0 for all interfaces)")
    p.add_argument("--port", type=int, default=cfg.port)
    p.add_argument("--token", default=cfg.token, help="require this bearer token (also GENSAS_TOKEN)")
    p.add_argument("--cors-origin", action="append", default=None, metavar="ORIGIN", help="allow a browser origin (repeatable; also GENSAS_CORS_ORIGINS)")
    p.add_argument("--static-dir", default=cfg.static_dir, help="serve a built web client at /")
    p.add_argument("--max-concurrent-runs", type=int, default=cfg.max_concurrent_runs)
    p.add_argument("--max-queued-runs", type=int, default=cfg.max_queued_runs)
    p.add_argument("--run-timeout", type=float, default=cfg.run_timeout, help="seconds before a run is stopped")
    p.add_argument("--max-source-bytes", type=int, default=cfg.max_source_bytes)
    p.add_argument("--work-dir", default=cfg.work_dir, help="where results are stored (default: a temporary directory)")
    args = p.parse_args(argv)
    for name in ("host", "port", "token", "static_dir", "max_concurrent_runs", "max_queued_runs", "run_timeout", "max_source_bytes", "work_dir"):
        setattr(cfg, name, getattr(args, name))
    if args.cors_origin:
        cfg.cors_origins = args.cors_origin
    if cfg.token is None and cfg.host not in ("127.0.0.1", "localhost", "::1"):
        print("warning: listening on %s without a token; anyone who can reach the port can run simulations" % cfg.host)
    uvicorn.run(create_app(cfg), host=cfg.host, port=cfg.port, log_level="info")


if __name__ == "__main__":
    main()
