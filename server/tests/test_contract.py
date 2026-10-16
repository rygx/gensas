from fastapi.routing import APIRoute
from openapi_spec_validator import validate as validate_spec

from contract import SPEC
from conftest import make_client


def test_spec_is_valid_openapi():
    validate_spec(SPEC)


def test_routes_match_the_spec():
    with make_client() as client:
        implemented = set()
        for route in client.app.routes:
            if isinstance(route, APIRoute):
                for method in route.methods:
                    implemented.add((method.lower(), route.path.replace("{run_id}", "{id}")))
    documented = {(method, path) for path, item in SPEC["paths"].items() for method in item if method in ("get", "post", "delete")}
    assert implemented == documented
