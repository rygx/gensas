"""The OpenAPI file as a test oracle: responses are validated against its schemas."""

import os

import jsonschema
import yaml

SPEC_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "api", "openapi.yaml")

with open(SPEC_PATH, encoding="utf-8") as f:
    SPEC = yaml.safe_load(f)

_resolver = jsonschema.RefResolver.from_schema(SPEC)


def validate(schema_name, instance):
    """Raises if `instance` does not match components/schemas/<schema_name>."""
    jsonschema.validate(instance, {"$ref": "#/components/schemas/" + schema_name}, resolver=_resolver)
