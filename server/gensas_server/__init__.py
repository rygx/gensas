"""GenSAS server: an HTTP API (api/openapi.yaml) over the GenSAS Python interface."""

import os
import sys

__version__ = "0.1.0"

try:
    import gensas  # noqa: F401
except ImportError:  # running from a source checkout: python/ next to server/
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "python"))
