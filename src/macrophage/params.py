"""Parameter loader — reads data/parameters.toml.

Every parameter the deterministic core uses comes from here. A missing
parameter raises rather than substituting a guess (anti-fabrication discipline,
borrowed from nature-skills).
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover
    import tomli as tomllib  # type: ignore

_PARAMS: dict[str, Any] | None = None


class MissingParameterError(KeyError):
    """Raised when a required parameter is absent or marked MISSING."""


def _default_path() -> Path:
    here = Path(__file__).resolve()
    # src/macrophage/params.py -> repo_root/data/parameters.toml
    return here.parents[2] / "data" / "parameters.toml"


def load(path: Path | None = None) -> dict[str, Any]:
    """Load (and memoize) parameters from data/parameters.toml."""
    global _PARAMS
    if _PARAMS is not None and path is None:
        return _PARAMS
    p = path or _default_path()
    with open(p, "rb") as fh:
        params = tomllib.load(fh)
    if path is None:
        _PARAMS = params
    return params


def get(*keys: str, params: dict[str, Any] | None = None) -> Any:
    """Look up a nested key; raise MissingParameterError on absent or MISSING."""
    cur = params or load()
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            raise MissingParameterError(f"missing parameter: {'.'.join(keys)}")
        cur = cur[k]
    if cur == "MISSING":
        raise MissingParameterError(
            f"parameter {'.'.join(keys)} is marked MISSING — fill it in "
            "data/parameters.toml with a citation, do not guess."
        )
    return cur
