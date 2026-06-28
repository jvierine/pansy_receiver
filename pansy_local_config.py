"""Local runtime configuration helpers for PANSY analysis scripts."""

from __future__ import annotations

import os
from pathlib import Path


DEFAULT_ANALYSIS_CONFIG = Path("~/.config/pansy-receiver/pansy-analysis.env").expanduser()


def analysis_config_path() -> Path:
    return Path(os.environ.get("PANSY_ANALYSIS_CONFIG", DEFAULT_ANALYSIS_CONFIG)).expanduser()


def load_env_config(path: Path | None = None) -> dict[str, str]:
    cfg_path = analysis_config_path() if path is None else Path(path).expanduser()
    if not cfg_path.exists():
        return {}
    values: dict[str, str] = {}
    for raw_line in cfg_path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip().strip("'\"")
        if key:
            values[key] = value
    return values


def config_value(name: str, default: str | None = None, path: Path | None = None) -> str | None:
    if name in os.environ:
        return os.environ[name]
    return load_env_config(path).get(name, default)


def config_path(name: str, default: str | Path | None = None, path: Path | None = None) -> Path | None:
    value = config_value(name, None, path)
    if value is None:
        value = None if default is None else str(default)
    return None if value is None else Path(value).expanduser()
