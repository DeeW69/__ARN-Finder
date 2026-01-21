"""Utility helpers for interacting with the filesystem."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


def ensure_dir(path: Path) -> Path:
    """Create the directory if it does not exist and return it."""

    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(path: Path, obj: Any) -> None:
    """Write data as JSON with UTF-8 encoding."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2)


def write_text(path: Path, text: str) -> None:
    """Write plain text content to disk."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_lines(path: Path, lines: Iterable[str]) -> None:
    """Write an iterable of strings as newline-delimited text."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for line in lines:
            handle.write(f"{line}\n")


def now_iso() -> str:
    """Return current UTC timestamp as ISO string."""

    return datetime.now(timezone.utc).isoformat()


def sha256_digest(value: str) -> str:
    """Return a deterministic SHA256 hex digest of the provided value."""

    return hashlib.sha256(value.encode("utf-8")).hexdigest()
