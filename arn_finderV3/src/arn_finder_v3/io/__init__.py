"""Utilitaires I/O V3 — logging Rich, helpers FASTA/CSV/JSON."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Generator, Iterable, Tuple

from rich.logging import RichHandler


# ── Logging ──────────────────────────────────────────────────────────────────

def configure_logging(verbose: bool = False) -> None:
    """Configure logging global avec RichHandler (sortie colorée, tracebacks)."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(rich_tracebacks=True, show_path=False)],
        force=True,
    )


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(name)


# ── Temps ────────────────────────────────────────────────────────────────────

def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


# ── Filesystem ───────────────────────────────────────────────────────────────

def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def sha256_digest(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def sha256_of_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


# ── JSON ─────────────────────────────────────────────────────────────────────

def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2, ensure_ascii=False)


def write_lines(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for line in lines:
            handle.write(f"{line}\n")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


# ── CSV ──────────────────────────────────────────────────────────────────────

def write_csv(path: Path, headers: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)


# ── FASTA ────────────────────────────────────────────────────────────────────

def read_fasta_records(path: Path) -> Generator[Tuple[str, str], None, None]:
    """Yield (header, sequence) depuis un fichier FASTA."""
    header: str | None = None
    seq_chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header is not None:
        yield header, "".join(seq_chunks)


def write_fasta_record(handle: Any, header: str, sequence: str, width: int = 80) -> None:
    """Écrit un enregistrement FASTA dans un handle ouvert."""
    handle.write(f">{header}\n")
    for idx in range(0, len(sequence), width):
        handle.write(sequence[idx: idx + width] + "\n")


def write_fasta(records: list[Tuple[str, str]], path: Path, width: int = 80) -> None:
    """Écrit une liste (header, sequence) dans un fichier FASTA."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for header, seq in records:
            write_fasta_record(f, header, seq, width)
