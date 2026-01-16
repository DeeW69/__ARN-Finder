"""Minimal FASTA parsing and writing utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Generator, Iterable, Tuple


def read_fasta_records(path: Path) -> Generator[Tuple[str, str], None, None]:
    """Yield (header, sequence) tuples from a FASTA file."""

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


def write_fasta_record(handle, header: str, sequence: str, width: int = 80) -> None:
    """Write a FASTA record to an open file handle."""

    handle.write(f">{header}\n")
    for idx in range(0, len(sequence), width):
        handle.write(sequence[idx : idx + width] + "\n")
