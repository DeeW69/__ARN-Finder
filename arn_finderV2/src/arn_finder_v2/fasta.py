"""Lightweight FASTA parsing utilities."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, List, Optional


@dataclass
class FastaRecord:
    identifier: str
    description: str
    sequence: str


def read_fasta(path: Path, uppercase: bool = True, convert_u_to_t: bool = False) -> List[FastaRecord]:
    """Read a FASTA file into a list of FastaRecords."""
    return list(parse_fasta_lines(path.read_text(encoding="utf-8").splitlines(), uppercase, convert_u_to_t))


def parse_fasta_lines(
    lines: Iterable[str],
    uppercase: bool = True,
    convert_u_to_t: bool = False,
) -> Iterator[FastaRecord]:
    header: Optional[str] = None
    seq_parts: List[str] = []

    def flush_record() -> Optional[FastaRecord]:
        nonlocal header, seq_parts
        if header is None:
            return None
        identifier, description = _split_header(header)
        sequence = "".join(seq_parts)
        if uppercase:
            sequence = sequence.upper()
        if convert_u_to_t:
            sequence = sequence.replace("U", "T")
            sequence = sequence.replace("u", "T")
        record = FastaRecord(identifier=identifier, description=description, sequence=sequence)
        header = None
        seq_parts = []
        return record

    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            record = flush_record()
            header = line[1:].strip()
            if record:
                yield record
        else:
            seq_parts.append(line)

    last_record = flush_record()
    if last_record:
        yield last_record


def _split_header(header: str) -> tuple[str, str]:
    if not header:
        return "", ""
    parts = header.split(maxsplit=1)
    if len(parts) == 1:
        return parts[0], ""
    return parts[0], parts[1]
