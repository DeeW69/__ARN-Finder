"""Extraction de motifs (k-mers) V3 — port V1 + progress Rich."""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from rich.progress import Progress, BarColumn, TextColumn, TaskProgressColumn, TimeElapsedColumn

from ..io import (
    ensure_dir, now_iso, read_fasta_records,
    write_csv, write_json, get_logger,
)

logger = get_logger("arnfinderv3.motifs")

COMPLEMENT_DNA = str.maketrans({"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"})
COMPLEMENT_RNA = str.maketrans({"A": "U", "U": "A", "C": "G", "G": "C", "N": "N"})


@dataclass(slots=True)
class MotifRequest:
    input_fasta: Path
    out_dir: Path
    k: int = 9
    top: int = 200
    min_count: int = 1
    ignore_n: bool = True
    canonical: bool = False
    per_record: bool = True
    alphabet: str = "AUTO"


@dataclass(slots=True)
class MotifResult:
    motifs_csv: Path
    motifs_by_record_csv: Optional[Path]
    summary_json: Path
    total_records: int
    total_kmers: int
    unique_kmers: int
    unique_kmers_after_filters: int


def compute_motifs(request: MotifRequest) -> MotifResult:
    out_dir = ensure_dir(request.out_dir)
    motifs_csv = out_dir / "motifs.csv"
    motifs_by_record_csv = out_dir / "motifs_by_record.csv"
    summary_json = out_dir / "motifs_summary.json"

    alphabet_mode = request.alphabet.upper()
    if alphabet_mode not in {"DNA", "RNA", "AUTO"}:
        raise ValueError(f"Alphabet non supporté : {request.alphabet}")

    canonical = request.canonical
    if canonical and alphabet_mode == "RNA":
        logger.warning("Mode canonique désactivé pour l'alphabet RNA.")
        canonical = False

    total_records = 0
    total_kmers = 0
    global_counts: Counter[str] = Counter()
    record_presence: Counter[str] = Counter()

    with Progress(
        TextColumn("[cyan]motifs[/cyan]"),
        BarColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        transient=True,
    ) as progress:
        task = progress.add_task("extraction k-mers...", total=None)

        for header, sequence in read_fasta_records(request.input_fasta):
            total_records += 1
            progress.advance(task)
            seq = sequence.strip().upper()
            _, per_record_counts = _extract_kmers(seq, request.k, request.ignore_n, canonical, alphabet_mode)
            if not per_record_counts:
                continue
            total_kmers += sum(per_record_counts.values())
            global_counts.update(per_record_counts)
            record_presence.update({kmer: 1 for kmer in per_record_counts})

    if total_kmers == 0:
        logger.warning("Aucun k-mer extrait ; vérifiez k ou les séquences d'entrée.")
        _write_empty_outputs(motifs_csv, motifs_by_record_csv, summary_json, request, total_records)
        return MotifResult(motifs_csv, None, summary_json, total_records, 0, 0, 0)

    unique_kmers = len(global_counts)
    filtered = [
        (kmer, count)
        for kmer, count in global_counts.items()
        if count >= request.min_count
    ]
    filtered.sort(key=lambda x: x[1], reverse=True)
    if request.top > 0:
        filtered = filtered[: request.top]
    whitelist = {kmer for kmer, _ in filtered}

    _write_global_csv(motifs_csv, filtered, record_presence, total_records, total_kmers)

    motifs_by_record_path: Optional[Path] = None
    if request.per_record and whitelist:
        motifs_by_record_path = motifs_by_record_csv
        _write_per_record_csv(
            request.input_fasta, motifs_by_record_path,
            request.k, request.ignore_n, canonical, alphabet_mode, whitelist,
        )

    summary = {
        "timestamp": now_iso(),
        "input_fasta": str(request.input_fasta),
        "params": {
            "k": request.k,
            "top": request.top,
            "min_count": request.min_count,
            "ignore_n": request.ignore_n,
            "canonical": canonical,
            "per_record": request.per_record,
            "alphabet": alphabet_mode,
        },
        "total_records": total_records,
        "total_kmers_considered": total_kmers,
        "unique_kmers": unique_kmers,
        "unique_kmers_after_filters": len(whitelist),
    }
    write_json(summary_json, summary)

    logger.info(
        "Motifs : %s records, %s k-mers totaux, %s uniques (après filtres : %s)",
        total_records, total_kmers, unique_kmers, len(whitelist),
    )

    return MotifResult(
        motifs_csv=motifs_csv,
        motifs_by_record_csv=motifs_by_record_path,
        summary_json=summary_json,
        total_records=total_records,
        total_kmers=total_kmers,
        unique_kmers=unique_kmers,
        unique_kmers_after_filters=len(whitelist),
    )


# ── Helpers ───────────────────────────────────────────────────────────────────

def _extract_kmers(
    sequence: str, k: int, ignore_n: bool, canonical: bool, alphabet_mode: str,
) -> tuple[list[str], Counter[str]]:
    counts: Counter[str] = Counter()
    if len(sequence) < k or k <= 0:
        return [], counts
    for idx in range(len(sequence) - k + 1):
        kmer = sequence[idx: idx + k]
        if ignore_n and "N" in kmer:
            continue
        if canonical:
            kmer = _canonical_form(kmer, alphabet_mode)
        counts[kmer] += 1
    return list(counts.keys()), counts


def _canonical_form(kmer: str, alphabet_mode: str) -> str:
    table = COMPLEMENT_DNA if alphabet_mode in {"AUTO", "DNA"} else COMPLEMENT_RNA
    revcomp = kmer.translate(table)[::-1]
    return min(kmer, revcomp)


def _write_global_csv(
    path: Path,
    filtered: list[tuple[str, int]],
    record_presence: Counter[str],
    total_records: int,
    total_kmers: int,
) -> None:
    rows = []
    for kmer, count in filtered:
        freq = count / total_kmers if total_kmers else 0.0
        n_records = record_presence.get(kmer, 0)
        idf = math.log((1 + total_records) / (1 + n_records)) if total_records else 0.0
        rows.append([kmer, str(count), f"{freq:.6f}", str(n_records), f"{idf:.6f}"])
    write_csv(path, ["kmer", "count", "freq", "n_records", "idf_like"], rows)


def _write_per_record_csv(
    fasta_path: Path, path: Path,
    k: int, ignore_n: bool, canonical: bool, alphabet_mode: str, whitelist: set[str],
) -> None:
    rows = []
    for header, sequence in read_fasta_records(fasta_path):
        record_id = header.split(" ", 1)[0]
        seq = sequence.strip().upper()
        _, counts = _extract_kmers(seq, k, ignore_n, canonical, alphabet_mode)
        total = sum(counts.values())
        if total == 0:
            continue
        for kmer, count in counts.items():
            if kmer not in whitelist:
                continue
            rows.append([record_id, kmer, str(count), f"{count / total:.6f}"])
    write_csv(path, ["record_id", "kmer", "count", "freq_in_record"], rows)


def _write_empty_outputs(
    motifs_csv: Path, per_record_csv: Path, summary_json: Path,
    request: MotifRequest, total_records: int
) -> None:
    write_csv(motifs_csv, ["kmer", "count", "freq", "n_records", "idf_like"], [])
    if request.per_record:
        write_csv(per_record_csv, ["record_id", "kmer", "count", "freq_in_record"], [])
    write_json(summary_json, {
        "timestamp": now_iso(),
        "input_fasta": str(request.input_fasta),
        "params": {
            "k": request.k, "top": request.top, "min_count": request.min_count,
            "ignore_n": request.ignore_n, "canonical": request.canonical,
            "per_record": request.per_record, "alphabet": request.alphabet,
        },
        "total_records": total_records,
        "total_kmers_considered": 0,
        "unique_kmers": 0,
        "unique_kmers_after_filters": 0,
        "notes": ["Aucun k-mer extrait"],
    })
