"""Filtrage qualité V3 — port V1 + filtre GC + progress Rich."""

from __future__ import annotations

import json
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Tuple

from rich.progress import Progress, BarColumn, TextColumn, TaskProgressColumn, TimeElapsedColumn

from ..io import (
    ensure_dir, now_iso, sha256_digest,
    read_fasta_records, write_fasta_record, write_csv, write_json,
    get_logger,
)

logger = get_logger("arnfinderv3.filter")

AMBIGUOUS_CHARS = set("RYKMSWBDHV")
DNA_ALPHABET = set("ACGTN")
RNA_ALPHABET = set("ACGUN")
ACCESSION_PATTERN = re.compile(r"([A-Z]{1,}[0-9]+(?:\.[0-9]+)?)")


@dataclass(slots=True)
class FilterRequest:
    input_fasta: Path
    out_dir: Path
    min_len: int = 200
    max_len: Optional[int] = None
    max_n_fraction: float = 0.05
    min_gc: float = 0.0       # V3 new
    max_gc: float = 1.0       # V3 new
    alphabet: str = "AUTO"
    dedupe: bool = False
    metadata_path: Optional[Path] = None
    write_report: bool = True
    max_per_organism: int = 0


@dataclass(slots=True)
class FilterResult:
    total_records: int
    kept_records: int
    dropped_by_reason: Dict[str, int]
    fasta_path: Path
    report_path: Optional[Path]
    metadata_path: Optional[Path]
    manifest_path: Path


def filter_sequences(request: FilterRequest) -> FilterResult:
    out_dir = ensure_dir(request.out_dir)
    fasta_out = out_dir / "filtered_sequences.fasta"
    report_out = out_dir / "filter_report.csv"
    metadata_out = out_dir / "filtered_metadata.jsonl"
    manifest_out = out_dir / "filter_manifest.json"

    total = kept = 0
    reason_counter: Counter[str] = Counter()
    kept_ids: set[str] = set()
    kept_accessions: set[str] = set()
    seen_sequences: set[str] = set()
    duplicate_map: dict[str, str] = {}
    organism_limit = max(request.max_per_organism, 0)
    organism_kept: Counter[str] = Counter()
    raw_organism_counts: Counter[str] = Counter()
    metadata_lookup = _load_metadata_organisms(request.metadata_path)
    unknown_warning_emitted = False

    alphabet_mode = request.alphabet.upper()
    if alphabet_mode not in {"DNA", "RNA", "AUTO"}:
        raise ValueError(f"Alphabet non supporté : {request.alphabet}")

    report_rows: list[list[str]] = []

    with fasta_out.open("w", encoding="utf-8") as fasta_handle, \
         Progress(
             TextColumn("[cyan]filter[/cyan]"),
             BarColumn(),
             TaskProgressColumn(),
             TimeElapsedColumn(),
             transient=True,
         ) as progress:

        task = progress.add_task("filtrage...", total=None)

        for header, raw_seq in read_fasta_records(request.input_fasta):
            total += 1
            progress.advance(task)
            record_id, accession = _extract_ids(header)
            organism, unknown_warning_emitted = _resolve_organism(
                accession, record_id, header, metadata_lookup, unknown_warning_emitted
            )
            raw_organism_counts[organism] += 1
            cleaned = _clean_sequence(raw_seq)
            detected_alphabet = _determine_alphabet(cleaned, alphabet_mode, record_id)
            normalized = _normalize_sequence(cleaned, detected_alphabet, record_id)
            quality = _compute_metrics(normalized)

            kept_flag, reason = _apply_filters(
                quality, normalized, request, organism, organism_limit,
                organism_kept, seen_sequences, duplicate_map, accession, record_id
            )

            if kept_flag:
                kept += 1
                kept_ids.add(record_id)
                if accession:
                    kept_accessions.add(accession)
                header_line = (
                    f"{record_id} | accession={accession or 'NA'} | "
                    f"len={quality['length']} | n_frac={quality['n_frac']:.3f} | "
                    f"gc={quality['gc_percent']:.1f}"
                )
                write_fasta_record(fasta_handle, header_line, normalized)

            reason_counter[reason] += 1
            if request.write_report:
                report_rows.append([
                    record_id,
                    accession or "",
                    str(quality["length"]),
                    str(quality["n_count"]),
                    f"{quality['n_frac']:.5f}",
                    f"{quality['gc_percent']:.2f}",
                    str(kept_flag).lower(),
                    reason,
                    header,
                ])

    report_path: Optional[Path] = None
    if request.write_report:
        write_csv(
            report_out,
            ["record_id", "accession", "length", "n_count", "n_frac",
             "gc_percent", "kept", "reason", "original_header"],
            report_rows,
        )
        report_path = report_out

    metadata_written: Optional[Path] = None
    if request.metadata_path and request.metadata_path.exists():
        metadata_written = metadata_out
        _filter_metadata(request.metadata_path, metadata_written, kept_accessions, kept_ids)

    manifest = {
        "timestamp": now_iso(),
        "inputs": {
            "fasta": str(request.input_fasta),
            "metadata": str(request.metadata_path) if request.metadata_path else None,
        },
        "outputs": {
            "sequences": str(fasta_out),
            "report": str(report_path) if report_path else None,
            "metadata": str(metadata_written) if metadata_written else None,
        },
        "params": {
            "min_len": request.min_len,
            "max_len": request.max_len,
            "max_n_frac": request.max_n_fraction,
            "min_gc": request.min_gc,
            "max_gc": request.max_gc,
            "alphabet": alphabet_mode,
            "dedupe": request.dedupe,
            "max_per_organism": organism_limit,
        },
        "total_records": total,
        "kept_records": kept,
        "dropped_records": total - kept,
        "reason_counts": dict(reason_counter),
    }
    write_json(manifest_out, manifest)

    logger.info("Filter : total=%s kept=%s dropped=%s", total, kept, total - kept)
    if kept == 0:
        logger.warning("Aucune séquence n'a passé le filtre.")

    return FilterResult(
        total_records=total,
        kept_records=kept,
        dropped_by_reason=dict(reason_counter),
        fasta_path=fasta_out,
        report_path=report_path,
        metadata_path=metadata_written,
        manifest_path=manifest_out,
    )


# ── Helpers ───────────────────────────────────────────────────────────────────

def _apply_filters(
    quality: dict,
    normalized: str,
    request: FilterRequest,
    organism: str,
    organism_limit: int,
    organism_kept: Counter[str],
    seen_sequences: set[str],
    duplicate_map: dict[str, str],
    accession: Optional[str],
    record_id: str,
) -> Tuple[bool, str]:
    if quality["length"] == 0:
        return False, "empty_after_clean"
    if quality["length"] < request.min_len:
        return False, "too_short"
    if request.max_len and quality["length"] > request.max_len:
        return False, "too_long"
    if quality["n_frac"] > request.max_n_fraction:
        return False, "too_many_N"
    gc = quality["gc_percent"] / 100.0
    if gc < request.min_gc:
        return False, "gc_too_low"
    if gc > request.max_gc:
        return False, "gc_too_high"
    if request.dedupe:
        seq_hash = sha256_digest(normalized)
        if seq_hash in seen_sequences:
            return False, "duplicate"
        seen_sequences.add(seq_hash)
        duplicate_map[seq_hash] = accession or record_id
    if organism_limit > 0 and organism_kept[organism] >= organism_limit:
        return False, "max_per_organism"
    organism_kept[organism] += 1
    return True, "ok"


def _extract_ids(header: str) -> Tuple[str, Optional[str]]:
    match = ACCESSION_PATTERN.search(header)
    accession = match.group(1) if match else None
    record_id = accession or header.split()[0]
    return record_id, accession


def _clean_sequence(seq: str) -> str:
    seq = "".join(ch for ch in seq.upper() if ch.isalpha())
    return "".join("N" if ch in AMBIGUOUS_CHARS else ch for ch in seq)


def _determine_alphabet(seq: str, mode: str, record_id: str) -> str:
    if mode in ("DNA", "RNA"):
        return mode
    has_u, has_t = "U" in seq, "T" in seq
    if has_u and not has_t:
        return "RNA"
    if has_t and has_u:
        logger.warning("Séquence %s : T et U présents, défaut DNA.", record_id)
    return "DNA"


def _normalize_sequence(seq: str, alphabet: str, record_id: str) -> str:
    if alphabet == "DNA":
        seq = seq.replace("U", "T")
        allowed = DNA_ALPHABET
    else:
        seq = seq.replace("T", "U")
        allowed = RNA_ALPHABET
    return "".join(ch if ch in allowed else "N" for ch in seq)


def _compute_metrics(seq: str) -> dict:
    length = len(seq)
    n_count = seq.count("N")
    gc_count = seq.count("G") + seq.count("C")
    denom = max(length - n_count, 1) if length else 1
    gc_percent = gc_count / denom * 100 if length else 0.0
    n_frac = n_count / length if length else 1.0
    return {"length": length, "n_count": n_count, "n_frac": n_frac, "gc_percent": gc_percent}


def _load_metadata_organisms(metadata_path: Optional[Path]) -> dict[str, str]:
    lookup: dict[str, str] = {}
    if not metadata_path or not metadata_path.exists():
        return lookup
    with metadata_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            payload = json.loads(line)
            organism = payload.get("organism") or payload.get("taxname")
            if not organism:
                continue
            for key in (
                payload.get("accession"),
                payload.get("accessionversion"),
                str(payload["uid"]) if payload.get("uid") is not None else None,
            ):
                if key:
                    lookup[str(key)] = organism
    return lookup


def _resolve_organism(
    accession: Optional[str],
    record_id: str,
    header: str,
    lookup: dict[str, str],
    unknown_warning_emitted: bool,
) -> Tuple[str, bool]:
    organism = (
        lookup.get(accession or "")
        or lookup.get(record_id)
        or _infer_organism_from_header(header)
    )
    if organism == "UNKNOWN" and not unknown_warning_emitted:
        logger.warning("Certaines séquences ont un organisme inconnu.")
        unknown_warning_emitted = True
    return organism, unknown_warning_emitted


def _infer_organism_from_header(header: str) -> str:
    if " " not in header:
        return "UNKNOWN"
    remainder = header.split(" ", 1)[1].strip()
    lower = remainder.lower()
    if " voucher" in lower:
        remainder = remainder[: lower.index(" voucher")].strip()
    return remainder or "UNKNOWN"


def _filter_metadata(
    metadata_path: Path, output_path: Path,
    kept_accessions: set[str], kept_ids: set[str]
) -> None:
    with metadata_path.open("r", encoding="utf-8") as source, \
         output_path.open("w", encoding="utf-8") as target:
        for line in source:
            line = line.strip()
            if not line:
                continue
            payload = json.loads(line)
            accession = payload.get("accession") or payload.get("accessionversion")
            uid = str(payload.get("uid")) if payload.get("uid") is not None else None
            if (accession and accession in kept_accessions) or (uid and uid in kept_ids):
                target.write(json.dumps(payload) + "\n")
