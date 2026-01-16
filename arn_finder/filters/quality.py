"""Sequence filtering utilities."""

from __future__ import annotations

import csv
import json
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple

from ..io.fasta import read_fasta_records, write_fasta_record
from ..io.paths import ensure_dir, now_iso, sha256_digest, write_json

AMBIGUOUS_CHARS = set("RYKMSWBDHV")
DNA_ALPHABET = set("ACGTN")
RNA_ALPHABET = set("ACGUN")
ACCESSION_PATTERN = re.compile(r"([A-Z]{1,}[0-9]+(?:\.[0-9]+)?)")


@dataclass(slots=True)
class FilterRequest:
    input_fasta: Path
    out_dir: Path
    min_len: int
    max_n_fraction: float
    alphabet: str = "AUTO"
    dedupe: bool = False
    metadata_path: Path | None = None
    write_report: bool = True
    max_per_organism: int = 0


@dataclass(slots=True)
class FilterResult:
    total_records: int
    kept_records: int
    dropped_by_reason: Dict[str, int]
    fasta_path: Path
    report_path: Path | None
    metadata_path: Path | None
    manifest_path: Path


def filter_sequences(request: FilterRequest, logger) -> FilterResult:
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
        raise ValueError(f"Unsupported alphabet: {request.alphabet}")

    report_writer = None
    report_file = None
    if request.write_report:
        report_file = report_out.open("w", newline="", encoding="utf-8")
        report_writer = csv.writer(report_file)
        report_writer.writerow(
            [
                "record_id",
                "accession",
                "length",
                "n_count",
                "n_frac",
                "gc_percent",
                "kept",
                "reason",
                "duplicate_of",
                "original_header",
            ]
        )

    if organism_limit > 0:
        logger.info("Applying max-per-organism limit: %s", organism_limit)

    with fasta_out.open("w", encoding="utf-8") as fasta_handle:
        for header, raw_seq in read_fasta_records(request.input_fasta):
            total += 1
            record_id, accession = _extract_ids(header, logger)
            organism, unknown_warning_emitted = _resolve_organism(
                accession, record_id, header, metadata_lookup, unknown_warning_emitted, logger
            )
            raw_organism_counts[organism] += 1
            cleaned = _clean_sequence(raw_seq)
            detected_alphabet = _determine_alphabet(cleaned, alphabet_mode, record_id, logger)
            normalized = _normalize_sequence(cleaned, detected_alphabet, record_id, logger)
            quality = _compute_metrics(normalized)

            if quality["length"] == 0:
                kept_flag = False
                reason = "empty_after_clean"
                duplicate_of = ""
            elif quality["length"] < request.min_len:
                kept_flag = False
                reason = "too_short"
                duplicate_of = ""
            elif quality["n_frac"] > request.max_n_fraction:
                kept_flag = False
                reason = "too_many_N"
                duplicate_of = ""
            else:
                seq_hash = sha256_digest(normalized) if request.dedupe else None
                if request.dedupe and seq_hash in seen_sequences:
                    kept_flag = False
                    reason = "duplicate"
                    duplicate_of = duplicate_map.get(seq_hash, "")
                else:
                    if organism_limit > 0 and organism_kept[organism] >= organism_limit:
                        kept_flag = False
                        reason = "max_per_organism"
                        duplicate_of = ""
                    else:
                        kept_flag = True
                        reason = "ok"
                        duplicate_of = ""
                        if request.dedupe and seq_hash is not None:
                            seen_sequences.add(seq_hash)
                            duplicate_map[seq_hash] = accession or record_id
                        organism_kept[organism] += 1

            if kept_flag:
                kept += 1
                kept_ids.add(record_id)
                if accession:
                    kept_accessions.add(accession)
                header_line = (
                    f"{record_id} | accession={accession or 'NA'} | "
                    f"len={quality['length']} | n_frac={quality['n_frac']:.3f} | gc={quality['gc_percent']:.1f}"
                )
                write_fasta_record(fasta_handle, header_line, normalized)
            reason_counter[reason] += 1

            if report_writer:
                metrics = quality
                report_writer.writerow(
                    [
                        record_id,
                        accession or "",
                        metrics["length"],
                        metrics["n_count"],
                        f"{metrics['n_frac']:.5f}",
                        f"{metrics['gc_percent']:.2f}",
                        str(kept_flag).lower(),
                        reason,
                        duplicate_of or "",
                        header,
                    ]
                )

    if report_file:
        report_file.close()
    else:
        report_out = None

    metadata_written = None
    if request.metadata_path and request.metadata_path.exists():
        metadata_written = metadata_out
        _filter_metadata(request.metadata_path, metadata_written, kept_accessions, kept_ids)
    else:
        metadata_written = None

    manifest = {
        "timestamp": now_iso(),
        "inputs": {
            "fasta": str(request.input_fasta),
            "metadata": str(request.metadata_path) if request.metadata_path else None,
        },
        "outputs": {
            "sequences": str(fasta_out),
            "report": str(report_out) if report_out else None,
            "metadata": str(metadata_written) if metadata_written else None,
        },
        "params": {
            "min_len": request.min_len,
            "max_n_frac": request.max_n_fraction,
            "alphabet": alphabet_mode,
            "dedupe": request.dedupe,
            "write_report": request.write_report,
            "max_per_organism": organism_limit,
        },
        "total_records": total,
        "kept_records": kept,
        "dropped_records": total - kept,
        "reason_counts": dict(reason_counter),
        "duplicates": reason_counter.get("duplicate", 0),
        "dropped_max_per_organism": reason_counter.get("max_per_organism", 0),
    }
    write_json(manifest_out, manifest)

    _log_summary(
        logger,
        total,
        kept,
        reason_counter,
        fasta_out,
        metadata_written,
        report_out,
        raw_organism_counts,
    )
    if kept == 0:
        logger.warning("No sequences passed filtering; outputs are empty.")

    return FilterResult(
        total_records=total,
        kept_records=kept,
        dropped_by_reason=dict(reason_counter),
        fasta_path=fasta_out,
        report_path=report_out,
        metadata_path=metadata_written,
        manifest_path=manifest_out,
    )


def _extract_ids(header: str, logger) -> Tuple[str, str | None]:
    accession = None
    match = ACCESSION_PATTERN.search(header)
    if match:
        accession = match.group(1)
    else:
        logger.warning("Unable to detect accession in header: %s", header)
    record_id = accession or header.split()[0]
    return record_id, accession


def _clean_sequence(seq: str) -> str:
    seq = "".join(ch for ch in seq.upper() if ch.isalpha())
    cleaned = []
    for ch in seq:
        if ch in AMBIGUOUS_CHARS:
            cleaned.append("N")
        else:
            cleaned.append(ch)
    return "".join(cleaned)


def _determine_alphabet(seq: str, mode: str, record_id: str, logger) -> str:
    if mode in ("DNA", "RNA"):
        return mode
    has_u = "U" in seq
    has_t = "T" in seq
    if has_u and not has_t:
        return "RNA"
    if has_t and not has_u:
        return "DNA"
    if has_t and has_u:
        logger.warning("Sequence %s contains both T and U; defaulting to DNA", record_id)
        return "DNA"
    return "DNA"


def _normalize_sequence(seq: str, alphabet: str, record_id: str, logger) -> str:
    if alphabet == "DNA":
        if "U" in seq:
            logger.warning("Converting U to T for record %s (DNA alphabet)", record_id)
        seq = seq.replace("U", "T")
        allowed = DNA_ALPHABET
    else:
        if "T" in seq:
            logger.warning("Converting T to U for record %s (RNA alphabet)", record_id)
        seq = seq.replace("T", "U")
        allowed = RNA_ALPHABET

    final_chars = []
    for ch in seq:
        final_chars.append(ch if ch in allowed else "N")
    return "".join(final_chars)


def _compute_metrics(seq: str) -> Dict[str, float]:
    length = len(seq)
    n_count = seq.count("N")
    gc_count = seq.count("G") + seq.count("C") if length else 0
    denom = length - n_count if length - n_count > 0 else length
    gc_percent = (gc_count / denom * 100) if denom else 0.0
    n_frac = (n_count / length) if length else 1.0
    return {
        "length": length,
        "n_count": n_count,
        "n_frac": n_frac,
        "gc_percent": gc_percent,
    }


def _load_metadata_organisms(metadata_path: Path | None) -> dict[str, str]:
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
            accession = payload.get("accession") or payload.get("accessionversion")
            if accession:
                lookup[str(accession)] = organism
            uid = payload.get("uid")
            if uid is not None:
                lookup[str(uid)] = organism
    return lookup


def _resolve_organism(
    accession: str | None,
    record_id: str,
    header: str,
    lookup: dict[str, str],
    unknown_warning_emitted: bool,
    logger,
) -> tuple[str, bool]:
    organism = None
    if accession and accession in lookup:
        organism = lookup[accession]
    elif record_id in lookup:
        organism = lookup[record_id]
    if not organism:
        organism = _infer_organism_from_header(header)
        if organism == "UNKNOWN" and not unknown_warning_emitted:
            logger.warning("Some sequences have unknown organism; defaulting to 'UNKNOWN'.")
            unknown_warning_emitted = True
    return organism, unknown_warning_emitted


def _infer_organism_from_header(header: str) -> str:
    if " " not in header:
        return "UNKNOWN"
    remainder = header.split(" ", 1)[1].strip()
    if not remainder:
        return "UNKNOWN"
    lower = remainder.lower()
    if " voucher" in lower:
        idx = lower.index(" voucher")
        remainder = remainder[:idx].strip()
    return remainder if remainder else "UNKNOWN"


def _filter_metadata(metadata_path: Path, output_path: Path, kept_accessions: set[str], kept_ids: set[str]) -> None:
    with metadata_path.open("r", encoding="utf-8") as source, output_path.open("w", encoding="utf-8") as target:
        for line in source:
            line = line.strip()
            if not line:
                continue
            payload = json.loads(line)
            accession = payload.get("accession") or payload.get("accessionversion")
            uid = str(payload.get("uid")) if payload.get("uid") is not None else None
            if (accession and accession in kept_accessions) or (uid and uid in kept_ids):
                target.write(json.dumps(payload) + "\n")


def _log_summary(
    logger,
    total: int,
    kept: int,
    reason_counter: Counter[str],
    fasta_path: Path,
    metadata_path: Path | None,
    report_path: Path | None,
    raw_organism_counts: Counter[str],
) -> None:
    dropped = total - kept
    logger.info("Filter processed %s sequences: kept=%s, dropped=%s", total, kept, dropped)
    if reason_counter:
        top_reasons = reason_counter.most_common(3)
        logger.info("Top exclusion reasons: %s", ", ".join(f"{reason}:{count}" for reason, count in top_reasons))
    if raw_organism_counts:
        top_orgs = raw_organism_counts.most_common(3)
        logger.info(
            "Top organisms (pre-filter): %s",
            ", ".join(f"{org}:{count}" for org, count in top_orgs),
        )
    logger.info("Filtered FASTA -> %s", fasta_path)
    if metadata_path:
        logger.info("Filtered metadata -> %s", metadata_path)
    if report_path:
        logger.info("Filter report -> %s", report_path)
