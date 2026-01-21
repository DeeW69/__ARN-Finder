"""Command-line interface for ARN Finder."""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Callable

import requests
from dotenv import load_dotenv

from .blast_client import (
    BlastClientError,
    BlastJobStatus,
    check as blast_check,
    fetch as blast_fetch,
    submit as blast_submit,
)
from .blast_parser import parse_blast_json
from .clustering.similarity import ClusterRequest, cluster_sequences
from .config import (
    CLUSTER_PLACEHOLDER,
    CLUSTER_DEFAULTS,
    CONSENSUS_DEFAULTS,
    CONSENSUS_PLACEHOLDER,
    BLAST_DEFAULTS,
    DEFAULT_FILTERED_DIR,
    DEFAULT_CONSENSUS_DIR,
    DEFAULT_RAW_DIR,
    FETCH_DEFAULTS,
    FILTER_DEFAULTS,
    FILTER_PLACEHOLDER,
    MOTIF_DEFAULTS,
    MOTIF_PLACEHOLDER,
)
from .consensus.runner import ConsensusRequest, build_consensus
from .fetch.ncbi_rna import FetchRequest, fetch_sequences
from .filters.quality import FilterRequest, filter_sequences
from .io import ensure_dir, now_iso, read_fasta_records, write_csv, write_json, write_lines
from .io.paths import sha256_digest
from .logging_utils import configure_logging, get_logger
from .motifs.kmer_scan import MotifRequest, compute_motifs

Handler = Callable[[argparse.Namespace], int]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="arnfinder",
        description="ARN Finder - exploratory RNA/DNA toolkit (MVP stubs).",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    _add_fetch_parser(subparsers)
    _add_filter_parser(subparsers)
    _add_motif_parser(subparsers)
    _add_cluster_parser(subparsers)
    _add_consensus_parser(subparsers)
    _add_blast_parser(subparsers)
    return parser


def _add_fetch_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser("fetch", help="Fetch public sequences via NCBI Entrez.")
    parser.add_argument("--query", required=True, help="NCBI Entrez query.")
    parser.add_argument(
        "--limit",
        type=int,
        default=FETCH_DEFAULTS.limit,
        help=f"Maximum sequences to fetch (default: {FETCH_DEFAULTS.limit}).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=FETCH_DEFAULTS.out_dir,
        help="Directory for raw data outputs.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=FETCH_DEFAULTS.cache_dir,
        help="Directory for local fetch cache.",
    )
    parser.add_argument("--db", default=FETCH_DEFAULTS.db, help="NCBI database to query (default: nuccore).")
    parser.add_argument("--rettype", default=FETCH_DEFAULTS.rettype, help="Rettype for efetch (default: fasta).")
    parser.add_argument("--retmode", default=FETCH_DEFAULTS.retmode, help="Retmode for efetch (default: text).")
    parser.add_argument(
        "--email",
        help="Contact email for NCBI (default: read from NCBI_EMAIL env if unset).",
    )
    parser.add_argument("--tool", default=FETCH_DEFAULTS.tool, help="Tool identifier reported to NCBI.")
    parser.add_argument(
        "--api-key",
        help="NCBI API key (default: read from NCBI_API_KEY env if unset).",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=None,
        help="Seconds to wait between requests (default: 0.34 or 0.1 when API key provided).",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=FETCH_DEFAULTS.timeout,
        help="HTTP timeout in seconds.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=FETCH_DEFAULTS.retries,
        help="Number of retries for transient network errors.",
    )
    parser.set_defaults(handler=_handle_fetch)


def _add_filter_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser("filter", help="Filter sequences by quality metrics.")
    parser.add_argument(
        "--in-fasta",
        type=Path,
        required=True,
        help="Input FASTA file path.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=FILTER_DEFAULTS.out_dir,
        help="Directory for filtered outputs.",
    )
    parser.add_argument(
        "--min-len",
        type=int,
        default=FILTER_DEFAULTS.min_len,
        help=f"Minimum sequence length (default: {FILTER_DEFAULTS.min_len}).",
    )
    parser.add_argument(
        "--max-n-frac",
        type=float,
        default=FILTER_DEFAULTS.max_n_fraction,
        help=f"Maximum fraction of ambiguous bases (default: {FILTER_DEFAULTS.max_n_fraction}).",
    )
    parser.add_argument(
        "--alphabet",
        choices=["DNA", "RNA", "AUTO"],
        default="AUTO",
        help="Alphabet assumption for sequences (default: AUTO).",
    )
    parser.add_argument(
        "--dedupe",
        action="store_true",
        help="Remove duplicate sequences, keeping first occurrence.",
    )
    parser.add_argument(
        "--in-metadata",
        type=Path,
        help="Metadata JSONL file to filter alongside sequences.",
    )
    parser.add_argument(
        "--write-report",
        dest="write_report",
        action="store_true",
        default=True,
        help="Write per-record CSV report (default: enabled).",
    )
    parser.add_argument(
        "--no-report",
        dest="write_report",
        action="store_false",
        help="Disable per-record CSV report generation.",
    )
    parser.add_argument(
        "--max-per-organism",
        type=int,
        default=0,
        help="Maximum sequences to keep per organism (0 disables).",
    )
    parser.set_defaults(handler=_handle_filter)


def _add_motif_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser("motifs", help="Extract k-mer motifs.")
    parser.add_argument(
        "--in-fasta",
        type=Path,
        required=True,
        help="Input FASTA file (typically filtered sequences).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=MOTIF_DEFAULTS.out_dir,
        help="Directory for motif outputs.",
    )
    parser.add_argument(
        "--k",
        type=int,
        default=MOTIF_DEFAULTS.k,
        help=f"k-mer length (default: {MOTIF_DEFAULTS.k}).",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=200,
        help="Number of motifs to retain (0 disables).",
    )
    parser.add_argument(
        "--min-count",
        type=int,
        default=1,
        help="Minimum global count for a motif to be retained.",
    )
    parser.add_argument(
        "--ignore-N",
        dest="ignore_n",
        action="store_true",
        default=True,
        help="Ignore k-mers containing 'N'.",
    )
    parser.add_argument(
        "--allow-N",
        dest="ignore_n",
        action="store_false",
        help="Include k-mers containing 'N'.",
    )
    parser.add_argument(
        "--canonical",
        action="store_true",
        help="Group DNA k-mers with their reverse complement.",
    )
    parser.add_argument(
        "--per-record",
        dest="per_record",
        action="store_true",
        default=True,
        help="Produce per-record motif counts.",
    )
    parser.add_argument(
        "--no-per-record",
        dest="per_record",
        action="store_false",
        help="Disable per-record motif counts.",
    )
    parser.add_argument(
        "--alphabet",
        choices=["DNA", "RNA", "AUTO"],
        default="AUTO",
        help="Alphabet assumption for canonicalization.",
    )
    parser.set_defaults(handler=_handle_motifs)


def _add_cluster_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser("cluster", help="Cluster sequences via k-mer Jaccard similarity.")
    parser.add_argument(
        "--in-fasta",
        type=Path,
        required=True,
        help="Filtered FASTA file path.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=CLUSTER_DEFAULTS.out_dir,
        help="Directory for cluster outputs.",
    )
    parser.add_argument(
        "--k",
        type=int,
        default=CLUSTER_DEFAULTS.k,
        help=f"k-mer size for similarity (default: {CLUSTER_DEFAULTS.k}).",
    )
    parser.add_argument(
        "--min-sim",
        type=float,
        default=CLUSTER_DEFAULTS.min_similarity,
        help=f"Minimum Jaccard similarity threshold (default: {CLUSTER_DEFAULTS.min_similarity}).",
    )
    parser.add_argument(
        "--ignore-N",
        dest="ignore_n",
        action="store_true",
        default=True,
        help="Ignore k-mers containing 'N'.",
    )
    parser.add_argument(
        "--allow-N",
        dest="ignore_n",
        action="store_false",
        help="Include k-mers containing 'N'.",
    )
    parser.add_argument(
        "--canonical",
        action="store_true",
        help="Group DNA k-mers with their reverse complement.",
    )
    parser.add_argument(
        "--alphabet",
        choices=["DNA", "RNA", "AUTO"],
        default="AUTO",
        help="Alphabet assumption for canonicalization.",
    )
    parser.add_argument(
        "--max-pairs",
        type=int,
        default=200_000,
        help="Maximum number of pairwise comparisons (guardrail).",
    )
    parser.add_argument(
        "--in-metadata",
        type=Path,
        help="Metadata JSONL for enrichment.",
    )
    parser.set_defaults(handler=_handle_cluster)


def _add_consensus_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser("consensus", help="Build exploratory cluster consensus sequences.")
    parser.add_argument(
        "--clusters",
        type=Path,
        required=True,
        help="clusters.jsonl produced by `arnfinder cluster`.",
    )
    parser.add_argument(
        "--in-fasta",
        type=Path,
        required=True,
        help="Filtered FASTA file path.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=CONSENSUS_DEFAULTS.out_dir,
        help="Directory for consensus outputs.",
    )
    parser.add_argument(
        "--min-cluster-size",
        type=int,
        default=CONSENSUS_DEFAULTS.min_cluster_size,
        help="Minimum cluster size to attempt consensus.",
    )
    parser.add_argument(
        "--method",
        choices=["simple", "kmer"],
        default=CONSENSUS_DEFAULTS.method,
        help="Consensus strategy (default: kmer).",
    )
    parser.add_argument(
        "--k",
        type=int,
        default=CONSENSUS_DEFAULTS.k,
        help="k-mer size used by kmer method.",
    )
    parser.add_argument(
        "--min-overlap",
        type=int,
        default=CONSENSUS_DEFAULTS.min_overlap,
        help="Minimum overlap length when merging sequences.",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=CONSENSUS_DEFAULTS.min_identity,
        help="Minimum identity across the overlap.",
    )
    parser.add_argument(
        "--max-n-frac",
        type=float,
        default=CONSENSUS_DEFAULTS.max_n_frac,
        help="Maximum allowed fraction of N in final consensus.",
    )
    parser.add_argument(
        "--canonical",
        dest="canonical",
        action="store_true",
        default=True,
        help="Use canonical reverse-complement heuristics where applicable.",
    )
    parser.add_argument(
        "--no-canonical",
        dest="canonical",
        action="store_false",
        help="Disable canonical heuristics.",
    )
    parser.add_argument(
        "--in-metadata",
        type=Path,
        help="Optional metadata JSONL for enrichment.",
    )
    parser.set_defaults(handler=_handle_consensus)


def _add_blast_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser("blast", help="Optional BLAST validation via NCBI.")
    parser.add_argument(
        "--in-fasta",
        type=Path,
        help="Input FASTA (defaults to consensus.fasta if present, else filtered_sequences.fasta).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=BLAST_DEFAULTS.out_dir,
        help="Directory for BLAST outputs.",
    )
    parser.add_argument(
        "--program",
        default=None,
        help="BLAST program (default: BLAST_PROGRAM env or blastn).",
    )
    parser.add_argument(
        "--db",
        default=None,
        help="BLAST database (default: BLAST_DB env or nt).",
    )
    parser.add_argument(
        "--format",
        choices=["JSON2", "Tabular", "XML"],
        default=BLAST_DEFAULTS.format_type,
        help="Output format retrieved from BLAST (default: JSON2).",
    )
    parser.add_argument(
        "--max-hits",
        type=int,
        default=BLAST_DEFAULTS.max_hits,
        help="Maximum number of hits to retain per query.",
    )
    parser.add_argument(
        "--poll-seconds",
        type=float,
        default=BLAST_DEFAULTS.poll_seconds,
        help="Sleep interval between BLAST status polls.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=BLAST_DEFAULTS.timeout_seconds,
        help="Maximum seconds to wait per query before giving up.",
    )
    parser.add_argument(
        "--rate-limit-seconds",
        type=float,
        default=BLAST_DEFAULTS.rate_limit_seconds,
        help="Minimum wait between HTTP requests to respect NCBI rate limits.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        help="Optional directory for caching BLAST results by sequence.",
    )
    parser.set_defaults(handler=_handle_blast)


def _handle_fetch(args: argparse.Namespace) -> int:
    logger = get_logger()
    ensure_dir(args.out_dir)
    ensure_dir(args.cache_dir)
    email = args.email or os.environ.get("NCBI_EMAIL")
    api_key = args.api_key or os.environ.get("NCBI_API_KEY")
    sleep_interval = args.sleep if args.sleep is not None else (0.1 if api_key else FETCH_DEFAULTS.sleep)
    request = FetchRequest(
        query=args.query,
        limit=args.limit,
        out_dir=args.out_dir,
        cache_dir=args.cache_dir,
        db=args.db,
        rettype=args.rettype,
        retmode=args.retmode,
        email=email,
        tool=args.tool,
        api_key=api_key,
        sleep=sleep_interval,
        timeout=args.timeout,
        retries=args.retries,
        batch_size=FETCH_DEFAULTS.batch_size,
    )
    result = fetch_sequences(request, logger)
    logger.info("Fetched %s/%s IDs", len(result.ids), result.total_count)
    logger.info("Sequences written to %s", result.fasta_path)
    logger.info("Metadata written to %s", result.metadata_path)
    logger.info("Manifest written to %s", result.manifest_path)
    return 0


def _handle_filter(args: argparse.Namespace) -> int:
    _require_input(args.in_fasta)
    if args.in_metadata:
        _require_input(args.in_metadata)
    request = FilterRequest(
        input_fasta=args.in_fasta,
        out_dir=args.out_dir,
        min_len=args.min_len,
        max_n_fraction=args.max_n_frac,
        alphabet=args.alphabet,
        dedupe=args.dedupe,
        metadata_path=args.in_metadata,
        write_report=args.write_report,
        max_per_organism=args.max_per_organism,
    )
    result = filter_sequences(request, get_logger())
    logger = get_logger()
    logger.info("Filter manifest -> %s", result.manifest_path)
    return 0


def _handle_motifs(args: argparse.Namespace) -> int:
    _require_input(args.in_fasta)
    request = MotifRequest(
        input_fasta=args.in_fasta,
        out_dir=args.out_dir,
        k=args.k,
        top=args.top,
        min_count=args.min_count,
        ignore_n=args.ignore_n,
        canonical=args.canonical,
        per_record=args.per_record,
        alphabet=args.alphabet,
    )
    compute_motifs(request, get_logger())
    return 0


def _handle_cluster(args: argparse.Namespace) -> int:
    _require_input(args.in_fasta)
    if args.in_metadata:
        _require_input(args.in_metadata)
    request = ClusterRequest(
        input_fasta=args.in_fasta,
        out_dir=args.out_dir,
        k=args.k,
        min_similarity=args.min_sim,
        ignore_n=args.ignore_n,
        canonical=args.canonical,
        alphabet=args.alphabet,
        max_pairs=args.max_pairs,
        metadata_path=args.in_metadata,
    )
    cluster_sequences(request, get_logger())
    return 0


def _handle_consensus(args: argparse.Namespace) -> int:
    _require_input(args.clusters)
    _require_input(args.in_fasta)
    if args.in_metadata:
        _require_input(args.in_metadata)
    request = ConsensusRequest(
        clusters_path=args.clusters,
        fasta_path=args.in_fasta,
        out_dir=args.out_dir,
        min_cluster_size=args.min_cluster_size,
        method=args.method,
        k=args.k,
        min_overlap=args.min_overlap,
        min_identity=args.min_identity,
        max_n_frac=args.max_n_frac,
        canonical=args.canonical,
        metadata_path=args.in_metadata,
    )
    build_consensus(request, get_logger())
    return 0


def _handle_blast(args: argparse.Namespace) -> int:
    logger = get_logger()
    input_fasta = args.in_fasta or _default_blast_input()
    _require_input(input_fasta)
    ensure_dir(args.out_dir)
    raw_dir = ensure_dir(args.out_dir / "raw")
    cache = BlastCache(args.cache_dir)

    program = args.program or os.environ.get("BLAST_PROGRAM") or BLAST_DEFAULTS.program
    db = args.db or os.environ.get("BLAST_DB") or BLAST_DEFAULTS.db
    format_type = args.format.upper()
    max_hits = args.max_hits
    email = os.environ.get("NCBI_EMAIL")
    if not email:
        logger.error("NCBI_EMAIL environment variable is required for BLAST submissions.")
        return 1
    tool = os.environ.get("NCBI_TOOL") or "arnfinder"
    api_key = os.environ.get("NCBI_API_KEY")
    sequences = list(read_fasta_records(input_fasta))
    if not sequences:
        logger.warning("No sequences found in %s", input_fasta)
    limiter = RateLimiter(args.rate_limit_seconds)
    session = requests.Session()

    hits_records: list[str] = []
    summary_rows: list[list[str]] = []
    success = 0
    failures = 0

    def _record(query_id: str, rid: str | None, top_hits: list[dict[str, Any]], raw_path: Path | None, error: str | None = None) -> None:
        entry: dict[str, Any] = {
            "query_id": query_id,
            "rid": rid,
            "top_hits": top_hits,
            "raw_path": str(raw_path) if raw_path else None,
        }
        if error:
            entry["error"] = error
        hits_records.append(json.dumps(entry))
        summary_rows.append(_blast_summary_row(query_id, top_hits))

    for idx, (header, sequence) in enumerate(sequences, start=1):
        query_id = _resolve_query_id(header, idx)
        cache_key = _blast_cache_key(sequence, program, db, format_type, max_hits)
        cached = cache.load(cache_key)
        if cached:
            rid = cached.get("rid")
            top_hits = cached.get("top_hits", [])
            raw_text: str = cached.get("raw_text", "")
            raw_path = _write_raw_result(raw_dir, idx, query_id, format_type, raw_text)
            _record(query_id, rid, top_hits, raw_path)
            success += 1
            logger.debug("Cache hit for %s (RID=%s)", query_id, rid)
            continue

        rid: str | None = None
        try:
            limiter.wait()
            rid, rtoe = blast_submit(
                sequence=sequence,
                program=program,
                db=db,
                session=session,
                email=email,
                tool=tool,
                api_key=api_key,
                timeout=args.timeout_seconds,
                hitlist_size=max_hits,
            )
            logger.info("Submitted BLAST job for %s (RID=%s, RTOE=%s)", query_id, rid, rtoe)
            deadline = time.monotonic() + args.timeout_seconds
            status = BlastJobStatus.WAITING
            while time.monotonic() < deadline:
                limiter.wait()
                status = blast_check(rid, session=session, api_key=api_key)
                if status == BlastJobStatus.READY:
                    break
                if status == BlastJobStatus.FAILED:
                    raise BlastClientError(f"BLAST job failed for RID {rid}")
                time.sleep(args.poll_seconds)
            if status != BlastJobStatus.READY:
                raise TimeoutError(f"BLAST job for RID {rid} timed out after {args.timeout_seconds} seconds.")
            limiter.wait()
            json_raw = blast_fetch(rid, format_type="JSON2", session=session, api_key=api_key, timeout=args.timeout_seconds)
            top_hits = parse_blast_json(json_raw, max_hits, len(sequence))
            if format_type == "JSON2":
                raw_text = json_raw
            else:
                limiter.wait()
                raw_text = blast_fetch(
                    rid,
                    format_type=format_type,
                    session=session,
                    api_key=api_key,
                    timeout=args.timeout_seconds,
                )
            raw_path = _write_raw_result(raw_dir, idx, query_id, format_type, raw_text)
            _record(query_id, rid, top_hits, raw_path)
            success += 1
            cache.save(
                cache_key,
                {
                    "query_id": query_id,
                    "rid": rid,
                    "format_type": format_type,
                    "program": program,
                    "database": db,
                    "max_hits": max_hits,
                    "top_hits": top_hits,
                    "timestamp": now_iso(),
                },
                raw_text,
            )
        except (requests.RequestException, BlastClientError, TimeoutError) as exc:
            failures += 1
            logger.error("BLAST failed for %s: %s", query_id, exc)
            _record(query_id, rid, [], None, error=str(exc))

    hits_path = args.out_dir / "blast_hits.jsonl"
    summary_path = args.out_dir / "blast_summary.csv"
    manifest_path = args.out_dir / "blast_manifest.json"
    write_lines(hits_path, hits_records)
    write_csv(
        summary_path,
        ["query_id", "accession", "pident", "coverage", "evalue", "bitscore"],
        summary_rows,
    )
    manifest = {
        "timestamp": now_iso(),
        "input_fasta": str(input_fasta),
        "out_dir": str(args.out_dir),
        "raw_dir": str(raw_dir),
        "program": program,
        "database": db,
        "format": format_type,
        "max_hits": max_hits,
        "poll_seconds": args.poll_seconds,
        "timeout_seconds": args.timeout_seconds,
        "rate_limit_seconds": args.rate_limit_seconds,
        "cache_dir": str(args.cache_dir) if args.cache_dir else None,
        "cache_enabled": bool(args.cache_dir),
        "email": email,
        "tool": tool,
        "api_key_provided": bool(api_key),
        "total_queries": len(sequences),
        "succeeded": success,
        "failed": failures,
    }
    write_json(manifest_path, manifest)
    logger.info(
        "BLAST complete for %s queries (%s ok / %s failed). Results -> %s",
        len(sequences),
        success,
        failures,
        args.out_dir,
    )
    return 0


def _write_placeholder(out_dir: Path, filename: str, payload: dict[str, Any]) -> Path:
    ensure_dir(out_dir)
    path = out_dir / filename
    data = dict(payload)
    data["timestamp"] = now_iso()
    write_json(path, data)
    return path


def _require_input(path: Path) -> None:
    if not Path(path).exists():
        raise FileNotFoundError(f"Input file not found: {path}")


def _default_blast_input() -> Path:
    consensus = DEFAULT_CONSENSUS_DIR / "consensus.fasta"
    if consensus.exists():
        return consensus
    return DEFAULT_FILTERED_DIR / "filtered_sequences.fasta"


def _resolve_query_id(header: str, index: int) -> str:
    if header:
        return header.split()[0]
    return f"query_{index:04d}"


def _blast_cache_key(sequence: str, program: str, db: str, format_type: str, max_hits: int) -> str:
    payload = f"{program}|{db}|{format_type}|{max_hits}|{sequence}"
    return sha256_digest(payload)


def _safe_filename(value: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in value)
    return cleaned[:100] or "query"


def _raw_extension(format_type: str) -> str:
    if format_type.upper() == "TABULAR":
        return "tsv"
    if format_type.upper() == "XML":
        return "xml"
    return "json"


def _write_raw_result(raw_dir: Path, index: int, query_id: str, format_type: str, content: str) -> Path:
    ensure_dir(raw_dir)
    ext = _raw_extension(format_type)
    filename = f"{index:04d}_{_safe_filename(query_id)}.{ext}"
    path = raw_dir / filename
    path.write_text(content, encoding="utf-8")
    return path


def _format_float(value: Any) -> str:
    try:
        return f"{float(value):.3f}"
    except (TypeError, ValueError):
        return ""


def _blast_summary_row(query_id: str, hits: list[dict[str, Any]]) -> list[str]:
    if hits:
        best = hits[0]
        return [
            query_id,
            str(best.get("accession", "")),
            _format_float(best.get("pident")),
            _format_float(best.get("coverage")),
            str(best.get("evalue", "")),
            str(best.get("bitscore", "")),
        ]
    return [query_id, "", "", "", "", ""]


class RateLimiter:
    """Enforce a minimum delay between HTTP operations."""

    def __init__(self, interval: float) -> None:
        self.interval = max(interval, 0.0)
        self._last_call = 0.0

    def wait(self) -> None:
        if self.interval <= 0:
            return
        now = time.monotonic()
        if self._last_call == 0.0:
            self._last_call = now
            return
        delta = now - self._last_call
        if delta < self.interval:
            time.sleep(self.interval - delta)
        self._last_call = time.monotonic()


class BlastCache:
    """Simple file-based cache for BLAST results."""

    def __init__(self, root: Path | None) -> None:
        self.root = root
        if self.root:
            ensure_dir(self.root)

    def load(self, key: str) -> dict[str, Any] | None:
        if not self.root:
            return None
        meta_path = self.root / f"{key}.json"
        raw_path = self.root / f"{key}.raw"
        if not meta_path.exists() or not raw_path.exists():
            return None
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        meta["raw_text"] = raw_path.read_text(encoding="utf-8")
        return meta

    def save(self, key: str, metadata: dict[str, Any], raw_text: str) -> None:
        if not self.root:
            return
        ensure_dir(self.root)
        meta_path = self.root / f"{key}.json"
        raw_path = self.root / f"{key}.raw"
        meta_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        raw_path.write_text(raw_text, encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    load_dotenv()
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(verbose=args.verbose)
    logger = get_logger()
    handler: Handler = args.handler

    try:
        return handler(args)
    except FileNotFoundError as exc:
        logger.error(str(exc))
        return 1
    except Exception:  # pragma: no cover - safety net
        logger.exception("Unexpected error")
        return 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
