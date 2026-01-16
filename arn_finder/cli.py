"""Command-line interface for ARN Finder."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Any, Callable

from .clustering.similarity import ClusterRequest, cluster_sequences
from .config import (
    CLUSTER_PLACEHOLDER,
    CLUSTER_DEFAULTS,
    CONSENSUS_DEFAULTS,
    CONSENSUS_PLACEHOLDER,
    DEFAULT_FILTERED_DIR,
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
from .io import ensure_dir, now_iso, write_json
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


def main(argv: list[str] | None = None) -> int:
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
