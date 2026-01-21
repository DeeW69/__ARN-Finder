"""Command-line interface for ARN Finder V2."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import time
from dataclasses import asdict
from pathlib import Path
from typing import Iterable, Optional

from . import __version__
from .blast_client import BlastClient, BlastClientError
from .blast_parser import BlastParseError, parse_json2_top_hits
from .config import collect_runtime_config, load_dotenv
from .exporter import export_ml_ready
from .fasta import FastaRecord, read_fasta

logger = logging.getLogger("arnfinderv2")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="arnfinderv2",
        description="Next generation tooling for ARN Finder experiments.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", metavar="command")

    blast_parser = subparsers.add_parser(
        "blast",
        help="Exécuter BLAST (NCBI) sur un fichier FASTA.",
        description="Soumet des séquences FASTA à BLAST (NCBI) et produit des résumés JSONL/CSV.",
    )
    blast_parser.add_argument("--in-fasta", required=True, help="Chemin vers un FASTA d'entrée.")
    blast_parser.add_argument("--out-dir", default="data/blast", help="Répertoire de sortie (défaut data/blast).")
    blast_parser.add_argument("--program", default=None, help="Programme BLAST (défaut ENV ou blastn).")
    blast_parser.add_argument("--db", default=None, help="Base BLAST (défaut ENV ou nt).")
    blast_parser.add_argument(
        "--format",
        default="JSON2",
        choices=["JSON2", "TABULAR", "XML"],
        type=lambda value: value.upper(),
        help="Format BLAST à récupérer (JSON2 recommandé).",
    )
    blast_parser.add_argument("--max-hits", type=int, default=10, help="Nombre maximum de hits par séquence.")
    blast_parser.add_argument("--poll-seconds", type=int, default=10, help="Intervalle entre deux polls (s).")
    blast_parser.add_argument("--timeout-seconds", type=int, default=600, help="Timeout global pour chaque RID (s).")
    blast_parser.add_argument(
        "--rate-limit-seconds",
        type=float,
        default=1.0,
        help="Délai minimal entre requêtes HTTP pour respecter les règles NCBI.",
    )
    blast_parser.add_argument(
        "--cache-dir",
        default="data/cache/blast",
        help="Répertoire de cache pour stocker les réponses BLAST brutes.",
    )
    blast_parser.set_defaults(func=_command_blast)

    export_parser = subparsers.add_parser(
        "export",
        help="Exporter les features/metadata en Parquet prêt pour le ML.",
        description="Charge les CSV/JSONL produits par ARN Finder et génère des Parquet consolidés.",
    )
    export_parser.add_argument(
        "--features-csv",
        default="data/filtered/filtered_features.csv",
        help="Chemin vers le CSV des features filtrées.",
    )
    export_parser.add_argument(
        "--metadata-jsonl",
        default="data/raw/metadata.jsonl",
        help="Chemin vers le JSONL metadata (une entrée par record).",
    )
    export_parser.add_argument(
        "--motifs-by-record-csv",
        default=None,
        help="CSV motifs par record (optionnel).",
    )
    export_parser.add_argument(
        "--out-dir",
        default="data/exports",
        help="Répertoire de sortie pour les Parquet exportés.",
    )
    export_parser.set_defaults(func=_command_export)

    return parser


def main(argv: Optional[Iterable[str]] = None) -> int:
    load_dotenv()  # Best-effort load before parsing to let args inspect env if needed.
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)

    if not hasattr(args, "func"):
        parser.print_help()
        return 0

    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    return args.func(args)


def _command_blast(args: argparse.Namespace) -> int:
    config = collect_runtime_config()
    missing = [k for k in ("ncbi_email", "ncbi_tool") if not getattr(config, k)]
    if missing:
        raise SystemExit(f"Variables manquantes dans le .env : {', '.join(missing)}")

    program = (args.program or config.blast_program or "blastn").lower()
    db = args.db or config.blast_db or "nt"
    format_type = args.format
    out_dir = Path(args.out_dir).resolve()
    raw_dir = out_dir / "raw"
    cache_dir = Path(args.cache_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    records = read_fasta(Path(args.in_fasta), uppercase=True, convert_u_to_t=program == "blastn")
    logger.info("FASTA chargé avec %d séquences.", len(records))

    client = BlastClient(
        email=config.ncbi_email,
        tool=config.ncbi_tool,
        api_key=config.ncbi_api_key,
        rate_limit_seconds=args.rate_limit_seconds,
    )

    hits_path = out_dir / "blast_hits.jsonl"
    summary_path = out_dir / "blast_summary.csv"
    manifest_path = out_dir / "blast_manifest.json"

    success = 0
    failed = 0
    start_ts = time.time()

    with hits_path.open("w", encoding="utf-8") as hits_fh, summary_path.open("w", newline="", encoding="utf-8") as summary_fh:
        summary_writer = csv.DictWriter(
            summary_fh,
            fieldnames=["query_id", "rid", "top_accession", "top_organism", "pident", "qcov", "evalue", "bitscore"],
        )
        summary_writer.writeheader()

        for record in records:
            result = _process_sequence(
                record=record,
                client=client,
                program=program,
                db=db,
                format_type=format_type,
                max_hits=args.max_hits,
                poll_seconds=args.poll_seconds,
                timeout_seconds=args.timeout_seconds,
                cache_dir=cache_dir,
                raw_dir=raw_dir,
            )
            if result["status"] == "ok":
                success += 1
            else:
                failed += 1

            hits_fh.write(json.dumps(result["hit_entry"], ensure_ascii=False) + "\n")

            top_hit = result["hit_entry"]["top_hits"][0] if result["hit_entry"]["top_hits"] else {}
            summary_writer.writerow(
                {
                    "query_id": record.identifier,
                    "rid": result["hit_entry"]["rid"],
                    "top_accession": top_hit.get("accession"),
                    "top_organism": top_hit.get("organism"),
                    "pident": top_hit.get("pident"),
                    "qcov": top_hit.get("qcov"),
                    "evalue": top_hit.get("evalue"),
                    "bitscore": top_hit.get("bitscore"),
                }
            )

    end_ts = time.time()
    manifest = {
        "params": {
            "program": program,
            "db": db,
            "format": format_type,
            "max_hits": args.max_hits,
            "poll_seconds": args.poll_seconds,
            "timeout_seconds": args.timeout_seconds,
            "rate_limit_seconds": args.rate_limit_seconds,
            "cache_dir": str(cache_dir),
            "out_dir": str(out_dir),
            "in_fasta": str(Path(args.in_fasta).resolve()),
        },
        "timestamps": {"started_at": _to_iso(start_ts), "finished_at": _to_iso(end_ts)},
        "nb_sequences": len(records),
        "ok": success,
        "failed": failed,
        "version": __version__,
        "config": _mask_config(config),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    logger.info("Traitement terminé : %d succès / %d échecs.", success, failed)
    return 0 if failed == 0 else 1


def _command_export(args: argparse.Namespace) -> int:
    try:
        export_ml_ready(
            features_csv=Path(args.features_csv),
            metadata_jsonl=Path(args.metadata_jsonl),
            motifs_csv=Path(args.motifs_by_record_csv) if args.motifs_by_record_csv else None,
            out_dir=Path(args.out_dir),
            cli_version=__version__,
        )
    except FileNotFoundError as exc:
        logger.error("%s", exc)
        return 2
    except Exception as exc:  # noqa: BLE001
        logger.error("Export ML-ready échoué : %s", exc)
        return 1

    logger.info("Export ML-ready terminé.")
    return 0


def _process_sequence(
    record: FastaRecord,
    client: BlastClient,
    program: str,
    db: str,
    format_type: str,
    max_hits: int,
    poll_seconds: int,
    timeout_seconds: int,
    cache_dir: Path,
    raw_dir: Path,
):
    rid = None
    try:
        rid = client.submit(query_seq=record.sequence, program=program, db=db)
        status = client.poll(rid=rid, poll_seconds=poll_seconds, timeout_seconds=timeout_seconds)
        if status.state == "FAILED":
            raise BlastClientError(status.message or "BLAST a renvoyé FAILED.")
        if status.state == "TIMEOUT":
            raise BlastClientError("Timeout sur le poll BLAST.")

        raw_text, _ = _load_or_fetch_result(client, rid, format_type, cache_dir)
        raw_file = _write_raw_result(record.identifier, rid, raw_text, raw_dir, format_type)

        top_hits = []
        if format_type.upper() == "JSON2":
            try:
                top_hits = parse_json2_top_hits(raw_text, max_hits=max_hits)
            except BlastParseError as exc:
                logger.warning("Parsing JSON2 impossible pour %s (%s): %s", record.identifier, rid, exc)
        hit_entry = {
            "query_id": record.identifier,
            "rid": rid,
            "program": program,
            "db": db,
            "top_hits": top_hits,
            "raw_result_path": str(raw_file),
        }
        return {"status": "ok", "hit_entry": hit_entry}
    except Exception as exc:  # noqa: BLE001
        logger.warning("Sequence %s impossible à traiter : %s", record.identifier, exc)
        hit_entry = {
            "query_id": record.identifier,
            "rid": rid,
            "program": program,
            "db": db,
            "top_hits": [],
            "raw_result_path": None,
            "error": str(exc),
        }
        return {"status": "failed", "hit_entry": hit_entry}


def _load_or_fetch_result(client: BlastClient, rid: str, format_type: str, cache_dir: Path):
    suffix = _format_suffix(format_type)
    cache_file = cache_dir / f"{rid}.{suffix}"
    if cache_file.exists():
        return cache_file.read_text(encoding="utf-8"), cache_file
    raw_text = client.fetch(rid=rid, format_type=format_type)
    cache_file.write_text(raw_text, encoding="utf-8")
    return raw_text, cache_file


def _write_raw_result(query_id: str, rid: str, content: str, raw_dir: Path, format_type: str) -> Path:
    safe_query = "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in query_id)[:120]
    suffix = _format_suffix(format_type)
    raw_file = raw_dir / f"{safe_query}_{rid}.{suffix}"
    raw_file.write_text(content, encoding="utf-8")
    return raw_file


def _format_suffix(format_type: str) -> str:
    mapping = {"JSON2": "json", "TABULAR": "txt", "XML": "xml"}
    return mapping.get(format_type.upper(), "txt")


def _to_iso(ts: float) -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(ts))


def _mask_config(config) -> dict:
    data = asdict(config)
    api_key = data.get("ncbi_api_key")
    if api_key:
        data["ncbi_api_key"] = f"{api_key[:4]}…"
    return data


if __name__ == "__main__":
    raise SystemExit(main())
