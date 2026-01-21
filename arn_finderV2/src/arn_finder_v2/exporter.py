"""Data export helpers for ML-ready parquet datasets."""

from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd

logger = logging.getLogger("arnfinderv2.export")


@dataclass
class ExportPaths:
    features_parquet: Path
    motifs_parquet: Optional[Path]
    joined_parquet: Path
    manifest_json: Path


def export_ml_ready(
    features_csv: Path,
    metadata_jsonl: Path,
    out_dir: Path,
    motifs_csv: Optional[Path] = None,
    cli_version: Optional[str] = None,
) -> ExportPaths:
    """Load CSV/JSONL inputs, join metadata, and write parquet exports."""
    features_csv = _ensure_file(features_csv)
    metadata_jsonl = _ensure_file(metadata_jsonl)
    if motifs_csv:
        motifs_csv = _ensure_file(motifs_csv)

    out_dir.mkdir(parents=True, exist_ok=True)
    features = _load_features(features_csv)
    metadata = _load_metadata(metadata_jsonl)
    joined, join_keys = _join_features_metadata(features, metadata)

    features_path = out_dir / "features.parquet"
    joined_path = out_dir / "joined_dataset.parquet"
    features.to_parquet(features_path, index=False)
    joined.to_parquet(joined_path, index=False)

    motifs_path = None
    motifs_rows = None
    if motifs_csv:
        motifs_df = pd.read_csv(motifs_csv)
        motifs_path = out_dir / "motifs_by_record.parquet"
        motifs_df.to_parquet(motifs_path, index=False)
        motifs_rows = len(motifs_df)

    manifest_path = out_dir / "export_manifest.json"
    manifest = {
        "params": {
            "features_csv": str(features_csv),
            "metadata_jsonl": str(metadata_jsonl),
            "motifs_by_record_csv": str(motifs_csv) if motifs_csv else None,
            "out_dir": str(out_dir),
        },
        "join_keys": join_keys,
        "counts": {
            "features_rows": int(len(features)),
            "metadata_rows": int(len(metadata)),
            "joined_rows": int(len(joined)),
            "motifs_rows": int(motifs_rows) if motifs_rows is not None else None,
        },
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    if cli_version:
        manifest["version"] = cli_version
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    logger.info("Export parquet terminé dans %s", out_dir)
    return ExportPaths(
        features_parquet=features_path,
        motifs_parquet=motifs_path,
        joined_parquet=joined_path,
        manifest_json=manifest_path,
    )


def _ensure_file(path: Path) -> Path:
    path = path.expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Fichier introuvable : {path}")
    return path


def _load_features(path: Path) -> pd.DataFrame:
    logger.info("Chargement des features depuis %s", path)
    return pd.read_csv(path)


def _load_metadata(path: Path) -> pd.DataFrame:
    logger.info("Chargement des metadata depuis %s", path)
    return pd.read_json(path, lines=True)


def _join_features_metadata(features: pd.DataFrame, metadata: pd.DataFrame) -> Tuple[pd.DataFrame, Optional[Tuple[str, str]]]:
    candidates = [
        ("record_id", "record_id"),
        ("record_id", "accession"),
        ("accession", "accession"),
    ]
    for left, right in candidates:
        if left in features.columns and right in metadata.columns:
            joined = features.merge(metadata, how="left", left_on=left, right_on=right, suffixes=("", "_meta"))
            logger.info("Jointure metadata via %s -> %s", left, right)
            return joined, (left, right)

    logger.warning("Aucune colonne compatible pour joindre les metadata. Le dataset reste inchangé.")
    return features.copy(), None
