"""Export V3 — Parquet ML-ready avec validation Pandera."""

from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd
import pandera as pa
from pandera import DataFrameSchema, Column, Check

logger = logging.getLogger("arnfinderv3.export")


# ── Schémas Pandera — validation stricte des données ─────────────────────────

FeaturesSchema = DataFrameSchema(
    {
        "record_id": Column(str, nullable=False),
        "length": Column(int, Check.greater_than(0)),
        "gc_frac": Column(float, Check.in_range(0.0, 1.0)),
        "n_frac": Column(float, Check.in_range(0.0, 1.0)),
        "cluster_id": Column(str, nullable=True),
    },
    coerce=True,
)

MetadataSchema = DataFrameSchema(
    {
        "record_id": Column(str, nullable=True),
        "accession": Column(str, nullable=True),
        "organism": Column(str, nullable=True),
    },
    coerce=True,
)


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
    validate: bool = True,
    cli_version: Optional[str] = None,
) -> ExportPaths:
    """
    Charge CSV/JSONL, valide les schémas, joint les metadata,
    et écrit les fichiers Parquet ML-ready.

    Nouveautés V3 :
    - Validation Pandera (optionnelle mais activée par défaut)
    - Rapport d'erreurs structuré si validation échoue
    - Schémas documentés et versionnés
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    features = pd.read_csv(features_csv)
    metadata = pd.read_json(metadata_jsonl, lines=True)

    if validate:
        try:
            features = FeaturesSchema.validate(features, lazy=True)
            metadata = MetadataSchema.validate(metadata, lazy=True)
        except pa.errors.SchemaErrors as exc:
            logger.error("Validation échouée :\n%s", exc.failure_cases)
            raise

    joined, join_keys = _join(features, metadata)

    features_path = out_dir / "features.parquet"
    joined_path = out_dir / "joined_dataset.parquet"
    features.to_parquet(features_path, index=False)
    joined.to_parquet(joined_path, index=False)

    motifs_path = None
    if motifs_csv:
        motifs_df = pd.read_csv(motifs_csv)
        motifs_path = out_dir / "motifs_by_record.parquet"
        motifs_df.to_parquet(motifs_path, index=False)

    manifest = {
        "version": cli_version or "3.0.0",
        "schema_version": "1",
        "params": {
            "features_csv": str(features_csv),
            "metadata_jsonl": str(metadata_jsonl),
            "motifs_csv": str(motifs_csv) if motifs_csv else None,
            "validate": validate,
        },
        "join_keys": list(join_keys) if join_keys else None,
        "counts": {
            "features_rows": len(features),
            "metadata_rows": len(metadata),
            "joined_rows": len(joined),
            "motifs_rows": len(pd.read_csv(motifs_csv)) if motifs_csv else None,
        },
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    manifest_path = out_dir / "export_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    logger.info("Export terminé dans %s", out_dir)
    return ExportPaths(features_path, motifs_path, joined_path, manifest_path)


def _join(
    features: pd.DataFrame, metadata: pd.DataFrame
) -> tuple[pd.DataFrame, Optional[tuple[str, str]]]:
    for left, right in [
        ("record_id", "record_id"),
        ("record_id", "accession"),
        ("accession", "accession"),
    ]:
        if left in features.columns and right in metadata.columns:
            joined = features.merge(metadata, how="left", left_on=left, right_on=right, suffixes=("", "_meta"))
            return joined, (left, right)
    logger.warning("Aucune clé de jointure trouvée.")
    return features.copy(), None
