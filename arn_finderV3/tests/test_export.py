"""Tests d'intégration pour l'export V3 (Parquet + validation Pandera)."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from arn_finder_v3.export import export_ml_ready


@pytest.fixture()
def sample_features(tmp_path: Path) -> Path:
    path = tmp_path / "features.csv"
    pd.DataFrame(
        {
            "record_id": ["SEQ001", "SEQ002"],
            "length": [300, 450],
            "gc_frac": [0.42, 0.55],
            "n_frac": [0.01, 0.0],
            "cluster_id": ["C1", "C1"],
        }
    ).to_csv(path, index=False)
    return path


@pytest.fixture()
def sample_metadata(tmp_path: Path) -> Path:
    path = tmp_path / "metadata.jsonl"
    lines = [
        json.dumps({"record_id": "SEQ001", "accession": "NC_001", "organism": "Homo sapiens"}),
        json.dumps({"record_id": "SEQ002", "accession": "NC_002", "organism": "Mus musculus"}),
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def test_export_creates_parquet_files(sample_features: Path, sample_metadata: Path, tmp_path: Path) -> None:
    out = tmp_path / "export"
    paths = export_ml_ready(sample_features, sample_metadata, out, validate=True)

    assert paths.features_parquet.exists()
    assert paths.joined_parquet.exists()
    assert paths.manifest_json.exists()

    joined = pd.read_parquet(paths.joined_parquet)
    assert len(joined) == 2
    assert "organism" in joined.columns


def test_manifest_schema_version(sample_features: Path, sample_metadata: Path, tmp_path: Path) -> None:
    paths = export_ml_ready(sample_features, sample_metadata, tmp_path / "out")
    manifest = json.loads(paths.manifest_json.read_text())
    assert manifest["schema_version"] == "1"
    assert manifest["counts"]["joined_rows"] == 2


def test_export_without_validation(sample_features: Path, sample_metadata: Path, tmp_path: Path) -> None:
    """L'export sans validation doit quand même fonctionner."""
    paths = export_ml_ready(sample_features, sample_metadata, tmp_path / "out", validate=False)
    assert paths.joined_parquet.exists()
