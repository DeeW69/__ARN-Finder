from pathlib import Path

import json
import pandas as pd
import pytest

pytest.importorskip("pyarrow")

from arn_finder_v2.exporter import export_ml_ready


def test_export_ml_ready_writes_parquet_and_manifest(tmp_path: Path):
    features_path = tmp_path / "features.csv"
    metadata_path = tmp_path / "metadata.jsonl"
    motifs_path = tmp_path / "motifs.csv"
    out_dir = tmp_path / "out"

    features_path.write_text("record_id,value\nrec1,1\nrec2,2\n", encoding="utf-8")
    metadata_entries = [
        {"accession": "rec1", "taxon": "Virus A"},
        {"accession": "rec2", "taxon": "Virus B"},
    ]
    metadata_path.write_text("\n".join(json.dumps(entry) for entry in metadata_entries), encoding="utf-8")
    motifs_path.write_text("record_id,motif\nrec1,M1\n", encoding="utf-8")

    export_paths = export_ml_ready(
        features_csv=features_path,
        metadata_jsonl=metadata_path,
        motifs_csv=motifs_path,
        out_dir=out_dir,
        cli_version="0.1.0",
    )

    assert export_paths.features_parquet.exists()
    assert export_paths.joined_parquet.exists()
    assert export_paths.motifs_parquet.exists()
    manifest = json.loads(export_paths.manifest_json.read_text(encoding="utf-8"))
    assert manifest["counts"]["features_rows"] == 2
    assert manifest["join_keys"] == ["record_id", "accession"]

    joined = pd.read_parquet(export_paths.joined_parquet)
    assert "taxon" in joined.columns
    assert joined.loc[joined["record_id"] == "rec1", "taxon"].iloc[0] == "Virus A"
