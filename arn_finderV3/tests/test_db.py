"""Tests DuckDB V3.1."""

from __future__ import annotations

import json
from pathlib import Path

import pytest


def _has_duckdb() -> bool:
    try:
        import duckdb  # noqa: F401
        return True
    except ImportError:
        return False


pytestmark = pytest.mark.skipif(not _has_duckdb(), reason="duckdb non installé")


@pytest.fixture()
def data_dir(tmp_path: Path) -> Path:
    """Crée une arborescence de données pipeline minimale."""
    # filter_report.csv
    filtered = tmp_path / "filtered"
    filtered.mkdir()
    (filtered / "filter_report.csv").write_text(
        "record_id,accession,length,n_count,n_frac,gc_percent,kept,reason,original_header\n"
        "SEQ1,NC_001,300,0,0.0,45.0,true,ok,>SEQ1\n"
        "SEQ2,NC_002,150,5,0.03,42.0,false,too_short,>SEQ2\n",
        encoding="utf-8",
    )
    (filtered / "filter_manifest.json").write_text(
        json.dumps({"total_records": 2, "kept_records": 1, "reason_counts": {"ok": 1, "too_short": 1}}),
        encoding="utf-8",
    )

    # motifs_by_record.csv
    motifs = tmp_path / "motifs"
    motifs.mkdir()
    (motifs / "motifs_by_record.csv").write_text(
        "record_id,kmer,count,freq_in_record\nSEQ1,ATGCATGCA,5,0.12\n",
        encoding="utf-8",
    )

    # structures.csv
    struct = tmp_path / "structure"
    struct.mkdir()
    (struct / "structures.csv").write_text(
        "record_id,length,gc_frac,dot_bracket,mfe_kcal_mol,backend,skipped,skip_reason\n"
        "SEQ1,300,0.45,(((...))),-5.2,nussinov_approx,False,\n",
        encoding="utf-8",
    )

    return tmp_path


def test_import_sequences(data_dir: Path, tmp_path: Path) -> None:
    from arn_finder_v3.db import ARNFinderDB

    with ARNFinderDB(tmp_path / "test.duckdb") as db:
        counts = db.import_all(data_dir)
        assert counts["sequences"] == 2
        df = db.query("SELECT * FROM sequences")
        assert len(df) == 2


def test_import_motifs(data_dir: Path, tmp_path: Path) -> None:
    from arn_finder_v3.db import ARNFinderDB

    with ARNFinderDB(tmp_path / "test.duckdb") as db:
        counts = db.import_all(data_dir)
        assert counts["motifs"] == 1
        df = db.query("SELECT * FROM motifs WHERE kmer = 'ATGCATGCA'")
        assert len(df) == 1


def test_import_structure(data_dir: Path, tmp_path: Path) -> None:
    from arn_finder_v3.db import ARNFinderDB

    with ARNFinderDB(tmp_path / "test.duckdb") as db:
        counts = db.import_all(data_dir)
        assert counts["structure"] == 1
        df = db.query("SELECT mfe FROM structure WHERE record_id = 'SEQ1'")
        assert abs(float(df.iloc[0]["mfe"]) - (-5.2)) < 0.01


def test_query_sql(data_dir: Path, tmp_path: Path) -> None:
    from arn_finder_v3.db import ARNFinderDB

    with ARNFinderDB(tmp_path / "test.duckdb") as db:
        db.import_all(data_dir)
        df = db.query("SELECT record_id FROM sequences WHERE kept = TRUE")
        assert list(df["record_id"]) == ["SEQ1"]


def test_tables_list(data_dir: Path, tmp_path: Path) -> None:
    from arn_finder_v3.db import ARNFinderDB

    with ARNFinderDB(tmp_path / "test.duckdb") as db:
        tables = db.tables()
        assert "sequences" in tables
        assert "motifs" in tables
        assert "structure" in tables


def test_export_parquet(data_dir: Path, tmp_path: Path) -> None:
    from arn_finder_v3.db import ARNFinderDB

    with ARNFinderDB(tmp_path / "test.duckdb") as db:
        db.import_all(data_dir)
        out = tmp_path / "sequences.parquet"
        db.export_parquet("sequences", out)
        assert out.exists()
