"""Tests du rapport HTML V3.1."""

from __future__ import annotations

import json
from pathlib import Path

import pytest


def _has_jinja2() -> bool:
    try:
        import jinja2  # noqa: F401
        return True
    except ImportError:
        return False


pytestmark = pytest.mark.skipif(not _has_jinja2(), reason="jinja2 non installé")


@pytest.fixture()
def data_dir(tmp_path: Path) -> Path:
    """Arborescence pipeline minimale pour le rapport."""
    (tmp_path / "filtered").mkdir()
    (tmp_path / "filtered" / "filter_manifest.json").write_text(
        json.dumps({
            "total_records": 10, "kept_records": 8,
            "reason_counts": {"ok": 8, "too_short": 2},
        }),
        encoding="utf-8",
    )

    (tmp_path / "motifs").mkdir()
    (tmp_path / "motifs" / "motifs.csv").write_text(
        "kmer,count,freq,n_records,idf_like\n"
        "ATGCATGCA,50,0.1,8,0.2\n"
        "GCATGCATG,30,0.06,6,0.3\n",
        encoding="utf-8",
    )

    (tmp_path / "clusters").mkdir()
    (tmp_path / "clusters" / "cluster_manifest.json").write_text(
        json.dumps({"cluster_count": 3, "params": {"min_similarity": 0.35}}),
        encoding="utf-8",
    )
    (tmp_path / "clusters" / "clusters.jsonl").write_text(
        json.dumps({
            "cluster_id": "C0001", "size": 5,
            "dominant_organism": "Homo sapiens",
            "mean_pairwise_sim": 0.72,
            "members": [{"record_id": "SEQ1", "accession": "SEQ1", "length": 300, "organism": "Homo sapiens"}],
        }) + "\n",
        encoding="utf-8",
    )

    (tmp_path / "structure").mkdir()
    (tmp_path / "structure" / "structure_manifest.json").write_text(
        json.dumps({"computed": 8, "skipped": 0, "backend": "nussinov_approx"}),
        encoding="utf-8",
    )
    (tmp_path / "structure" / "structures.csv").write_text(
        "record_id,length,gc_frac,dot_bracket,mfe_kcal_mol,backend,skipped,skip_reason\n"
        "SEQ1,300,0.45,(((...))),-5.2,nussinov_approx,False,\n",
        encoding="utf-8",
    )

    return tmp_path


def test_report_generated(data_dir: Path, tmp_path: Path) -> None:
    from arn_finder_v3.report import generate_report

    out = tmp_path / "report.html"
    result = generate_report(data_dir=data_dir, query="BRCA1 test", out_path=out)
    assert result.exists()
    assert result == out


def test_report_contains_query(data_dir: Path, tmp_path: Path) -> None:
    from arn_finder_v3.report import generate_report

    out = tmp_path / "report.html"
    generate_report(data_dir=data_dir, query="BRCA1[Gene Name]", out_path=out)
    content = out.read_text(encoding="utf-8")
    assert "BRCA1" in content


def test_report_contains_stats(data_dir: Path, tmp_path: Path) -> None:
    from arn_finder_v3.report import generate_report

    out = tmp_path / "report.html"
    generate_report(data_dir=data_dir, query="test", out_path=out)
    content = out.read_text(encoding="utf-8")
    # Stats filtre
    assert "10" in content  # total records
    assert "8" in content   # kept records
    # Backend structure
    assert "nussinov_approx" in content


def test_report_default_path(data_dir: Path) -> None:
    from arn_finder_v3.report import generate_report

    result = generate_report(data_dir=data_dir, query="test")
    assert result == data_dir / "report.html"
    assert result.exists()
