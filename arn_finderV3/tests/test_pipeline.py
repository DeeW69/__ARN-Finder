"""Tests d'intégration du pipeline complet V3 (sans réseau)."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from arn_finder_v3.pipeline import PipelineConfig, run_pipeline


FASTA_CONTENT = """\
>NC_001 Homo sapiens
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
>NC_002 Homo sapiens
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
>NC_003 Mus musculus
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
"""

METADATA_CONTENT = """\
{"accession": "NC_001", "organism": "Homo sapiens", "uid": 1}
{"accession": "NC_002", "organism": "Homo sapiens", "uid": 2}
{"accession": "NC_003", "organism": "Mus musculus", "uid": 3}
"""


@pytest.fixture()
def mock_fetch_result(tmp_path: Path):
    """Simule le résultat d'un fetch NCBI sans réseau."""
    raw_dir = tmp_path / "data" / "raw"
    raw_dir.mkdir(parents=True)
    fasta = raw_dir / "sequences.fasta"
    fasta.write_text(FASTA_CONTENT, encoding="utf-8")
    meta = raw_dir / "metadata.jsonl"
    meta.write_text(METADATA_CONTENT, encoding="utf-8")

    from arn_finder_v3.fetch import FetchResult
    return FetchResult(fasta_path=fasta, fetched=3, elapsed=0.1, ids=["NC_001", "NC_002", "NC_003"])


def test_pipeline_without_fetch(mock_fetch_result, tmp_path: Path) -> None:
    """Teste filter → motifs → cluster → consensus → export sans fetch réseau."""
    from arn_finder_v3.fetch import FetchResult

    with patch("arn_finder_v3.pipeline.fetch_sequences", return_value=mock_fetch_result):
        cfg = PipelineConfig(
            query="test query",
            out_base=tmp_path / "data",
            run_fetch=True,
            run_filter=True,
            run_motifs=True,
            run_cluster=True,
            run_consensus=True,
            run_blast=False,
            run_export=False,   # export requiert features CSV valide
            filter_min_len=200,
        )
        result = run_pipeline(cfg)

    # filter et motifs doivent avoir tourné
    assert "filter" in result.stage_results
    assert "motifs" in result.stage_results
    assert "cluster" in result.stage_results
    assert result.stage_results["filter"].kept_records > 0


def test_pipeline_stops_on_error(tmp_path: Path) -> None:
    """stop_on_error=True doit interrompre après la première erreur."""
    from arn_finder_v3.fetch import FetchResult

    bad_result = FetchResult(fasta_path=tmp_path / "nonexistent.fasta", fetched=0, elapsed=0.0)

    with patch("arn_finder_v3.pipeline.fetch_sequences", return_value=bad_result):
        cfg = PipelineConfig(
            query="test",
            out_base=tmp_path / "data",
            run_fetch=True,
            run_filter=True,
            run_motifs=True,
            run_blast=False,
            run_export=False,
            stop_on_error=True,
        )
        result = run_pipeline(cfg)

    assert not result.success
    # motifs ne devrait pas avoir tourné (erreur en filter)
    assert "motifs" not in result.stage_results
