"""Tests du module consensus V3."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from arn_finder_v3.consensus import ConsensusRequest, build_consensus
from arn_finder_v3.consensus.merge import merge_sequences, _overlap_identity


FASTA = """\
>SEQ1
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
>SEQ2
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGT
"""

CLUSTERS_JSONL = json.dumps({
    "cluster_id": "C0001",
    "size": 2,
    "members": [
        {"record_id": "SEQ1", "accession": "SEQ1", "length": 80, "organism": "Test"},
        {"record_id": "SEQ2", "accession": "SEQ2", "length": 80, "organism": "Test"},
    ],
    "dominant_organism": "Test",
    "n_unique_organisms": 1,
    "mean_len": 80.0,
    "max_pairwise_sim": 0.9,
    "mean_pairwise_sim": 0.9,
})


@pytest.fixture()
def data_dir(tmp_path: Path) -> Path:
    (tmp_path / "seqs.fasta").write_text(FASTA, encoding="utf-8")
    (tmp_path / "clusters.jsonl").write_text(CLUSTERS_JSONL + "\n", encoding="utf-8")
    return tmp_path


def test_build_consensus_creates_fasta(data_dir: Path, tmp_path: Path) -> None:
    req = ConsensusRequest(
        clusters_path=data_dir / "clusters.jsonl",
        fasta_path=data_dir / "seqs.fasta",
        out_dir=tmp_path / "out",
        min_cluster_size=2,
    )
    manifest = build_consensus(req)
    assert (tmp_path / "out" / "consensus.fasta").exists()
    assert manifest["nb_clusters_attempted"] == 1


def test_min_cluster_size_skips(data_dir: Path, tmp_path: Path) -> None:
    req = ConsensusRequest(
        clusters_path=data_dir / "clusters.jsonl",
        fasta_path=data_dir / "seqs.fasta",
        out_dir=tmp_path / "out",
        min_cluster_size=5,  # > 2 donc sauté
    )
    manifest = build_consensus(req)
    assert manifest["nb_clusters_attempted"] == 0


def test_merge_sequences_overlap() -> None:
    a = "ATGCATGCATGC"
    b = "ATGCATGCTTTT"
    result = merge_sequences(a, b, min_overlap=8, min_identity=0.85)
    assert result is not None
    assert len(result.merged) < len(a) + len(b)


def test_merge_sequences_no_overlap() -> None:
    a = "AAAAAAAAAA"
    b = "TTTTTTTTTT"
    result = merge_sequences(a, b, min_overlap=5, min_identity=0.99)
    assert result is None


def test_overlap_identity_perfect() -> None:
    assert _overlap_identity("ATGC", "ATGC") == 1.0


def test_overlap_identity_with_n() -> None:
    # N ne compte pas dans le dénominateur
    identity = _overlap_identity("ATNG", "ATCG")
    assert identity == 1.0  # 2 matchs (A,T,G) sur 3 positions effectives... wait
    # A==A, T==T, N ignoré, G==G → 3/3
    assert abs(identity - 1.0) < 1e-6
