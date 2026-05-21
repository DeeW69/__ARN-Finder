"""Tests du module clustering V3."""

from __future__ import annotations

from pathlib import Path

import pytest

from arn_finder_v3.clustering import ClusterRequest, cluster_sequences, _jaccard
from arn_finder_v3.clustering.union_find import UnionFind


# Deux séquences similaires + une séquence très différente
FASTA = """\
>SEQ1
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
>SEQ2
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGT
>SEQ3
GGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGG
GGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGG
"""


@pytest.fixture()
def fasta_file(tmp_path: Path) -> Path:
    p = tmp_path / "seqs.fasta"
    p.write_text(FASTA, encoding="utf-8")
    return p


def test_cluster_produces_outputs(fasta_file: Path, tmp_path: Path) -> None:
    req = ClusterRequest(input_fasta=fasta_file, out_dir=tmp_path / "out", k=4, min_similarity=0.3)
    res = cluster_sequences(req)
    assert res.clusters_path.exists()
    assert res.edges_path.exists()
    assert res.stats_path.exists()
    assert res.total_records == 3


def test_similar_sequences_grouped(fasta_file: Path, tmp_path: Path) -> None:
    """SEQ1 et SEQ2 sont très similaires, ils devraient être dans le même cluster."""
    req = ClusterRequest(input_fasta=fasta_file, out_dir=tmp_path / "out", k=4, min_similarity=0.5)
    res = cluster_sequences(req)
    assert res.edges_kept >= 1  # au moins une edge entre SEQ1 et SEQ2


def test_jaccard_identical() -> None:
    s = {"ATGC", "TGCA", "GCAT"}
    assert _jaccard(s, s) == 1.0


def test_jaccard_disjoint() -> None:
    a = {"AAAA", "TTTT"}
    b = {"CCCC", "GGGG"}
    assert _jaccard(a, b) == 0.0


def test_union_find_basic() -> None:
    uf = UnionFind(5)
    uf.union(0, 1)
    uf.union(1, 2)
    assert uf.find(0) == uf.find(2)
    assert uf.find(3) != uf.find(0)


def test_union_find_path_compression() -> None:
    uf = UnionFind(10)
    for i in range(9):
        uf.union(i, i + 1)
    root = uf.find(0)
    for i in range(10):
        assert uf.find(i) == root
