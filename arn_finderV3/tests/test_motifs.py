"""Tests du module motifs V3."""

from __future__ import annotations

from pathlib import Path

import pytest

from arn_finder_v3.motifs import MotifRequest, compute_motifs, _extract_kmers, _canonical_form


FASTA = """\
>SEQ1 test
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
>SEQ2 test
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
"""


@pytest.fixture()
def fasta_file(tmp_path: Path) -> Path:
    p = tmp_path / "seq.fasta"
    p.write_text(FASTA, encoding="utf-8")
    return p


def test_compute_motifs_basic(fasta_file: Path, tmp_path: Path) -> None:
    req = MotifRequest(input_fasta=fasta_file, out_dir=tmp_path / "out", k=4, top=50)
    res = compute_motifs(req)
    assert res.total_records == 2
    assert res.total_kmers > 0
    assert res.motifs_csv.exists()


def test_per_record_csv_created(fasta_file: Path, tmp_path: Path) -> None:
    req = MotifRequest(input_fasta=fasta_file, out_dir=tmp_path / "out", k=4, per_record=True)
    res = compute_motifs(req)
    assert res.motifs_by_record_csv is not None
    assert res.motifs_by_record_csv.exists()


def test_extract_kmers_count() -> None:
    _, counts = _extract_kmers("ATGCATGC", 4, False, False, "DNA")
    assert counts["ATGC"] > 0
    assert sum(counts.values()) == 5  # len 8 - k 4 + 1


def test_canonical_form_symmetry() -> None:
    kmer = "ATGCATGCA"
    canon = _canonical_form(kmer, "DNA")
    revcomp = "TGCATGCAT"
    assert canon == min(kmer, revcomp)


def test_ignore_n_kmers() -> None:
    _, counts = _extract_kmers("ATGNATGC", 4, True, False, "DNA")
    for kmer in counts:
        assert "N" not in kmer
