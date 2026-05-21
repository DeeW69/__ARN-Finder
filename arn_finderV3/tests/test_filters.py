"""Tests du module filters V3."""

from __future__ import annotations

from pathlib import Path

import pytest

from arn_finder_v3.filters import FilterRequest, filter_sequences, _compute_metrics, _clean_sequence


FASTA_CONTENT = """\
>NC_001 Homo sapiens chromosome 1
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
>NC_002 Mus musculus short
ATGCATGCATGC
>NC_003 organism_n_heavy
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
"""


@pytest.fixture()
def fasta_file(tmp_path: Path) -> Path:
    p = tmp_path / "test.fasta"
    p.write_text(FASTA_CONTENT, encoding="utf-8")
    return p


def test_filter_keeps_valid(fasta_file: Path, tmp_path: Path) -> None:
    req = FilterRequest(input_fasta=fasta_file, out_dir=tmp_path / "out", min_len=200)
    res = filter_sequences(req)
    assert res.kept_records == 1
    assert res.total_records == 3
    assert res.dropped_by_reason["too_short"] == 1
    assert res.fasta_path.exists()


def test_filter_report_created(fasta_file: Path, tmp_path: Path) -> None:
    req = FilterRequest(input_fasta=fasta_file, out_dir=tmp_path / "out", min_len=200, write_report=True)
    res = filter_sequences(req)
    assert res.report_path is not None
    assert res.report_path.exists()
    lines = res.report_path.read_text().splitlines()
    assert len(lines) == 4  # header + 3 records


def test_filter_gc_bounds(tmp_path: Path) -> None:
    """Le filtre GC doit rejeter les séquences hors bornes."""
    # Séquence 100% GC
    fasta = tmp_path / "gc.fasta"
    fasta.write_text(">SEQ1 test\n" + "GC" * 150 + "\n", encoding="utf-8")
    req = FilterRequest(
        input_fasta=fasta, out_dir=tmp_path / "out",
        min_len=200, max_gc=0.6,
    )
    res = filter_sequences(req)
    assert res.kept_records == 0
    assert res.dropped_by_reason.get("gc_too_high", 0) == 1


def test_filter_dedupe(tmp_path: Path) -> None:
    """Les séquences dupliquées doivent être détectées."""
    seq = "ATGC" * 75  # 300 nt
    fasta = tmp_path / "dup.fasta"
    fasta.write_text(f">SEQ1 a\n{seq}\n>SEQ2 b\n{seq}\n", encoding="utf-8")
    req = FilterRequest(input_fasta=fasta, out_dir=tmp_path / "out", min_len=200, dedupe=True)
    res = filter_sequences(req)
    assert res.kept_records == 1
    assert res.dropped_by_reason.get("duplicate", 0) == 1


def test_compute_metrics_basic() -> None:
    seq = "ATGCNNNN"
    m = _compute_metrics(seq)
    assert m["length"] == 8
    assert m["n_count"] == 4
    assert abs(m["n_frac"] - 0.5) < 1e-6


def test_clean_sequence_replaces_ambiguous() -> None:
    cleaned = _clean_sequence("ATGCRYWWSS")
    assert "R" not in cleaned
    assert "N" in cleaned
