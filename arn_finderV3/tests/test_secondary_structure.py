"""Tests du module structure secondaire V3.1."""

from __future__ import annotations

from pathlib import Path

import pytest

from arn_finder_v3.secondary_structure import (
    StructureRequest, compute_secondary_structure,
    _predict_nussinov, _gc_frac, _detect_backend,
)

# FASTA court (< 2000 nt) → doit être traité
SHORT_FASTA = """\
>SEQ1 test
AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUG
>SEQ2 test_dna
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT
>SEQ3 long_skip
""" + "A" * 2100 + "\n"


@pytest.fixture()
def fasta_file(tmp_path: Path) -> Path:
    p = tmp_path / "seqs.fasta"
    p.write_text(SHORT_FASTA, encoding="utf-8")
    return p


def test_structure_nussinov_backend(fasta_file: Path, tmp_path: Path) -> None:
    """Force le backend Nussinov (toujours disponible)."""
    req = StructureRequest(
        input_fasta=fasta_file,
        out_dir=tmp_path / "out",
        max_len=2000,
        backend="nussinov_approx",
    )
    res = compute_secondary_structure(req)
    assert res.total == 3
    assert res.computed == 2     # SEQ1 et SEQ2
    assert res.skipped == 1     # SEQ3 trop long
    assert res.structures_csv.exists()
    assert res.structures_jsonl.exists()
    assert res.manifest_path.exists()


def test_structure_dna_converted_to_rna(fasta_file: Path, tmp_path: Path) -> None:
    """Les séquences DNA doivent être converties en RNA (T→U)."""
    req = StructureRequest(
        input_fasta=fasta_file,
        out_dir=tmp_path / "out",
        backend="nussinov_approx",
    )
    res = compute_secondary_structure(req)
    for rec in res.records:
        if not rec.skipped:
            assert "T" not in rec.sequence_rna


def test_nussinov_dot_bracket_valid() -> None:
    """La structure Nussinov doit avoir autant de '(' que de ')'."""
    seq = "AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGC"
    dot_bracket, mfe = _predict_nussinov(seq)
    assert len(dot_bracket) == len(seq)
    assert dot_bracket.count("(") == dot_bracket.count(")")


def test_nussinov_mfe_negative_on_paired() -> None:
    """Si des paires sont formées, MFE doit être ≤ 0."""
    seq = "AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGC"
    _, mfe = _predict_nussinov(seq)
    assert mfe <= 0.0


def test_nussinov_empty_sequence() -> None:
    """Séquence vide → structure vide, mfe=0."""
    dot_bracket, mfe = _predict_nussinov("")
    assert dot_bracket == ""
    assert mfe == 0.0


def test_gc_frac() -> None:
    assert abs(_gc_frac("GCGC") - 1.0) < 1e-9
    assert abs(_gc_frac("AUAU") - 0.0) < 1e-9
    assert abs(_gc_frac("AUGC") - 0.5) < 1e-9
    assert _gc_frac("") == 0.0


def test_detect_backend_returns_valid() -> None:
    backend = _detect_backend()
    assert backend in ("viennarna", "rnafold_subprocess", "nussinov_approx")


def test_max_len_skip(tmp_path: Path) -> None:
    """Séquences plus longues que max_len doivent être ignorées."""
    fasta = tmp_path / "long.fasta"
    fasta.write_text(">LONG\n" + "A" * 500 + "\n", encoding="utf-8")
    req = StructureRequest(input_fasta=fasta, out_dir=tmp_path / "out", max_len=100)
    res = compute_secondary_structure(req)
    assert res.skipped == 1
    assert res.computed == 0
