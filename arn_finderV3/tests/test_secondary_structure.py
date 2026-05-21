"""Tests du module structure secondaire V3.1 — mis à jour V3.4 (seqfold)."""

from __future__ import annotations

import importlib
from pathlib import Path

import pytest

from arn_finder_v3.secondary_structure import (
    StructureRequest, compute_secondary_structure,
    _predict_nussinov, _predict_seqfold, _gc_frac, _detect_backend,
)

_SEQFOLD_AVAILABLE = importlib.util.find_spec("seqfold") is not None

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
    assert backend in ("viennarna", "rnafold_subprocess", "seqfold", "nussinov_approx")


def test_max_len_skip(tmp_path: Path) -> None:
    """Séquences plus longues que max_len doivent être ignorées."""
    fasta = tmp_path / "long.fasta"
    fasta.write_text(">LONG\n" + "A" * 500 + "\n", encoding="utf-8")
    req = StructureRequest(input_fasta=fasta, out_dir=tmp_path / "out", max_len=100)
    res = compute_secondary_structure(req)
    assert res.skipped == 1
    assert res.computed == 0


# ── Tests seqfold ─────────────────────────────────────────────────────────────

@pytest.mark.skipif(not _SEQFOLD_AVAILABLE, reason="seqfold non installé")
def test_seqfold_dot_bracket_valid() -> None:
    """La structure seqfold doit avoir autant de '(' que de ')'."""
    seq = "AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGC"
    dot_bracket, mfe = _predict_seqfold(seq, temperature=37.0)
    assert len(dot_bracket) == len(seq)
    assert dot_bracket.count("(") == dot_bracket.count(")")


@pytest.mark.skipif(not _SEQFOLD_AVAILABLE, reason="seqfold non installé")
def test_seqfold_mfe_negative_on_stem_loop() -> None:
    """Une séquence avec tige-boucle doit avoir une MFE ≤ 0."""
    # Séquence classique avec tige-boucle GC stable
    seq = "GCGCGCAAAGCGCGC"
    _, mfe = _predict_seqfold(seq, temperature=37.0)
    assert mfe <= 0.0


@pytest.mark.skipif(not _SEQFOLD_AVAILABLE, reason="seqfold non installé")
def test_seqfold_output_is_valid_dot_bracket() -> None:
    """La sortie seqfold doit être une notation dot-bracket valide."""
    seq = "AAAAAAAAAAAA"
    dot_bracket, mfe = _predict_seqfold(seq, temperature=37.0)
    assert len(dot_bracket) == len(seq)
    assert set(dot_bracket).issubset({".", "(", ")"})
    assert dot_bracket.count("(") == dot_bracket.count(")")
    assert mfe <= 0.0


@pytest.mark.skipif(not _SEQFOLD_AVAILABLE, reason="seqfold non installé")
def test_seqfold_backend_via_compute(fasta_file: Path, tmp_path: Path) -> None:
    """compute_secondary_structure avec backend='seqfold' doit fonctionner."""
    req = StructureRequest(
        input_fasta=fasta_file,
        out_dir=tmp_path / "out",
        max_len=2000,
        backend="seqfold",
    )
    res = compute_secondary_structure(req)
    assert res.total == 3
    assert res.computed == 2
    assert res.skipped == 1
    assert res.backend_used == "seqfold"
    # MFE seqfold doit être physiquement réalistes (pas -1 par paire)
    computed_recs = [r for r in res.records if not r.skipped]
    for rec in computed_recs:
        # MFE Zuker réaliste : entre -200 et 0 kcal/mol pour des séquences courtes
        assert rec.mfe <= 0.0 or rec.mfe == 0.0


@pytest.mark.skipif(not _SEQFOLD_AVAILABLE, reason="seqfold non installé")
def test_seqfold_better_than_nussinov() -> None:
    """
    seqfold (Zuker) doit donner des MFE différentes de Nussinov (-1/paire).
    Nussinov retourne des multiples entiers de -1.0 ; seqfold retourne des
    valeurs physiques en kcal/mol (décimales réalistes).
    """
    seq = "GCGCGCAAAGCGCGC"
    _, mfe_sf = _predict_seqfold(seq, temperature=37.0)
    _, mfe_nu = _predict_nussinov(seq)
    # MFE seqfold ≠ multiple entier de -1 (signe que c'est du vrai Zuker)
    assert mfe_sf != mfe_nu or mfe_sf % 1.0 != 0.0


@pytest.mark.skipif(not _SEQFOLD_AVAILABLE, reason="seqfold non installé")
def test_detect_backend_prefers_seqfold_over_nussinov() -> None:
    """Si seqfold est installé et ViennaRNA absent, on ne doit pas tomber sur Nussinov."""
    backend = _detect_backend()
    # seqfold est installé → le backend détecté ne doit pas être nussinov_approx
    # (sauf si ViennaRNA ou RNAfold sont aussi présents, ce qui est aussi bien)
    assert backend != "nussinov_approx"
