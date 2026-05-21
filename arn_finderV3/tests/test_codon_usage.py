"""Tests du module codon_usage V3.2."""

from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture()
def fasta_file(tmp_path: Path) -> Path:
    """FASTA minimal avec séquences codantes connues."""
    p = tmp_path / "seqs.fasta"
    # Séquence 1 : ATG (Met) + AAA AAA (Lys×2) + TAA (stop) → 4 codons, biaisée Lys
    # Séquence 2 : 30 codons variés pour un ENc raisonnable
    p.write_text(
        ">SEQ1 test sequence 1\n"
        "ATGAAAAAAAAAAAAAAAAAAAAAAAAAATAA\n"
        ">SEQ2 test sequence 2\n"
        "ATGAAAGCCTTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATAA\n",
        encoding="utf-8",
    )
    return p


def test_compute_codon_usage_runs(fasta_file: Path, tmp_path: Path) -> None:
    from arn_finder_v3.codon_usage import CodonUsageRequest, compute_codon_usage

    req = CodonUsageRequest(input_fasta=fasta_file, out_dir=tmp_path / "cu", min_codons=3)
    result = compute_codon_usage(req)
    assert result.computed >= 1
    assert result.stats_csv.exists()
    assert result.rscu_csv.exists()
    assert result.manifest_path.exists()


def test_stats_csv_has_expected_columns(fasta_file: Path, tmp_path: Path) -> None:
    import csv
    from arn_finder_v3.codon_usage import CodonUsageRequest, compute_codon_usage

    req = CodonUsageRequest(input_fasta=fasta_file, out_dir=tmp_path / "cu", min_codons=3)
    result = compute_codon_usage(req)

    with result.stats_csv.open(encoding="utf-8") as f:
        headers = next(csv.reader(f))
    for col in ("record_id", "gc1", "gc2", "gc3", "enc", "total_codons"):
        assert col in headers, f"Colonne manquante : {col}"


def test_rscu_csv_has_expected_columns(fasta_file: Path, tmp_path: Path) -> None:
    import csv
    from arn_finder_v3.codon_usage import CodonUsageRequest, compute_codon_usage

    req = CodonUsageRequest(input_fasta=fasta_file, out_dir=tmp_path / "cu", min_codons=3)
    result = compute_codon_usage(req)

    with result.rscu_csv.open(encoding="utf-8") as f:
        headers = next(csv.reader(f))
    for col in ("record_id", "codon", "amino_acid", "rscu"):
        assert col in headers, f"Colonne manquante : {col}"


def test_enc_in_valid_range(fasta_file: Path, tmp_path: Path) -> None:
    import csv
    from arn_finder_v3.codon_usage import CodonUsageRequest, compute_codon_usage

    req = CodonUsageRequest(input_fasta=fasta_file, out_dir=tmp_path / "cu", min_codons=3)
    result = compute_codon_usage(req)

    with result.stats_csv.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            enc = float(row["enc"])
            assert 20.0 <= enc <= 61.0, f"ENc hors range [20,61] : {enc}"


def test_rscu_neutral_for_single_codon_aa() -> None:
    """Met (AUG seul) doit avoir RSCU = 1.0 quelle que soit la séquence."""
    from arn_finder_v3.codon_usage import AA_CODON_GROUPS, _rscu

    # Séquence avec plusieurs Met
    counts = {c: 0 for c in ["AUG", "UAA", "UUU", "UUC"]}
    counts["AUG"] = 5
    counts["UUU"] = 3
    counts["UUC"] = 7
    counts["UAA"] = 1

    rscu = _rscu(counts)
    assert rscu["AUG"] == 1.0, f"RSCU(AUG) attendu 1.0, obtenu {rscu['AUG']}"


def test_rscu_biased_codon() -> None:
    """Un codon très utilisé parmi ses synonymes doit avoir RSCU > 1."""
    from arn_finder_v3.codon_usage import _rscu

    counts = {c: 0 for c in ["UUU", "UUC"]}
    counts["UUU"] = 100   # très préféré
    counts["UUC"] = 0

    rscu = _rscu(counts)
    assert rscu["UUU"] == 2.0, f"RSCU(UUU) attendu 2.0, obtenu {rscu['UUU']}"
    assert rscu["UUC"] == 0.0


def test_gc_by_position() -> None:
    from arn_finder_v3.codon_usage import _gc_by_position

    # GCG : G(pos1)=GC, C(pos2)=GC, G(pos3)=GC → gc1=1, gc2=1, gc3=1
    gc1, gc2, gc3 = _gc_by_position(["GCG", "GCG"])
    assert gc1 == 1.0
    assert gc2 == 1.0
    assert gc3 == 1.0

    # AUU : A=non-GC, U=non-GC, U=non-GC → tous à 0
    gc1, gc2, gc3 = _gc_by_position(["AUU"])
    assert gc1 == 0.0
    assert gc2 == 0.0
    assert gc3 == 0.0


def test_dna_converted_to_rna(tmp_path: Path) -> None:
    """Les séquences ADN (T) doivent être converties en ARN (U) avant calcul."""
    from arn_finder_v3.codon_usage import CodonUsageRequest, compute_codon_usage
    import csv

    fasta = tmp_path / "dna.fasta"
    # ATG (ADN) = AUG (ARN) → Met
    fasta.write_text(">DNA_SEQ\nATGAAAAAAATAA\n", encoding="utf-8")

    req = CodonUsageRequest(input_fasta=fasta, out_dir=tmp_path / "cu", min_codons=2)
    result = compute_codon_usage(req)
    assert result.computed == 1

    with result.rscu_csv.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    # AUG doit apparaître dans les résultats
    aug_row = next((r for r in rows if r["codon"] == "AUG"), None)
    assert aug_row is not None, "AUG absent du CSV RSCU"


def test_min_codons_filter(tmp_path: Path) -> None:
    """Séquences avec moins de min_codons codons doivent être ignorées."""
    from arn_finder_v3.codon_usage import CodonUsageRequest, compute_codon_usage

    fasta = tmp_path / "short.fasta"
    fasta.write_text(">SHORT\nATGTAA\n", encoding="utf-8")  # 2 codons seulement

    req = CodonUsageRequest(input_fasta=fasta, out_dir=tmp_path / "cu", min_codons=5)
    result = compute_codon_usage(req)
    assert result.computed == 0
    assert result.skipped == 1
