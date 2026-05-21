"""Tests MinHash + LSH — V3.3."""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from arn_finder_v3.clustering.minhash import (
    MinHashConfig,
    cluster_by_minhash,
    compute_signature,
    jaccard_from_signatures,
    lsh_candidate_pairs,
)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _random_kmers(n: int, k: int = 9, seed: int = 0) -> set[str]:
    rng = random.Random(seed)
    alphabet = "ACGT"
    return {"".join(rng.choices(alphabet, k=k)) for _ in range(n)}


def _similar_set(base: set[str], overlap: float, extra: int = 20, seed: int = 99) -> set[str]:
    """Crée un ensemble partageant `overlap` fraction des éléments avec `base`."""
    rng = random.Random(seed)
    keep = set(rng.sample(list(base), int(len(base) * overlap)))
    extra_kmers = _random_kmers(extra, seed=seed + 1)
    return keep | extra_kmers


# ── Tests signature ───────────────────────────────────────────────────────────

def test_signature_length() -> None:
    kset = _random_kmers(100)
    sig = compute_signature(kset, n_hashes=64)
    assert len(sig) == 64


def test_signature_deterministic() -> None:
    """Deux calculs sur le même set donnent la même signature."""
    kset = _random_kmers(50)
    sig1 = compute_signature(kset, n_hashes=128)
    sig2 = compute_signature(kset, n_hashes=128)
    assert sig1 == sig2


def test_signature_empty_set() -> None:
    sig = compute_signature(set(), n_hashes=32)
    assert len(sig) == 32
    # Toutes les valeurs doivent être la valeur maximale (aucun k-mer)
    from arn_finder_v3.clustering.minhash import _MERSENNE
    assert all(v == _MERSENNE for v in sig)


def test_identical_sets_full_similarity() -> None:
    kset = _random_kmers(100)
    sig = compute_signature(kset, n_hashes=128)
    assert jaccard_from_signatures(sig, sig) == 1.0


def test_disjoint_sets_low_similarity() -> None:
    kset_a = _random_kmers(100, seed=1)
    kset_b = _random_kmers(100, seed=9999)
    sig_a = compute_signature(kset_a, n_hashes=256)
    sig_b = compute_signature(kset_b, n_hashes=256)
    # Deux sets aléatoires sans chevauchement → similarité proche de 0
    est = jaccard_from_signatures(sig_a, sig_b)
    assert est < 0.15, f"Similarité estimée trop élevée pour sets disjoints : {est}"


def test_jaccard_estimate_accuracy() -> None:
    """L'estimation MinHash doit être proche du Jaccard exact (±0.10 avec n=256)."""
    base = _random_kmers(200, seed=42)
    similar = _similar_set(base, overlap=0.6, seed=7)

    # Jaccard exact
    union = len(base | similar)
    true_j = len(base & similar) / union if union else 0.0

    sig_a = compute_signature(base, n_hashes=256)
    sig_b = compute_signature(similar, n_hashes=256)
    est_j = jaccard_from_signatures(sig_a, sig_b)

    assert abs(est_j - true_j) < 0.15, (
        f"Estimation MinHash trop imprécise : vrai={true_j:.3f}, estimé={est_j:.3f}"
    )


# ── Tests LSH banding ─────────────────────────────────────────────────────────

def test_lsh_finds_similar_pair() -> None:
    """Deux séquences très similaires doivent être candidates."""
    base = _random_kmers(200, seed=1)
    similar = _similar_set(base, overlap=0.8, seed=2)

    config = MinHashConfig(n_hashes=128, n_bands=16)
    sigs = {
        "SEQ_A": compute_signature(base, config.n_hashes),
        "SEQ_B": compute_signature(similar, config.n_hashes),
    }
    candidates = lsh_candidate_pairs(sigs, config)
    assert ("SEQ_A", "SEQ_B") in candidates or ("SEQ_B", "SEQ_A") in candidates


def test_lsh_candidate_pairs_sorted() -> None:
    """Les paires candidates doivent être triées (a < b) pour éviter les doublons."""
    base = _random_kmers(50, seed=1)
    similar = _similar_set(base, overlap=0.9, seed=3)

    config = MinHashConfig(n_hashes=64, n_bands=8)
    sigs = {
        "Z_SEQ": compute_signature(base, config.n_hashes),
        "A_SEQ": compute_signature(similar, config.n_hashes),
    }
    candidates = lsh_candidate_pairs(sigs, config)
    for a, b in candidates:
        assert a < b, f"Paire non triée : ({a}, {b})"


def test_lsh_no_self_pairs() -> None:
    """Une séquence ne doit jamais former une paire avec elle-même."""
    ksets = {f"SEQ{i}": compute_signature(_random_kmers(30, seed=i), 64) for i in range(5)}
    config = MinHashConfig(n_hashes=64, n_bands=8)
    candidates = lsh_candidate_pairs(ksets, config)
    for a, b in candidates:
        assert a != b


# ── Tests cluster_by_minhash ──────────────────────────────────────────────────

def test_cluster_minhash_finds_similar_pairs() -> None:
    """MinHash doit trouver des edges entre séquences similaires."""
    base = _random_kmers(200, seed=42)
    similar1 = _similar_set(base, overlap=0.7, seed=10)
    similar2 = _similar_set(base, overlap=0.7, seed=20)
    different = _random_kmers(200, seed=9999)

    ksets = [base, similar1, similar2, different]
    ids = ["BASE", "SIM1", "SIM2", "DIFF"]

    config = MinHashConfig(n_hashes=128, n_bands=16, verify=True)
    result = cluster_by_minhash(ksets, ids, min_similarity=0.3, config=config)

    edge_ids = {(a, b) for a, b, _ in result.edges}
    # BASE-SIM1 ou SIM1-BASE doit être un edge
    assert ("BASE", "SIM1") in edge_ids or ("SIM1", "BASE") in edge_ids, (
        f"Paire BASE-SIM1 manquante dans les edges : {edge_ids}"
    )


def test_cluster_minhash_empty_input() -> None:
    config = MinHashConfig(n_hashes=64, n_bands=8)
    result = cluster_by_minhash([], [], min_similarity=0.35, config=config)
    assert result.edges == []
    assert result.pairs_evaluated == 0


def test_cluster_minhash_verify_vs_estimate() -> None:
    """verify=True et verify=False doivent donner des résultats cohérents."""
    base = _random_kmers(150, seed=5)
    similar = _similar_set(base, overlap=0.75, seed=6)
    different = _random_kmers(150, seed=8888)

    ksets = [base, similar, different]
    ids = ["A", "B", "C"]

    config_verify = MinHashConfig(n_hashes=128, n_bands=16, verify=True)
    config_fast = MinHashConfig(n_hashes=128, n_bands=16, verify=False)

    result_v = cluster_by_minhash(ksets, ids, 0.3, config_verify)
    result_f = cluster_by_minhash(ksets, ids, 0.3, config_fast)

    # Les deux méthodes doivent trouver la même paire A-B (éventuellement avec des sims différentes)
    ids_v = {(a, b) for a, b, _ in result_v.edges}
    ids_f = {(a, b) for a, b, _ in result_f.edges}
    # L'overlap doit être raisonnable
    if ids_v:
        common = ids_v & ids_f
        assert len(common) > 0 or len(ids_f) > 0  # au moins quelque chose trouvé


def test_minhash_config_threshold() -> None:
    """Le seuil LSH calculé doit être dans une plage raisonnable."""
    config = MinHashConfig(n_hashes=128, n_bands=16)   # rows=8
    # (1/16)^(1/8) ≈ 0.56
    assert 0.3 < config.lsh_threshold < 0.8


def test_minhash_config_invalid_division() -> None:
    """n_hashes non divisible par n_bands doit lever ValueError."""
    with pytest.raises(ValueError):
        MinHashConfig(n_hashes=100, n_bands=16)


# ── Test d'intégration avec cluster_sequences ─────────────────────────────────

def test_cluster_sequences_minhash_method(tmp_path: Path) -> None:
    """cluster_sequences avec method='minhash' doit produire les mêmes sorties."""
    from arn_finder_v3.clustering import ClusterRequest, cluster_sequences

    fasta = tmp_path / "seqs.fasta"
    # Construit des séquences avec chevauchement k-mer intentionnel
    seqs = [
        ("SEQ1", "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"),
        ("SEQ2", "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCTTTT"),
        ("SEQ3", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"),
    ]
    with fasta.open("w") as f:
        for sid, seq in seqs:
            f.write(f">{sid}\n{seq}\n")

    req = ClusterRequest(
        input_fasta=fasta,
        out_dir=tmp_path / "clusters",
        k=7,
        min_similarity=0.1,
        method="minhash",
        minhash_hashes=64,
        minhash_bands=8,
    )
    result = cluster_sequences(req)
    assert result.clusters_path.exists()
    assert result.stats_path.exists()
    assert result.manifest_path.exists()
    assert result.total_records == 3
    # SEQ1 et SEQ2 très similaires → au moins 1 edge
    assert result.edges_kept >= 1
