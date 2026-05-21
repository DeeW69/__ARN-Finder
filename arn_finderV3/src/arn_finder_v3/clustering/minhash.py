"""
MinHash + LSH — clustering large-scale V3.3

Remplace le Jaccard exact O(n²) par une approche approximative :

  1. Signature MinHash — résumé probabiliste d'un ensemble de k-mers
     Utilise le double-hashing (h1 + i*h2 mod prime) pour dériver n_hashes
     valeurs indépendantes à partir d'un seul blake2b par k-mer. C'est O(|S| × n_hashes)
     avec des opérations arithmétiques simples au lieu de n_hashes calculs de hash.

  2. LSH banding — localité sensible au seuil de similarité
     Découpe la signature en n_bands bandes de rows_per_band lignes.
     Deux séquences sont candidates si elles sont identiques sur AU MOINS une bande
     complète. Probabilité P(candidat) ≈ 1 - (1 - J^r)^b
     Avec n_hashes=128, n_bands=16, r=8 : seuil ~0.45 (P=0.5)

  3. Vérification optionnelle — Jaccard exact sur les paires candidates
     Par défaut activée (verify=True) pour une précision maximale.
     Désactivable (verify=False) pour une performance maximale sur très grands datasets.

Performance attendue :
  - 1k  séquences : jaccard ~0.1s   | minhash ~0.05s
  - 10k séquences : jaccard >60s    | minhash ~2s
  - 100k séquences : jaccard >>     | minhash ~30s

Aucune dépendance externe — 100% stdlib Python 3.11+.

Référence : Broder (1997), Indyk & Motwani (1998)
"""

from __future__ import annotations

import hashlib
import struct
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Iterator

# ── Constantes ────────────────────────────────────────────────────────────────

_MERSENNE: int = (1 << 61) - 1   # Premier de Mersenne M61 — borne pour les hashes
_MAX_HASH: int = _MERSENNE        # Valeur initiale (∞) des signatures


# ── Paramètres MinHash ────────────────────────────────────────────────────────

@dataclass
class MinHashConfig:
    """Paramètres de la signature MinHash et du banding LSH."""
    n_hashes: int = 128        # Taille de la signature (plus = plus précis, plus lent)
    n_bands: int = 16          # Nombre de bandes LSH
    verify: bool = True        # Vérification Jaccard exact après LSH

    @property
    def rows_per_band(self) -> int:
        return self.n_hashes // self.n_bands

    @property
    def lsh_threshold(self) -> float:
        """Seuil de similarité approximatif où P(détection) ≈ 0.5."""
        return (1.0 / self.n_bands) ** (1.0 / self.rows_per_band)

    def __post_init__(self) -> None:
        if self.n_hashes % self.n_bands != 0:
            raise ValueError(
                f"n_hashes ({self.n_hashes}) doit être divisible par n_bands ({self.n_bands})"
            )


# ── Double-hashing MinHash ────────────────────────────────────────────────────

def _double_hash(kmer: str) -> tuple[int, int]:
    """
    Calcule deux valeurs de hash indépendantes pour un k-mer via blake2b 128 bits.
    Déterministe — même kmer → même (h1, h2) quelle que soit la session Python.
    """
    digest = hashlib.blake2b(kmer.encode(), digest_size=16).digest()
    h1, h2 = struct.unpack("<QQ", digest)
    return h1, h2


def compute_signature(kmer_set: set[str], n_hashes: int) -> list[int]:
    """
    Calcule la signature MinHash d'un ensemble de k-mers.

    Algorithme double-hashing (Kirsch & Mitzenmacher 2008) :
      h_i(x) = (h1(x) + i * h2(x)) mod M61

    Complexité : O(|kmer_set| × n_hashes) — opérations arithmétiques rapides.
    """
    sig = [_MAX_HASH] * n_hashes
    for kmer in kmer_set:
        h1, h2 = _double_hash(kmer)
        for i in range(n_hashes):
            h = (h1 + i * h2) % _MERSENNE
            if h < sig[i]:
                sig[i] = h
    return sig


def jaccard_from_signatures(sig_a: list[int], sig_b: list[int]) -> float:
    """
    Estimation de la similarité de Jaccard depuis deux signatures MinHash.
    J_hat = |{i: sig_a[i] == sig_b[i]}| / n_hashes
    Variance attendue : ≈ J(1-J) / n_hashes
    """
    matches = sum(1 for a, b in zip(sig_a, sig_b) if a == b)
    return matches / len(sig_a)


# ── LSH Banding ───────────────────────────────────────────────────────────────

def lsh_candidate_pairs(
    signatures: dict[str, list[int]],
    config: MinHashConfig,
) -> set[tuple[str, str]]:
    """
    Trouve les paires candidates via LSH banding.

    Pour chaque bande b :
      - Hache les rows_per_band valeurs de signature en une clé de bucket
      - Toute paire dans le même bucket sur n'importe quelle bande est candidate

    Complexité : O(n × n_bands) — linéaire en nombre de séquences.

    Retourne un ensemble de paires (id_a, id_b) triées (id_a < id_b).
    """
    r = config.rows_per_band
    candidates: set[tuple[str, str]] = set()

    for b in range(config.n_bands):
        start = b * r
        end = start + r
        buckets: dict[tuple, list[str]] = defaultdict(list)

        for seq_id, sig in signatures.items():
            band_key = tuple(sig[start:end])
            buckets[band_key].append(seq_id)

        for bucket in buckets.values():
            if len(bucket) < 2:
                continue
            # Toutes les paires dans ce bucket
            for i in range(len(bucket)):
                for j in range(i + 1, len(bucket)):
                    a, b_id = bucket[i], bucket[j]
                    if a > b_id:
                        a, b_id = b_id, a
                    candidates.add((a, b_id))

    return candidates


# ── Pipeline MinHash complet ──────────────────────────────────────────────────

@dataclass
class MinHashResult:
    """Résultat du clustering par MinHash — même structure que le clustering Jaccard."""
    edges: list[tuple[str, str, float]]   # (id_a, id_b, sim)
    pairs_evaluated: int
    candidates_found: int
    lsh_threshold: float


def cluster_by_minhash(
    kmer_sets: list[set[str]],
    record_ids: list[str],
    min_similarity: float,
    config: MinHashConfig,
) -> MinHashResult:
    """
    Lance le clustering MinHash + LSH sur une liste de k-mer sets.

    Étapes :
      1. Calcul des signatures MinHash pour chaque séquence
      2. LSH banding pour trouver les paires candidates
      3. Calcul du Jaccard (exact si verify=True, estimé sinon)
      4. Filtre les edges au seuil min_similarity

    Retourne les edges retenus avec leurs similarités.
    """
    n = len(record_ids)
    if n == 0:
        return MinHashResult([], 0, 0, config.lsh_threshold)

    # Étape 1 — Signatures
    id_to_idx = {rid: i for i, rid in enumerate(record_ids)}
    signatures: dict[str, list[int]] = {}
    for i, rid in enumerate(record_ids):
        signatures[rid] = compute_signature(kmer_sets[i], config.n_hashes)

    # Étape 2 — LSH candidates
    candidates = lsh_candidate_pairs(signatures, config)
    candidates_found = len(candidates)

    # Étape 3+4 — Scoring et filtrage
    edges: list[tuple[str, str, float]] = []
    pairs_evaluated = 0

    for id_a, id_b in candidates:
        pairs_evaluated += 1
        if config.verify:
            # Jaccard exact sur les k-mer sets — plus précis
            idx_a = id_to_idx[id_a]
            idx_b = id_to_idx[id_b]
            s_a = kmer_sets[idx_a]
            s_b = kmer_sets[idx_b]
            union = len(s_a | s_b)
            sim = len(s_a & s_b) / union if union else 0.0
        else:
            # Estimation depuis la signature — plus rapide
            sim = jaccard_from_signatures(signatures[id_a], signatures[id_b])

        if sim >= min_similarity:
            edges.append((id_a, id_b, sim))

    edges.sort(key=lambda e: e[2], reverse=True)

    return MinHashResult(
        edges=edges,
        pairs_evaluated=pairs_evaluated,
        candidates_found=candidates_found,
        lsh_threshold=config.lsh_threshold,
    )
