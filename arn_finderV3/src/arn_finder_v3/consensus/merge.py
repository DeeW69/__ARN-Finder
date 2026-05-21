"""Fusion de séquences pour la construction de consensus."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class MergeResult:
    merged: str
    overlap: int
    identity: float


def merge_sequences(seq_a: str, seq_b: str, min_overlap: int, min_identity: float) -> MergeResult | None:
    """Tente de fusionner B dans A ; retourne la séquence fusionnée si le chevauchement est suffisant."""
    best = _try_merge(seq_a, seq_b, min_overlap, min_identity)
    if best:
        return best
    return _try_merge(seq_b, seq_a, min_overlap, min_identity)


def _try_merge(left: str, right: str, min_overlap: int, min_identity: float) -> MergeResult | None:
    max_overlap = min(len(left), len(right))
    for overlap in range(max_overlap, min_overlap - 1, -1):
        suffix = left[-overlap:]
        prefix = right[:overlap]
        identity = _overlap_identity(suffix, prefix)
        if identity >= min_identity:
            resolved = _resolve_conflicts(suffix, prefix)
            merged = left[:-overlap] + resolved + right[overlap:]
            return MergeResult(merged=merged, overlap=overlap, identity=identity)
    return None


def _overlap_identity(suffix: str, prefix: str) -> float:
    effective = sum(1 for a, b in zip(suffix, prefix) if not (a == "N" or b == "N"))
    if effective == 0:
        return 0.0
    matches = sum(1 for a, b in zip(suffix, prefix) if a != "N" and b != "N" and a == b)
    return matches / effective


def _resolve_conflicts(suffix: str, prefix: str) -> str:
    result = []
    for a, b in zip(suffix, prefix):
        if a == b:
            result.append(a)
        elif a == "N":
            result.append(b)
        elif b == "N":
            result.append(a)
        else:
            result.append("N")
    return "".join(result)
