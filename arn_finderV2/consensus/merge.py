"""Sequence merging utilities for consensus building."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class MergeResult:
    merged: str
    overlap: int
    identity: float


def merge_sequences(
    seq_a: str,
    seq_b: str,
    min_overlap: int,
    min_identity: float,
) -> MergeResult | None:
    """Attempt to merge B into A; returns merged sequence if overlap meets criteria."""

    best_result: MergeResult | None = None
    max_overlap = min(len(seq_a), len(seq_b))
    for overlap in range(max_overlap, min_overlap - 1, -1):
        suffix = seq_a[-overlap:]
        prefix = seq_b[:overlap]
        identity = _overlap_identity(suffix, prefix)
        if identity >= min_identity:
            resolved = _resolve_conflicts(suffix, prefix)
            merged = seq_a[:-overlap] + resolved + seq_b[overlap:]
            best_result = MergeResult(merged=merged, overlap=overlap, identity=identity)
            break
    if best_result:
        return best_result

    # Try opposite direction (B -> A)
    for overlap in range(max_overlap, min_overlap - 1, -1):
        suffix = seq_b[-overlap:]
        prefix = seq_a[:overlap]
        identity = _overlap_identity(suffix, prefix)
        if identity >= min_identity:
            resolved = _resolve_conflicts(suffix, prefix)
            merged = seq_b[:-overlap] + resolved + seq_a[overlap:]
            return MergeResult(merged=merged, overlap=overlap, identity=identity)
    return None


def _overlap_identity(suffix: str, prefix: str) -> float:
    matches = 0
    length = len(suffix)
    if length == 0:
        return 0.0
    for a, b in zip(suffix, prefix):
        if a == "N" or b == "N":
            continue
        if a == b:
            matches += 1
    effective_length = sum(1 for a, b in zip(suffix, prefix) if not (a == "N" or b == "N"))
    if effective_length == 0:
        return 0.0
    return matches / effective_length


def _resolve_conflicts(suffix: str, prefix: str) -> str:
    resolved = []
    for a, b in zip(suffix, prefix):
        if a == b:
            resolved.append(a)
        elif a == "N":
            resolved.append(b)
        elif b == "N":
            resolved.append(a)
        else:
            resolved.append("N")
    return "".join(resolved)
