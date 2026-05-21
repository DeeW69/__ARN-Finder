"""Clustering V3 — Jaccard k-mer + progress Rich."""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

from rich.progress import Progress, BarColumn, TextColumn, TaskProgressColumn, TimeElapsedColumn

from ..io import ensure_dir, now_iso, read_fasta_records, write_csv, write_json, get_logger
from ..motifs import _canonical_form
from .union_find import UnionFind

logger = get_logger("arnfinderv3.cluster")


@dataclass(slots=True)
class ClusterRequest:
    input_fasta: Path
    out_dir: Path
    k: int = 9
    min_similarity: float = 0.35
    ignore_n: bool = True
    canonical: bool = False
    alphabet: str = "AUTO"
    max_pairs: int = 200_000
    metadata_path: Optional[Path] = None


@dataclass(slots=True)
class ClusterMember:
    record_id: str
    accession: Optional[str]
    length: int


@dataclass(slots=True)
class ClusterResult:
    edges_path: Path
    clusters_path: Path
    stats_path: Path
    manifest_path: Path
    total_records: int
    pairs_evaluated: int
    edges_kept: int
    cluster_count: int


def cluster_sequences(request: ClusterRequest) -> ClusterResult:
    out_dir = ensure_dir(request.out_dir)
    edges_path = out_dir / "cluster_edges.csv"
    clusters_path = out_dir / "clusters.jsonl"
    stats_path = out_dir / "cluster_stats.csv"
    manifest_path = out_dir / "cluster_manifest.json"

    canonical = request.canonical
    alphabet_mode = request.alphabet.upper()
    if canonical and alphabet_mode == "RNA":
        logger.warning("Mode canonique désactivé pour RNA.")
        canonical = False

    records, kmer_sets = _build_kmer_sets(
        request.input_fasta, request.k, request.ignore_n, canonical, alphabet_mode
    )
    total_records = len(records)

    if total_records <= 1:
        logger.warning("Pas assez de records pour le clustering.")
        _write_empty_outputs(edges_path, clusters_path, stats_path)
        manifest = _build_manifest(request, total_records, 0, 0, total_records,
                                   {"edges": str(edges_path), "clusters": str(clusters_path), "stats": str(stats_path)})
        write_json(manifest_path, manifest)
        return ClusterResult(edges_path, clusters_path, stats_path, manifest_path, total_records, 0, 0, total_records)

    metadata_lookup = _load_metadata(request.metadata_path)
    inverted_index = _build_inverted_index(kmer_sets)
    edges, pairs_evaluated, max_hit = _compute_edges(kmer_sets, request.min_similarity, inverted_index, request.max_pairs)
    if max_hit:
        logger.warning("Limite de paires atteinte (%s) — clustering peut être incomplet.", request.max_pairs)

    uf = UnionFind(total_records)
    for i, j, _ in edges:
        uf.union(i, j)

    clusters = _collect_clusters(uf, records, metadata_lookup, edges)

    _write_clusters_json(clusters_path, clusters)
    _write_edges_csv(edges_path, edges, records)
    _write_cluster_stats(stats_path, clusters)

    manifest = _build_manifest(
        request, total_records, pairs_evaluated, len(edges), len(clusters),
        {"edges": str(edges_path), "clusters": str(clusters_path), "stats": str(stats_path)},
    )
    write_json(manifest_path, manifest)

    logger.info(
        "Clustering : %s records, k=%s, min_sim=%.2f → %s paires évaluées, %s edges, %s clusters",
        total_records, request.k, request.min_similarity, pairs_evaluated, len(edges), len(clusters),
    )
    return ClusterResult(edges_path, clusters_path, stats_path, manifest_path, total_records, pairs_evaluated, len(edges), len(clusters))


# ── Helpers ───────────────────────────────────────────────────────────────────

def _build_kmer_sets(
    fasta_path: Path, k: int, ignore_n: bool, canonical: bool, alphabet: str
) -> tuple[list[ClusterMember], list[set[str]]]:
    records: list[ClusterMember] = []
    kmer_sets: list[set[str]] = []

    with Progress(
        TextColumn("[cyan]cluster[/cyan]"),
        BarColumn(), TaskProgressColumn(), TimeElapsedColumn(), transient=True
    ) as progress:
        task = progress.add_task("construction k-mer sets...", total=None)

        for header, sequence in read_fasta_records(fasta_path):
            progress.advance(task)
            record_id = header.split(" ", 1)[0]
            seq = sequence.strip().upper()
            records.append(ClusterMember(record_id=record_id, accession=record_id, length=len(seq)))
            kset: set[str] = set()
            if len(seq) >= k:
                for idx in range(len(seq) - k + 1):
                    kmer = seq[idx: idx + k]
                    if ignore_n and "N" in kmer:
                        continue
                    if canonical:
                        kmer = _canonical_form(kmer, alphabet)
                    kset.add(kmer)
            kmer_sets.append(kset)

    return records, kmer_sets


def _build_inverted_index(kmer_sets: list[set[str]]) -> Dict[str, list[int]]:
    index: Dict[str, list[int]] = defaultdict(list)
    for idx, kset in enumerate(kmer_sets):
        for kmer in kset:
            index[kmer].append(idx)
    return index


def _compute_edges(
    kmer_sets: list[set[str]],
    min_similarity: float,
    inverted_index: Dict[str, list[int]],
    max_pairs: int,
) -> tuple[list[tuple[int, int, float]], int, bool]:
    edges: list[tuple[int, int, float]] = []
    pairs_evaluated = 0
    max_hit = False

    for i in range(len(kmer_sets)):
        candidates: set[int] = set()
        for kmer in kmer_sets[i]:
            for other in inverted_index.get(kmer, []):
                if other > i:
                    candidates.add(other)
        for j in sorted(candidates):
            if pairs_evaluated >= max_pairs:
                max_hit = True
                break
            pairs_evaluated += 1
            sim = _jaccard(kmer_sets[i], kmer_sets[j])
            if sim >= min_similarity:
                edges.append((i, j, sim))
        if max_hit:
            break

    edges.sort(key=lambda x: x[2], reverse=True)
    return edges, pairs_evaluated, max_hit


def _jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 1.0
    union = len(a | b)
    return len(a & b) / union if union else 0.0


def _collect_clusters(
    uf: UnionFind,
    records: list[ClusterMember],
    metadata_lookup: dict[str, dict],
    edges: list[tuple[int, int, float]],
) -> list[dict]:
    groups: Dict[int, list[int]] = defaultdict(list)
    for idx in range(len(records)):
        groups[uf.find(idx)].append(idx)

    edge_lookup = {(min(i, j), max(i, j)): sim for i, j, sim in edges}
    result: list[dict] = []
    cluster_id = 1

    for indices in groups.values():
        members = []
        organisms: list[str] = []
        lengths = []
        for idx in indices:
            m = records[idx]
            meta = metadata_lookup.get(m.record_id) or metadata_lookup.get(m.accession or "", {})
            organism = meta.get("organism") if meta else None
            organisms.append(organism or "UNKNOWN")
            members.append({
                "record_id": m.record_id,
                "accession": m.accession,
                "length": m.length,
                "organism": organism,
                "taxid": meta.get("taxid") if meta else None,
            })
            lengths.append(m.length)

        sims = [
            edge_lookup[(min(indices[a], indices[b]), max(indices[a], indices[b]))]
            for a in range(len(indices))
            for b in range(a + 1, len(indices))
            if (min(indices[a], indices[b]), max(indices[a], indices[b])) in edge_lookup
        ]
        result.append({
            "cluster_id": f"C{cluster_id:04d}",
            "size": len(indices),
            "members": members,
            "mean_len": sum(lengths) / len(lengths),
            "dominant_organism": Counter(organisms).most_common(1)[0][0],
            "n_unique_organisms": len(set(organisms)),
            "max_pairwise_sim": max(sims) if sims else 0.0,
            "mean_pairwise_sim": sum(sims) / len(sims) if sims else 0.0,
        })
        cluster_id += 1

    result.sort(key=lambda c: c["size"], reverse=True)
    for idx, c in enumerate(result, start=1):
        c["cluster_id"] = f"C{idx:04d}"
    return result


def _write_clusters_json(path: Path, clusters: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for c in clusters:
            f.write(json.dumps(c) + "\n")


def _write_edges_csv(path: Path, edges: list[tuple[int, int, float]], records: list[ClusterMember]) -> None:
    rows = [[records[i].record_id, records[j].record_id, f"{sim:.4f}"] for i, j, sim in edges]
    write_csv(path, ["record_id_1", "record_id_2", "jaccard"], rows)


def _write_cluster_stats(path: Path, clusters: list[dict]) -> None:
    rows = [
        [c["cluster_id"], str(c["size"]), f"{c['mean_len']:.2f}",
         c["dominant_organism"] or "", str(c["n_unique_organisms"]),
         f"{c['max_pairwise_sim']:.4f}", f"{c['mean_pairwise_sim']:.4f}"]
        for c in clusters
    ]
    write_csv(path, ["cluster_id", "size", "mean_len", "dominant_organism",
                     "n_unique_organisms", "max_pairwise_sim", "mean_pairwise_sim"], rows)


def _build_manifest(
    request: ClusterRequest, total_records: int, pairs_evaluated: int,
    edges_kept: int, cluster_count: int, outputs: dict
) -> dict:
    return {
        "timestamp": now_iso(),
        "params": {
            "k": request.k, "min_similarity": request.min_similarity,
            "ignore_n": request.ignore_n, "canonical": request.canonical,
            "alphabet": request.alphabet, "max_pairs": request.max_pairs,
        },
        "inputs": {"fasta": str(request.input_fasta), "metadata": str(request.metadata_path) if request.metadata_path else None},
        "total_records": total_records,
        "pairs_evaluated": pairs_evaluated,
        "edges_kept": edges_kept,
        "cluster_count": cluster_count,
        "outputs": outputs,
    }


def _load_metadata(path: Optional[Path]) -> dict[str, dict]:
    mapping: dict[str, dict] = {}
    if not path or not path.exists():
        return mapping
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            payload = json.loads(line)
            key = payload.get("accession") or str(payload.get("uid"))
            if key:
                mapping[str(key)] = {
                    "organism": payload.get("organism") or payload.get("taxname"),
                    "taxid": payload.get("taxid"),
                }
    return mapping


def _write_empty_outputs(edges_path: Path, clusters_path: Path, stats_path: Path) -> None:
    write_csv(edges_path, ["record_id_1", "record_id_2", "jaccard"], [])
    clusters_path.parent.mkdir(parents=True, exist_ok=True)
    clusters_path.write_text("", encoding="utf-8")
    write_csv(stats_path, ["cluster_id", "size", "mean_len", "dominant_organism",
                            "n_unique_organisms", "max_pairwise_sim", "mean_pairwise_sim"], [])
