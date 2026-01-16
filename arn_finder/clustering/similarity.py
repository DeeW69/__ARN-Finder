"""K-mer based clustering via Jaccard similarity."""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

from ..io import ensure_dir, now_iso, read_fasta_records, write_csv, write_json
from ..motifs.kmer_scan import _canonical_form  # reuse canonical helper
from .union_find import UnionFind


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
    metadata_path: Path | None = None


@dataclass(slots=True)
class ClusterMember:
    record_id: str
    accession: str | None
    length: int
    organism: str | None = None
    taxid: int | None = None


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


def cluster_sequences(request: ClusterRequest, logger) -> ClusterResult:
    out_dir = ensure_dir(request.out_dir)
    edges_path = out_dir / "cluster_edges.csv"
    clusters_path = out_dir / "clusters.jsonl"
    stats_path = out_dir / "cluster_stats.csv"
    manifest_path = out_dir / "cluster_manifest.json"

    canonical = request.canonical
    alphabet_mode = request.alphabet.upper()
    if canonical and alphabet_mode == "RNA":
        logger.warning("Canonical reverse-complement disabled for RNA alphabet.")
        canonical = False

    records, kmer_sets = _build_kmer_sets(
        request.input_fasta,
        request.k,
        request.ignore_n,
        canonical,
        alphabet_mode,
        logger,
    )
    total_records = len(records)

    if total_records <= 1:
        logger.warning("Not enough records for clustering. Outputs will be minimal.")
        _write_empty_outputs(edges_path, clusters_path, stats_path)
        manifest = _build_manifest(
            request,
            total_records,
            pairs_evaluated=0,
            edges_kept=0,
            cluster_count=total_records,
            outputs={
                "edges": str(edges_path),
                "clusters": str(clusters_path),
                "stats": str(stats_path),
            },
        )
        write_json(manifest_path, manifest)
        return ClusterResult(edges_path, clusters_path, stats_path, manifest_path, total_records, 0, 0, total_records)

    metadata_lookup = _load_metadata(request.metadata_path)
    inverted_index = _build_inverted_index(kmer_sets)
    uf = UnionFind(total_records)
    edges, pairs_evaluated, max_pairs_hit = _compute_edges(
        kmer_sets,
        request.min_similarity,
        inverted_index,
        request.max_pairs,
    )
    if max_pairs_hit:
        logger.warning("Reached max pair limit (%s); clustering may be incomplete.", request.max_pairs)

    for i, j, _ in edges:
        uf.union(i, j)

    clusters = _collect_clusters(uf, records, metadata_lookup, edges)
    clusters.sort(key=lambda c: len(c["members"]), reverse=True)

    _write_clusters_json(clusters_path, clusters)
    _write_edges_csv(edges_path, edges, records)
    _write_cluster_stats(stats_path, clusters)

    manifest = _build_manifest(
        request,
        total_records,
        pairs_evaluated=pairs_evaluated,
        edges_kept=len(edges),
        cluster_count=len(clusters),
        outputs={
            "edges": str(edges_path),
            "clusters": str(clusters_path),
            "stats": str(stats_path),
        },
    )
    write_json(manifest_path, manifest)

    logger.info(
        "Clustering: %s records, k=%s, min_sim=%.2f, evaluated=%s pairs, kept=%s edges, clusters=%s",
        total_records,
        request.k,
        request.min_similarity,
        pairs_evaluated,
        len(edges),
        len(clusters),
    )
    logger.info("Clusters JSONL -> %s", clusters_path)
    logger.info("Cluster edges CSV -> %s", edges_path)
    logger.info("Cluster stats -> %s", stats_path)

    return ClusterResult(edges_path, clusters_path, stats_path, manifest_path, total_records, pairs_evaluated, len(edges), len(clusters))


def _build_kmer_sets(
    fasta_path: Path,
    k: int,
    ignore_n: bool,
    canonical: bool,
    alphabet: str,
    logger,
) -> tuple[list[ClusterMember], list[set[str]]]:
    records: list[ClusterMember] = []
    kmer_sets: list[set[str]] = []
    for header, sequence in read_fasta_records(fasta_path):
        record_id = header.split(" ", 1)[0]
        accession = record_id
        seq = sequence.strip().upper()
        record = ClusterMember(record_id=record_id, accession=accession, length=len(seq))
        records.append(record)
        kset = set()
        if len(seq) < k:
            logger.warning("Sequence %s shorter than k=%s; produces no k-mers.", record_id, k)
        else:
            for idx in range(len(seq) - k + 1):
                kmer = seq[idx : idx + k]
                if ignore_n and "N" in kmer:
                    continue
                if canonical:
                    kmer = _canonical_form(kmer, alphabet)
                kset.add(kmer)
        if not kset:
            logger.warning("Sequence %s has no valid k-mers.", record_id)
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
    n = len(kmer_sets)
    pairs_evaluated = 0
    max_pairs_hit = False

    for i in range(n):
        candidates: set[int] = set()
        for kmer in kmer_sets[i]:
            for other in inverted_index.get(kmer, []):
                if other <= i:
                    continue
                candidates.add(other)
        for j in sorted(candidates):
            if pairs_evaluated >= max_pairs:
                max_pairs_hit = True
                break
            pairs_evaluated += 1
            sim = _jaccard(kmer_sets[i], kmer_sets[j])
            if sim >= min_similarity:
                edges.append((i, j, sim))
        if max_pairs_hit:
            break
    edges.sort(key=lambda x: x[2], reverse=True)
    return edges, pairs_evaluated, max_pairs_hit


def _jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 1.0
    intersection = len(a & b)
    union = len(a | b)
    return intersection / union if union else 0.0


def _collect_clusters(
    uf: UnionFind,
    records: list[ClusterMember],
    metadata_lookup: dict[str, dict],
    edges: list[tuple[int, int, float]],
) -> list[dict]:
    groups: Dict[int, list[int]] = defaultdict(list)
    for idx in range(len(records)):
        root = uf.find(idx)
        groups[root].append(idx)
    edge_lookup = {(min(i, j), max(i, j)): sim for i, j, sim in edges}
    result: list[dict] = []
    cluster_id = 1
    for indices in groups.values():
        members = []
        organisms: list[str] = []
        lengths = []
        for idx in indices:
            member = records[idx]
            meta = metadata_lookup.get(member.record_id) or metadata_lookup.get(member.accession or "", {})
            organism = meta.get("organism") if meta else None
            taxid = meta.get("taxid") if meta else None
            if organism:
                organisms.append(organism)
            else:
                organisms.append("UNKNOWN")
            member_dict = {
                "record_id": member.record_id,
                "accession": member.accession,
                "length": member.length,
                "organism": organism,
                "taxid": taxid,
            }
            members.append(member_dict)
            lengths.append(member.length)
        sims = [
            edge_lookup[(min(i, j), max(i, j))]
            for idx_i in range(len(indices))
            for idx_j in range(idx_i + 1, len(indices))
            for i, j in [(indices[idx_i], indices[idx_j])]
            if (min(i, j), max(i, j)) in edge_lookup
        ]
        cluster_info = {
            "cluster_id": f"C{cluster_id:04d}",
            "size": len(indices),
            "members": members,
            "mean_len": sum(lengths) / len(lengths) if lengths else 0,
            "dominant_organism": _dominant_organism(organisms),
            "n_unique_organisms": len(set(organisms)),
            "max_pairwise_sim": max(sims) if sims else 0.0,
            "mean_pairwise_sim": sum(sims) / len(sims) if sims else 0.0,
        }
        result.append(cluster_info)
        cluster_id += 1
    result.sort(key=lambda c: c["size"], reverse=True)
    for idx, cluster in enumerate(result, start=1):
        cluster["cluster_id"] = f"C{idx:04d}"
    return result


def _dominant_organism(organisms: list[str]) -> str | None:
    if not organisms:
        return None
    counts = Counter(organisms)
    organism, _ = counts.most_common(1)[0]
    return organism


def _write_clusters_json(path: Path, clusters: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for cluster in clusters:
            handle.write(json.dumps(cluster) + "\n")


def _write_edges_csv(path: Path, edges: list[tuple[int, int, float]], records: list[ClusterMember]) -> None:
    rows = [
        [records[i].record_id, records[j].record_id, f"{sim:.4f}"]
        for i, j, sim in edges
    ]
    write_csv(path, ["record_id_1", "record_id_2", "jaccard"], rows)


def _write_cluster_stats(path: Path, clusters: list[dict]) -> None:
    rows = []
    for cluster in clusters:
        rows.append(
            [
                cluster["cluster_id"],
                str(cluster["size"]),
                f"{cluster['mean_len']:.2f}",
                cluster["dominant_organism"] or "",
                str(cluster["n_unique_organisms"]),
                f"{cluster['max_pairwise_sim']:.4f}",
                f"{cluster['mean_pairwise_sim']:.4f}",
            ]
        )
    write_csv(
        path,
        [
            "cluster_id",
            "size",
            "mean_len",
            "dominant_organism",
            "n_unique_organisms",
            "max_pairwise_sim",
            "mean_pairwise_sim",
        ],
        rows,
    )


def _build_manifest(
    request: ClusterRequest,
    total_records: int,
    pairs_evaluated: int,
    edges_kept: int,
    cluster_count: int,
    outputs: dict[str, str],
) -> dict:
    return {
        "timestamp": now_iso(),
        "params": {
            "k": request.k,
            "min_similarity": request.min_similarity,
            "ignore_n": request.ignore_n,
            "canonical": request.canonical,
            "alphabet": request.alphabet,
            "max_pairs": request.max_pairs,
        },
        "inputs": {
            "fasta": str(request.input_fasta),
            "metadata": str(request.metadata_path) if request.metadata_path else None,
        },
        "total_records": total_records,
        "pairs_evaluated": pairs_evaluated,
        "edges_kept": edges_kept,
        "cluster_count": cluster_count,
        "outputs": outputs,
    }


def _load_metadata(path: Path | None) -> dict[str, dict]:
    mapping: dict[str, dict] = {}
    if not path or not path.exists():
        return mapping
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
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
    write_csv(
        stats_path,
        [
            "cluster_id",
            "size",
            "mean_len",
            "dominant_organism",
            "n_unique_organisms",
            "max_pairwise_sim",
            "mean_pairwise_sim",
        ],
        [],
    )
