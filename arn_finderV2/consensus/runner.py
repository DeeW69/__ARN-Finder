"""Consensus builder orchestrating cluster-based merging."""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict

from ..io import ensure_dir, now_iso, read_fasta_records, write_csv, write_fasta_record, write_json
from .merge import merge_sequences


@dataclass(slots=True)
class ConsensusRequest:
    clusters_path: Path
    fasta_path: Path
    out_dir: Path
    min_cluster_size: int = 2
    method: str = "kmer"
    k: int = 9
    min_overlap: int = 80
    min_identity: float = 0.97
    max_n_frac: float = 0.05
    alphabet: str = "AUTO"
    canonical: bool = True
    metadata_path: Path | None = None


def build_consensus(request: ConsensusRequest, logger) -> dict:
    out_dir = ensure_dir(request.out_dir)
    fasta_out = out_dir / "consensus.fasta"
    stats_out = out_dir / "consensus_stats.csv"
    manifest_out = out_dir / "consensus_manifest.json"

    sequence_lookup = _load_sequences(request.fasta_path)
    clusters = _load_clusters(request.clusters_path)
    metadata_lookup = _load_metadata(request.metadata_path)

    stats_rows = []
    total_clusters = len(clusters)
    attempted = ok = failed = 0

    with fasta_out.open("w", encoding="utf-8") as fasta_handle:
        for cluster in clusters:
            cluster_id = cluster["cluster_id"]
            members = cluster.get("members", [])
            if len(members) < request.min_cluster_size:
                logger.info("Skipping cluster %s (size=%s < min_cluster_size)", cluster_id, len(members))
                continue
            attempted += 1
            sequences = []
            for member in members:
                record_id = member["record_id"]
                seq = sequence_lookup.get(record_id)
                if not seq:
                    logger.warning("Sequence %s not found in FASTA; skipping in consensus.", record_id)
                    continue
                sequences.append(seq)
            if not sequences:
                logger.warning("No sequences available for cluster %s", cluster_id)
                failed += 1
                stats_rows.append(_stats_row(cluster_id, len(members), "", "failed", "no_sequences", request))
                continue
            consensus_seq, status, note = _build_cluster_consensus(sequences, request, logger)
            if status == "failed":
                failed += 1
            else:
                ok += 1
            header = _build_header(
                cluster_id,
                len(members),
                cluster.get("dominant_organism"),
                members,
                metadata_lookup,
                request,
            )
            write_fasta_record(fasta_handle, header, consensus_seq)
            stats_rows.append(_stats_row(cluster_id, len(members), consensus_seq, status, note, request))

    write_csv(
        stats_out,
        [
            "cluster_id",
            "cluster_size",
            "consensus_len",
            "n_count",
            "n_frac",
            "gc_percent",
            "method",
            "k",
            "min_overlap",
            "min_identity",
            "status",
            "note",
        ],
        stats_rows,
    )

    manifest = {
        "timestamp": now_iso(),
        "params": {
            "method": request.method,
            "k": request.k,
            "min_overlap": request.min_overlap,
            "min_identity": request.min_identity,
            "min_cluster_size": request.min_cluster_size,
            "max_n_frac": request.max_n_frac,
            "canonical": request.canonical,
        },
        "inputs": {
            "clusters": str(request.clusters_path),
            "fasta": str(request.fasta_path),
            "metadata": str(request.metadata_path) if request.metadata_path else None,
        },
        "nb_clusters_total": total_clusters,
        "nb_clusters_attempted": attempted,
        "nb_consensus_ok": ok,
        "nb_consensus_failed": failed,
        "outputs": {
            "consensus_fasta": str(fasta_out),
            "consensus_stats": str(stats_out),
        },
    }
    write_json(manifest_out, manifest)

    logger.info(
        "Consensus: total=%s, attempted=%s, ok=%s, failed=%s",
        total_clusters,
        attempted,
        ok,
        failed,
    )
    logger.info("Consensus FASTA -> %s", fasta_out)
    logger.info("Consensus stats -> %s", stats_out)

    return manifest


def _load_sequences(fasta_path: Path) -> dict[str, str]:
    return {header.split(" ", 1)[0]: seq.strip().upper() for header, seq in read_fasta_records(fasta_path)}


def _load_clusters(path: Path) -> list[dict]:
    clusters = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                clusters.append(json.loads(line))
    return clusters


def _load_metadata(path: Path | None) -> dict[str, dict]:
    if not path or not path.exists():
        return {}
    data: dict[str, dict] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            payload = json.loads(line)
            key = payload.get("accession") or payload.get("uid")
            if key:
                data[str(key)] = payload
    return data


def _build_cluster_consensus(sequences: list[str], request: ConsensusRequest, logger):
    sequences.sort(key=len, reverse=True)
    if request.method == "simple":
        consensus = sequences[0]
        status = "ok"
        note = "representative_longest"
    else:
        consensus = sequences[0]
        merged_any = False
        for seq in sequences[1:]:
            result = merge_sequences(consensus, seq, request.min_overlap, request.min_identity)
            if not result:
                logger.debug("Failed to merge sequence of len %s into cluster", len(seq))
                continue
            consensus = result.merged
            merged_any = True
        if not merged_any and len(sequences) == 1:
            consensus = sequences[0]
        if merged_any or len(sequences) == 1:
            status = "ok"
            note = ""
        else:
            status = "failed"
            note = "insufficient_overlap"
    n_frac = consensus.count("N") / len(consensus) if consensus else 1.0
    if n_frac > request.max_n_frac:
        return consensus, "low_quality", "high_N_fraction"
    return consensus, status, note


def _stats_row(cluster_id: str, size: int, sequence: str, status: str, note: str, request: ConsensusRequest) -> list[str]:
    length = len(sequence)
    n_count = sequence.count("N")
    n_frac = n_count / length if length else 0.0
    gc_percent = 0.0
    if length:
        gc = sequence.count("G") + sequence.count("C")
        gc_percent = gc / length * 100
    return [
        cluster_id,
        str(size),
        str(length),
        str(n_count),
        f"{n_frac:.4f}",
        f"{gc_percent:.2f}",
        request.method,
        str(request.k),
        str(request.min_overlap),
        f"{request.min_identity:.2f}",
        status,
        note,
    ]


def _build_header(
    cluster_id: str,
    size: int,
    dominant_organism: str | None,
    members: list[dict],
    metadata_lookup: dict[str, dict],
    request: ConsensusRequest,
) -> str:
    organism = dominant_organism
    if not organism and metadata_lookup:
        counts = Counter()
        for member in members:
            rid = member["record_id"]
            meta = metadata_lookup.get(rid)
            if meta and (org := meta.get("organism")):
                counts[org] += 1
        if counts:
            organism = counts.most_common(1)[0][0]
    organism = organism or "UNKNOWN"
    return (
        f"CONS_{cluster_id} size={size} dominant_organism=\"{organism}\" "
        f"method={request.method} k={request.k}"
    )
