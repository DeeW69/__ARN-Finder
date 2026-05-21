"""
Rapport HTML — V3.1

Génère un rapport `report.html` dans le répertoire de données du run.
Lit les sorties pipeline existantes (CSV/JSONL) — pas de dépendance à DuckDB.

Usage :
    from arn_finder_v3.report import generate_report
    generate_report(data_dir=Path("data"), query="BRCA1", out_path=Path("data/report.html"))
"""

from __future__ import annotations

import csv
import json
import uuid
from pathlib import Path
from typing import Any, Optional

from ..io import get_logger, now_iso

logger = get_logger("arnfinderv3.report")

# Chemin vers le template HTML (même package)
_TEMPLATE = Path(__file__).parent / "template.html"


def generate_report(
    data_dir: Path,
    query: str,
    out_path: Optional[Path] = None,
    top_kmers: int = 20,
    top_clusters: int = 15,
    top_structures: int = 20,
) -> Path:
    """
    Lit les sorties du pipeline dans `data_dir` et génère un rapport HTML complet.

    Paramètres
    ----------
    data_dir    : répertoire contenant raw/, filtered/, motifs/, clusters/, etc.
    query       : requête NCBI d'origine (pour l'affichage)
    out_path    : chemin du HTML généré (défaut : data_dir/report.html)
    top_kmers   : nombre de k-mers à afficher
    top_clusters: nombre de clusters à afficher
    top_structures: nombre de structures à afficher

    Retourne
    --------
    Path vers le fichier HTML généré.
    """
    try:
        from jinja2 import Environment, FileSystemLoader, select_autoescape
    except ImportError:
        raise ImportError("jinja2 non installé. Lance : pip install jinja2>=3.1") from None

    data_dir = Path(data_dir)
    out_path = out_path or (data_dir / "report.html")
    run_id = str(uuid.uuid4())[:8]
    generated_at = now_iso()

    # ── Collecte des données ──────────────────────────────────────────────────
    ctx = _build_context(
        data_dir, query, run_id, generated_at,
        top_kmers, top_clusters, top_structures,
    )

    # ── Rendu Jinja2 ─────────────────────────────────────────────────────────
    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE.parent)),
        autoescape=select_autoescape(["html"]),
    )
    template = env.get_template(_TEMPLATE.name)
    html = template.render(**ctx)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")
    logger.info("Rapport HTML généré → %s", out_path)
    return out_path


# ── Construction du contexte ──────────────────────────────────────────────────

def _build_context(
    data_dir: Path, query: str, run_id: str, generated_at: str,
    n_kmers: int, n_clusters: int, n_structures: int,
) -> dict[str, Any]:
    stats = _collect_stats(data_dir)
    pipeline_stages = _collect_pipeline_stages(data_dir)
    filter_reasons = _collect_filter_reasons(data_dir)
    top_kmer_rows = _collect_top_kmers(data_dir, n_kmers)
    cluster_rows = _collect_clusters(data_dir, n_clusters)
    structure_rows = _collect_structures(data_dir, n_structures)

    return {
        "query": query,
        "run_id": run_id,
        "generated_at": generated_at,
        "stats": stats,
        "pipeline_stages": pipeline_stages,
        "filter_reasons": filter_reasons,
        "top_kmers": top_kmer_rows,
        "clusters": cluster_rows,
        "structure_records": structure_rows,
    }


def _collect_stats(data_dir: Path) -> dict[str, Any]:
    stats: dict[str, Any] = {
        "fetched": 0, "kept": 0, "filter_pct": 0,
        "clusters": 0, "consensus": 0,
        "structures": 0, "structure_backend": "—",
        "min_sim": "—",
    }

    # Fetch count
    ids_file = data_dir / "raw" / "ids.txt"
    if ids_file.exists():
        stats["fetched"] = sum(1 for _ in ids_file.open())

    # Filter manifest
    filt_manifest = data_dir / "filtered" / "filter_manifest.json"
    if filt_manifest.exists():
        m = json.loads(filt_manifest.read_text())
        stats["fetched"] = m.get("total_records", stats["fetched"])
        stats["kept"] = m.get("kept_records", 0)
        total = stats["fetched"] or 1
        stats["filter_pct"] = round(stats["kept"] / total * 100, 1)

    # Cluster manifest
    cluster_manifest = data_dir / "clusters" / "cluster_manifest.json"
    if cluster_manifest.exists():
        m = json.loads(cluster_manifest.read_text())
        stats["clusters"] = m.get("cluster_count", 0)
        stats["min_sim"] = m.get("params", {}).get("min_similarity", "—")

    # Consensus manifest
    cons_manifest = data_dir / "consensus" / "consensus_manifest.json"
    if cons_manifest.exists():
        m = json.loads(cons_manifest.read_text())
        stats["consensus"] = m.get("nb_consensus_ok", 0)

    # Structure manifest
    struct_manifest = data_dir / "structure" / "structure_manifest.json"
    if struct_manifest.exists():
        m = json.loads(struct_manifest.read_text())
        stats["structures"] = m.get("computed", 0)
        stats["structure_backend"] = m.get("backend", "—")

    return stats


def _collect_pipeline_stages(data_dir: Path) -> list[dict]:
    stage_manifests = [
        ("fetch",     data_dir / "raw" / "manifest.json"),
        ("filter",    data_dir / "filtered" / "filter_manifest.json"),
        ("motifs",    data_dir / "motifs" / "motifs_summary.json"),
        ("cluster",   data_dir / "clusters" / "cluster_manifest.json"),
        ("consensus", data_dir / "consensus" / "consensus_manifest.json"),
        ("structure", data_dir / "structure" / "structure_manifest.json"),
        ("blast",     data_dir / "blast" / "blast_manifest.json"),
        ("export",    data_dir / "export" / "export_manifest.json"),
    ]
    stages = []
    for name, path in stage_manifests:
        if path.exists():
            m = json.loads(path.read_text())
            elapsed = m.get("elapsed_s") or m.get("elapsed") or ""
            elapsed_str = f"{elapsed:.1f}s" if isinstance(elapsed, (int, float)) else str(elapsed)
            stages.append({"name": name, "status": "ok", "detail": "", "elapsed": elapsed_str})
        else:
            stages.append({"name": name, "status": "skipped", "detail": "", "elapsed": "—"})
    return stages


def _collect_filter_reasons(data_dir: Path) -> list[dict]:
    manifest = data_dir / "filtered" / "filter_manifest.json"
    if not manifest.exists():
        return []
    m = json.loads(manifest.read_text())
    reasons = m.get("reason_counts", {})
    total = max(sum(reasons.values()), 1)
    rows = []
    for reason, count in sorted(reasons.items(), key=lambda x: -x[1]):
        if reason == "ok":
            continue
        rows.append({
            "reason": reason,
            "count": count,
            "pct": round(count / total * 100, 1),
        })
    return rows


def _collect_top_kmers(data_dir: Path, n: int) -> list[dict]:
    csv_path = data_dir / "motifs" / "motifs.csv"
    if not csv_path.exists():
        return []
    rows = []
    with csv_path.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            if i >= n:
                break
            rows.append({
                "kmer": row.get("kmer", ""),
                "count": row.get("count", "0"),
                "freq": f"{float(row.get('freq', 0)):.4f}",
                "n_records": row.get("n_records", "0"),
                "presence_pct": min(
                    round(int(row.get("n_records", 0)) / max(1, 1) * 100, 1), 100
                ),
            })
    return rows


def _collect_clusters(data_dir: Path, n: int) -> list[dict]:
    jsonl = data_dir / "clusters" / "clusters.jsonl"
    struct_data = _load_structure_map(data_dir)
    if not jsonl.exists():
        return []
    rows = []
    with jsonl.open(encoding="utf-8") as f:
        for i, line in enumerate(f):
            if i >= n:
                break
            line = line.strip()
            if not line:
                continue
            c = json.loads(line)
            # MFE moyen si structure disponible
            mfes = [
                struct_data[m["record_id"]]["mfe"]
                for m in c.get("members", [])
                if m["record_id"] in struct_data
            ]
            mean_mfe = round(sum(mfes) / len(mfes), 2) if mfes else None
            rows.append({
                "cluster_id": c["cluster_id"],
                "size": c["size"],
                "organism": c.get("dominant_organism") or "—",
                "mean_sim": f"{c.get('mean_pairwise_sim', 0):.3f}",
                "mean_mfe": mean_mfe,
            })
    return rows


def _collect_structures(data_dir: Path, n: int) -> list[dict]:
    csv_path = data_dir / "structure" / "structures.csv"
    if not csv_path.exists():
        return []
    rows = []
    with csv_path.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            if i >= n:
                break
            if row.get("skipped", "False") == "True":
                continue
            try:
                mfe = float(row.get("mfe_kcal_mol", 0))
            except ValueError:
                mfe = 0.0
            rows.append({
                "record_id": row["record_id"],
                "length": row["length"],
                "gc_frac": f"{float(row.get('gc_frac', 0)):.3f}",
                "mfe": round(mfe, 2),
                "dot_bracket": row.get("dot_bracket", ""),
            })
    return rows


def _load_structure_map(data_dir: Path) -> dict[str, dict]:
    csv_path = data_dir / "structure" / "structures.csv"
    result: dict[str, dict] = {}
    if not csv_path.exists():
        return result
    with csv_path.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                result[row["record_id"]] = {"mfe": float(row.get("mfe_kcal_mol", 0))}
            except (KeyError, ValueError):
                pass
    return result
