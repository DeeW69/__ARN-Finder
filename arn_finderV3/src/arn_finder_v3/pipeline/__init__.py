"""Orchestrateur de pipeline V3 — enchaîne toutes les étapes automatiquement."""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
from rich.table import Table

from ..io import get_logger

logger = get_logger("arnfinderv3.pipeline")
console = Console()


@dataclass
class PipelineResult:
    stage_results: dict[str, Any] = field(default_factory=dict)
    errors: dict[str, Exception] = field(default_factory=dict)
    elapsed: dict[str, float] = field(default_factory=dict)

    @property
    def success(self) -> bool:
        return len(self.errors) == 0


@dataclass
class PipelineConfig:
    # Requête NCBI
    query: str

    # Répertoire de base
    out_base: Path = Path("data")

    # Étapes activées
    run_fetch: bool = True
    run_filter: bool = True
    run_motifs: bool = True
    run_cluster: bool = True
    run_consensus: bool = True
    run_blast: bool = False     # Désactivé par défaut (très long)
    run_export: bool = True

    # Paramètres fetch
    fetch_limit: int = 200
    fetch_db: str = "nucleotide"

    # Paramètres filter
    filter_min_len: int = 200
    filter_max_n_frac: float = 0.05
    filter_min_gc: float = 0.0
    filter_max_gc: float = 1.0
    filter_dedupe: bool = False

    # Paramètres motifs
    motif_k: int = 9
    motif_top: int = 200
    motif_canonical: bool = False

    # Paramètres cluster
    cluster_k: int = 9
    cluster_min_sim: float = 0.35
    cluster_method: str = "jaccard"   # "jaccard" | "minhash" (V3.3)
    cluster_minhash_hashes: int = 128
    cluster_minhash_bands: int = 16

    # Paramètres consensus
    consensus_min_cluster_size: int = 2
    consensus_method: str = "kmer"

    # Structure secondaire (V3.1)
    run_structure: bool = True
    structure_max_len: int = 2000
    structure_temperature: float = 37.0

    # Rapport HTML (V3.1)
    run_report: bool = True

    # Codon usage (V3.2)
    run_codon_usage: bool = True
    codon_reading_frame: int = 0
    codon_min_codons: int = 10

    # Comportement
    stop_on_error: bool = True


def run_pipeline(cfg: PipelineConfig) -> PipelineResult:
    """
    Lance le pipeline complet étape par étape.

    Killer-feature V3 : une seule commande `arnfinderv3 run "..."` remplace
    les 6 commandes manuelles de V1/V2.
    """
    from ..fetch import FetchRequest, fetch_sequences
    from ..filters import FilterRequest, filter_sequences
    from ..motifs import MotifRequest, compute_motifs
    from ..clustering import ClusterRequest, cluster_sequences
    from ..consensus import ConsensusRequest, build_consensus
    from ..blast import BlastRequest, run_blast
    from ..secondary_structure import StructureRequest, compute_secondary_structure
    from ..codon_usage import CodonUsageRequest, compute_codon_usage
    from ..report import generate_report
    from ..export import export_ml_ready

    # Imports config (facultatif, peut crasher si .env absent)
    try:
        from ..config import ncbi, blast as blast_cfg, dirs
        email = ncbi.email
        api_key = ncbi.api_key
        tool = ncbi.tool
    except Exception:
        email = ""
        api_key = None
        tool = "arnfinder_v3"

    result = PipelineResult()
    base = cfg.out_base

    # Chemins intermédiaires résolus dynamiquement au fil du pipeline
    raw_fasta: Optional[Path] = None
    filtered_fasta: Optional[Path] = None
    filtered_metadata: Optional[Path] = None
    motifs_by_record: Optional[Path] = None
    clusters_path: Optional[Path] = None
    consensus_fasta: Optional[Path] = None
    features_csv: Optional[Path] = None
    structure_csv: Optional[Path] = None

    stages = [
        ("fetch",        cfg.run_fetch),
        ("filter",       cfg.run_filter),
        ("motifs",       cfg.run_motifs),
        ("cluster",      cfg.run_cluster),
        ("consensus",    cfg.run_consensus),
        ("structure",    cfg.run_structure),
        ("codon_usage",  cfg.run_codon_usage),
        ("blast",        cfg.run_blast),
        ("export",       cfg.run_export),
        ("report",       cfg.run_report),
    ]

    for stage_name, enabled in stages:
        if not enabled:
            console.print(f"  [dim]⊘  {stage_name} ignoré[/dim]")
            continue

        t0 = time.perf_counter()
        console.print(f"  [bold cyan]▶  {stage_name}[/bold cyan]")

        try:
            if stage_name == "fetch":
                req = FetchRequest(
                    query=cfg.query,
                    db=cfg.fetch_db,
                    limit=cfg.fetch_limit,
                    out_dir=base / "raw",
                )
                res = fetch_sequences(req)
                raw_fasta = res.fasta_path
                result.stage_results["fetch"] = res

            elif stage_name == "filter":
                if not raw_fasta or not raw_fasta.exists():
                    raise FileNotFoundError(f"Pas de FASTA fetch en entrée : {raw_fasta}")
                req = FilterRequest(
                    input_fasta=raw_fasta,
                    out_dir=base / "filtered",
                    min_len=cfg.filter_min_len,
                    max_n_fraction=cfg.filter_max_n_frac,
                    min_gc=cfg.filter_min_gc,
                    max_gc=cfg.filter_max_gc,
                    dedupe=cfg.filter_dedupe,
                    metadata_path=raw_fasta.parent / "metadata.jsonl",
                )
                res = filter_sequences(req)
                filtered_fasta = res.fasta_path
                filtered_metadata = res.metadata_path
                result.stage_results["filter"] = res

            elif stage_name == "motifs":
                if not filtered_fasta or not filtered_fasta.exists():
                    raise FileNotFoundError(f"Pas de FASTA filtré en entrée : {filtered_fasta}")
                req = MotifRequest(
                    input_fasta=filtered_fasta,
                    out_dir=base / "motifs",
                    k=cfg.motif_k,
                    top=cfg.motif_top,
                    canonical=cfg.motif_canonical,
                )
                res = compute_motifs(req)
                motifs_by_record = res.motifs_by_record_csv
                features_csv = base / "motifs" / "motifs_by_record.csv"
                result.stage_results["motifs"] = res

            elif stage_name == "cluster":
                if not filtered_fasta or not filtered_fasta.exists():
                    raise FileNotFoundError(f"Pas de FASTA filtré en entrée : {filtered_fasta}")
                req = ClusterRequest(
                    input_fasta=filtered_fasta,
                    out_dir=base / "clusters",
                    k=cfg.cluster_k,
                    min_similarity=cfg.cluster_min_sim,
                    metadata_path=filtered_metadata,
                    method=cfg.cluster_method,
                    minhash_hashes=cfg.cluster_minhash_hashes,
                    minhash_bands=cfg.cluster_minhash_bands,
                )
                res = cluster_sequences(req)
                clusters_path = res.clusters_path
                result.stage_results["cluster"] = res

            elif stage_name == "consensus":
                if not filtered_fasta or not clusters_path:
                    raise FileNotFoundError("FASTA filtré ou clusters.jsonl manquant.")
                req = ConsensusRequest(
                    clusters_path=clusters_path,
                    fasta_path=filtered_fasta,
                    out_dir=base / "consensus",
                    min_cluster_size=cfg.consensus_min_cluster_size,
                    method=cfg.consensus_method,
                    metadata_path=filtered_metadata,
                )
                res = build_consensus(req)
                consensus_fasta = base / "consensus" / "consensus.fasta"
                result.stage_results["consensus"] = res

            elif stage_name == "structure":
                target = consensus_fasta or filtered_fasta
                if not target or not target.exists():
                    raise FileNotFoundError(f"Pas de FASTA pour la prédiction de structure : {target}")
                req = StructureRequest(
                    input_fasta=target,
                    out_dir=base / "structure",
                    max_len=cfg.structure_max_len,
                    temperature=cfg.structure_temperature,
                )
                res = compute_secondary_structure(req)
                structure_csv = res.structures_csv
                result.stage_results["structure"] = res

            elif stage_name == "codon_usage":
                target = filtered_fasta
                if not target or not target.exists():
                    raise FileNotFoundError(f"Pas de FASTA filtré pour codon usage : {target}")
                req = CodonUsageRequest(
                    input_fasta=target,
                    out_dir=base / "codon_usage",
                    reading_frame=cfg.codon_reading_frame,
                    min_codons=cfg.codon_min_codons,
                )
                res = compute_codon_usage(req)
                result.stage_results["codon_usage"] = res

            elif stage_name == "blast":
                target = consensus_fasta or filtered_fasta
                if not target or not target.exists():
                    raise FileNotFoundError(f"Pas de FASTA pour BLAST : {target}")
                req = BlastRequest(
                    input_fasta=target,
                    out_dir=base / "blast",
                    email=email,
                    tool=tool,
                    api_key=api_key,
                )
                res = run_blast(req)
                result.stage_results["blast"] = res

            elif stage_name == "export":
                # On utilise filter_report.csv comme features si motifs_by_record absent
                feat = features_csv or (base / "filtered" / "filter_report.csv")
                meta = filtered_metadata or (base / "raw" / "metadata.jsonl")
                if not feat or not feat.exists():
                    raise FileNotFoundError(f"Pas de features CSV pour export : {feat}")
                if not meta or not meta.exists():
                    raise FileNotFoundError(f"Pas de metadata JSONL pour export : {meta}")
                res = export_ml_ready(
                    features_csv=feat,
                    metadata_jsonl=meta,
                    out_dir=base / "export",
                    motifs_csv=motifs_by_record,
                    validate=False,   # On ne force pas Pandera sur le run pipeline
                )
                result.stage_results["export"] = res

            elif stage_name == "report":
                report_path = generate_report(
                    data_dir=base,
                    query=cfg.query,
                    out_path=base / "report.html",
                )
                result.stage_results["report"] = {"path": str(report_path)}

            elapsed = time.perf_counter() - t0
            result.elapsed[stage_name] = elapsed
            console.print(f"  [green]✓  {stage_name}[/green] [dim]{elapsed:.1f}s[/dim]")

        except Exception as exc:
            elapsed = time.perf_counter() - t0
            result.elapsed[stage_name] = elapsed
            result.errors[stage_name] = exc
            console.print(f"  [red]✗  {stage_name} : {exc}[/red]")
            logger.exception("Erreur étape '%s'", stage_name)
            if cfg.stop_on_error:
                break

    return result


def print_pipeline_summary(result: PipelineResult) -> None:
    """Affiche un tableau Rich de résumé du pipeline."""
    table = Table(title="Résumé du pipeline", show_lines=True)
    table.add_column("Étape", style="cyan")
    table.add_column("Statut")
    table.add_column("Durée (s)", justify="right")

    for stage, elapsed in result.elapsed.items():
        if stage in result.errors:
            status = f"[red]✗ {result.errors[stage]}[/red]"
        else:
            status = "[green]✓ OK[/green]"
        table.add_row(stage, status, f"{elapsed:.1f}")

    console.print(table)
