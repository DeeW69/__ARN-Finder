"""CLI V3 — Typer + Rich, toutes les commandes câblées."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.console import Console

from . import __version__
from .io import configure_logging

app = typer.Typer(
    name="arnfinderv3",
    help="ARN Finder V3 — pipeline bioinformatique ML-ready.",
    add_completion=True,
    no_args_is_help=True,
)
console = Console()


# ── Version ───────────────────────────────────────────────────────────────────

def _version_cb(value: bool) -> None:
    if value:
        console.print(f"arnfinderv3 v{__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Annotated[Optional[bool], typer.Option(
        "--version", "-v", callback=_version_cb, is_eager=True, help="Affiche la version."
    )] = None,
    verbose: Annotated[bool, typer.Option("--verbose", help="Active les logs DEBUG.")] = False,
) -> None:
    configure_logging(verbose)


# ── run : pipeline complet ────────────────────────────────────────────────────

@app.command("run")
def cmd_run(
    query: Annotated[str, typer.Argument(help="Requête NCBI Entrez")],
    out: Annotated[Path, typer.Option("--out", "-o", help="Répertoire de sortie")] = Path("data"),
    limit: Annotated[int, typer.Option(help="Nombre max de séquences")] = 200,
    min_len: Annotated[int, typer.Option(help="Longueur minimale")] = 200,
    max_n_frac: Annotated[float, typer.Option(help="Fraction N max")] = 0.05,
    min_gc: Annotated[float, typer.Option(help="GC min (0-1)")] = 0.0,
    max_gc: Annotated[float, typer.Option(help="GC max (0-1)")] = 1.0,
    dedupe: Annotated[bool, typer.Option("--dedupe/--no-dedupe", help="Déduplication")] = False,
    k: Annotated[int, typer.Option(help="Taille des k-mers")] = 9,
    min_sim: Annotated[float, typer.Option(help="Seuil Jaccard clustering")] = 0.35,
    cluster_method: Annotated[str, typer.Option(help="Méthode clustering : jaccard | minhash")] = "jaccard",
    blast: Annotated[bool, typer.Option("--blast/--no-blast", help="Active BLAST")] = False,
    no_structure: Annotated[bool, typer.Option("--no-structure", help="Désactive structure secondaire")] = False,
    no_codon_usage: Annotated[bool, typer.Option("--no-codon-usage", help="Désactive le calcul codon usage")] = False,
    no_embeddings: Annotated[bool, typer.Option("--no-embeddings", help="Désactive le calcul d'embeddings")] = False,
    embedding_backend: Annotated[str, typer.Option(help="Backend embeddings : auto | kmer_tfidf | nucleotide_transformer | esm2")] = "auto",
    no_export: Annotated[bool, typer.Option("--no-export", help="Désactive l'export Parquet")] = False,
    no_report: Annotated[bool, typer.Option("--no-report", help="Désactive le rapport HTML")] = False,
    stop_on_error: Annotated[bool, typer.Option("--stop-on-error/--continue-on-error")] = True,
) -> None:
    """Lance le pipeline complet (fetch → filter → motifs → cluster → consensus → structure → [blast] → export → report)."""
    from .pipeline import PipelineConfig, run_pipeline, print_pipeline_summary

    cfg = PipelineConfig(
        query=query,
        out_base=out,
        fetch_limit=limit,
        filter_min_len=min_len,
        filter_max_n_frac=max_n_frac,
        filter_min_gc=min_gc,
        filter_max_gc=max_gc,
        filter_dedupe=dedupe,
        motif_k=k,
        cluster_k=k,
        cluster_min_sim=min_sim,
        cluster_method=cluster_method,
        run_blast=blast,
        run_structure=not no_structure,
        run_codon_usage=not no_codon_usage,
        run_embeddings=not no_embeddings,
        embedding_backend=embedding_backend,
        run_export=not no_export,
        run_report=not no_report,
        stop_on_error=stop_on_error,
    )

    console.rule("[bold cyan]ARN Finder V3 — Pipeline[/bold cyan]")
    result = run_pipeline(cfg)
    print_pipeline_summary(result)
    raise typer.Exit(0 if result.success else 1)


# ── fetch ─────────────────────────────────────────────────────────────────────

@app.command("fetch")
def cmd_fetch(
    query: Annotated[str, typer.Option("--query", "-q", help="Requête NCBI Entrez")],
    out: Annotated[Path, typer.Option("-o")] = Path("data/raw"),
    limit: Annotated[int, typer.Option()] = 200,
    db: Annotated[str, typer.Option()] = "nucleotide",
) -> None:
    """Télécharge des séquences depuis NCBI Entrez (aiohttp, async)."""
    from .fetch import FetchRequest, fetch_sequences

    req = FetchRequest(query=query, db=db, limit=limit, out_dir=out)
    res = fetch_sequences(req)
    console.print(f"[green]✓[/green] {res.fetched} séquences → {res.fasta_path} ({res.elapsed:.1f}s)")


# ── filter ────────────────────────────────────────────────────────────────────

@app.command("filter")
def cmd_filter(
    fasta: Annotated[Path, typer.Argument(help="FASTA en entrée")],
    out: Annotated[Path, typer.Option("-o")] = Path("data/filtered"),
    min_len: Annotated[int, typer.Option()] = 200,
    max_len: Annotated[Optional[int], typer.Option()] = None,
    max_n_frac: Annotated[float, typer.Option()] = 0.05,
    min_gc: Annotated[float, typer.Option(help="GC min (0-1)")] = 0.0,
    max_gc: Annotated[float, typer.Option(help="GC max (0-1)")] = 1.0,
    dedupe: Annotated[bool, typer.Option("--dedupe/--no-dedupe")] = False,
    metadata: Annotated[Optional[Path], typer.Option(help="JSONL metadata NCBI")] = None,
    max_per_organism: Annotated[int, typer.Option()] = 0,
) -> None:
    """Filtre les séquences par qualité (longueur, N, GC, dédoublonnage)."""
    from .filters import FilterRequest, filter_sequences

    req = FilterRequest(
        input_fasta=fasta, out_dir=out,
        min_len=min_len, max_len=max_len,
        max_n_fraction=max_n_frac,
        min_gc=min_gc, max_gc=max_gc,
        dedupe=dedupe,
        metadata_path=metadata,
        max_per_organism=max_per_organism,
    )
    res = filter_sequences(req)
    console.print(
        f"[green]✓[/green] total={res.total_records} kept={res.kept_records} "
        f"dropped={res.total_records - res.kept_records} → {res.fasta_path}"
    )


# ── motifs ────────────────────────────────────────────────────────────────────

@app.command("motifs")
def cmd_motifs(
    fasta: Annotated[Path, typer.Argument()],
    out: Annotated[Path, typer.Option("-o")] = Path("data/motifs"),
    k: Annotated[int, typer.Option()] = 9,
    top: Annotated[int, typer.Option()] = 200,
    min_count: Annotated[int, typer.Option()] = 1,
    canonical: Annotated[bool, typer.Option("--canonical/--no-canonical")] = False,
) -> None:
    """Extrait les k-mers (motifs) des séquences."""
    from .motifs import MotifRequest, compute_motifs

    req = MotifRequest(input_fasta=fasta, out_dir=out, k=k, top=top, min_count=min_count, canonical=canonical)
    res = compute_motifs(req)
    console.print(
        f"[green]✓[/green] {res.total_records} records, {res.unique_kmers_after_filters} k-mers "
        f"→ {res.motifs_csv}"
    )


# ── cluster ───────────────────────────────────────────────────────────────────

@app.command("cluster")
def cmd_cluster(
    fasta: Annotated[Path, typer.Argument()],
    out: Annotated[Path, typer.Option("-o")] = Path("data/clusters"),
    k: Annotated[int, typer.Option()] = 9,
    min_sim: Annotated[float, typer.Option(help="Seuil Jaccard")] = 0.35,
    max_pairs: Annotated[int, typer.Option()] = 200_000,
    method: Annotated[str, typer.Option(help="jaccard (exact) | minhash (approx LSH, grands datasets)")] = "jaccard",
    minhash_hashes: Annotated[int, typer.Option(help="Taille signature MinHash")] = 128,
    minhash_bands: Annotated[int, typer.Option(help="Nombre de bandes LSH")] = 16,
    metadata: Annotated[Optional[Path], typer.Option()] = None,
) -> None:
    """Regroupe les séquences similaires — Jaccard exact ou MinHash LSH (grands datasets)."""
    from .clustering import ClusterRequest, cluster_sequences

    req = ClusterRequest(
        input_fasta=fasta, out_dir=out, k=k,
        min_similarity=min_sim, max_pairs=max_pairs, metadata_path=metadata,
        method=method, minhash_hashes=minhash_hashes, minhash_bands=minhash_bands,
    )
    res = cluster_sequences(req)
    console.print(
        f"[green]✓[/green] {res.total_records} records → {res.cluster_count} clusters "
        f"({res.edges_kept} edges) → {res.clusters_path}"
    )


# ── consensus ─────────────────────────────────────────────────────────────────

@app.command("consensus")
def cmd_consensus(
    fasta: Annotated[Path, typer.Argument(help="FASTA filtré")],
    clusters: Annotated[Path, typer.Option(help="clusters.jsonl")],
    out: Annotated[Path, typer.Option("-o")] = Path("data/consensus"),
    min_cluster_size: Annotated[int, typer.Option()] = 2,
    method: Annotated[str, typer.Option(help="kmer | simple")] = "kmer",
    metadata: Annotated[Optional[Path], typer.Option()] = None,
) -> None:
    """Construit les séquences consensus par cluster."""
    from .consensus import ConsensusRequest, build_consensus

    req = ConsensusRequest(
        clusters_path=clusters, fasta_path=fasta, out_dir=out,
        min_cluster_size=min_cluster_size, method=method, metadata_path=metadata,
    )
    manifest = build_consensus(req)
    console.print(
        f"[green]✓[/green] ok={manifest['nb_consensus_ok']} "
        f"échoués={manifest['nb_consensus_failed']} → {out / 'consensus.fasta'}"
    )


# ── blast ─────────────────────────────────────────────────────────────────────

@app.command("blast")
def cmd_blast(
    fasta: Annotated[Path, typer.Argument()],
    out: Annotated[Path, typer.Option("-o")] = Path("data/blast"),
    program: Annotated[str, typer.Option()] = "blastn",
    db: Annotated[str, typer.Option()] = "nt",
    max_hits: Annotated[int, typer.Option()] = 10,
    timeout: Annotated[int, typer.Option(help="Timeout BLAST (s)")] = 600,
) -> None:
    """Soumet des séquences à NCBI BLAST (polling asynchrone, tenacity retry)."""
    from .blast import BlastRequest, run_blast

    try:
        from .config import ncbi
        email, api_key, tool = ncbi.email, ncbi.api_key, ncbi.tool
    except Exception:
        email = typer.prompt("Email NCBI (requis pour BLAST)")
        api_key, tool = None, "arnfinder_v3"

    req = BlastRequest(
        input_fasta=fasta, out_dir=out,
        program=program, db=db, max_hits=max_hits, timeout=timeout,
        email=email, api_key=api_key, tool=tool,
    )
    res = run_blast(req)
    console.print(
        f"[green]✓[/green] {res.total_queries} queries, {res.total_hits} hits "
        f"en {res.elapsed:.1f}s → {res.hits_jsonl}"
    )


# ── export ────────────────────────────────────────────────────────────────────

@app.command("export")
def cmd_export(
    features: Annotated[Path, typer.Argument(help="CSV features (filter_report ou motifs_by_record)")],
    metadata: Annotated[Path, typer.Argument(help="JSONL metadata NCBI")],
    out: Annotated[Path, typer.Option("-o")] = Path("data/export"),
    motifs: Annotated[Optional[Path], typer.Option(help="motifs_by_record.csv")] = None,
    validate: Annotated[bool, typer.Option("--validate/--no-validate", help="Validation Pandera")] = True,
) -> None:
    """Exporte les données au format Parquet ML-ready (avec validation Pandera optionnelle)."""
    from .export import export_ml_ready

    paths = export_ml_ready(
        features_csv=features, metadata_jsonl=metadata,
        out_dir=out, motifs_csv=motifs, validate=validate,
    )
    console.print(f"[green]✓[/green] Export → {out}")
    console.print(f"  features  : {paths.features_parquet}")
    console.print(f"  joined    : {paths.joined_parquet}")
    if paths.motifs_parquet:
        console.print(f"  motifs    : {paths.motifs_parquet}")
    console.print(f"  manifest  : {paths.manifest_json}")


# ── structure ─────────────────────────────────────────────────────────────────

@app.command("structure")
def cmd_structure(
    fasta: Annotated[Path, typer.Argument(help="FASTA en entrée (DNA ou RNA)")],
    out: Annotated[Path, typer.Option("-o")] = Path("data/structure"),
    max_len: Annotated[int, typer.Option(help="Longueur max (séquences plus longues ignorées)")] = 2000,
    temperature: Annotated[float, typer.Option(help="Température de repliement (°C)")] = 37.0,
    backend: Annotated[Optional[str], typer.Option(help="viennarna | rnafold_subprocess | seqfold | nussinov_approx")] = None,
) -> None:
    """Prédit la structure secondaire ARN (ViennaRNA > RNAfold > Nussinov approximatif)."""
    from .secondary_structure import StructureRequest, compute_secondary_structure

    req = StructureRequest(
        input_fasta=fasta, out_dir=out,
        max_len=max_len, temperature=temperature,
        backend=backend,  # type: ignore[arg-type]
    )
    res = compute_secondary_structure(req)
    console.print(
        f"[green]✓[/green] {res.computed} structures calculées, "
        f"{res.skipped} ignorées, backend=[bold]{res.backend_used}[/bold] ({res.elapsed:.1f}s)"
    )
    console.print(f"  CSV  : {res.structures_csv}")
    console.print(f"  JSONL: {res.structures_jsonl}")


# ── db ────────────────────────────────────────────────────────────────────────

@app.command("db")
def cmd_db(
    data_dir: Annotated[Path, typer.Argument(help="Répertoire de données du run")],
    db_path: Annotated[Path, typer.Option(help="Chemin du fichier DuckDB")] = Path("data/arnfinder.duckdb"),
    query: Annotated[Optional[str], typer.Option("-q", help="Requête SQL à exécuter après import")] = None,
    export_table: Annotated[Optional[str], typer.Option(help="Table à exporter en Parquet")] = None,
    export_path: Annotated[Optional[Path], typer.Option(help="Chemin Parquet de sortie")] = None,
) -> None:
    """Importe toutes les sorties pipeline dans DuckDB et permet des requêtes SQL."""
    from .db import ARNFinderDB
    from rich.table import Table as RichTable

    with ARNFinderDB(db_path) as db:
        counts = db.import_all(data_dir)

        # Affiche les compteurs d'import
        t = RichTable(title=f"DuckDB : {db_path}", show_lines=True)
        t.add_column("Table", style="cyan")
        t.add_column("Lignes importées", justify="right")
        for table, count in counts.items():
            t.add_row(table, str(count))
        console.print(t)

        # Requête SQL optionnelle
        if query:
            console.print(f"\n[bold]SQL :[/bold] {query}")
            df = db.query(query)
            console.print(df.to_string(index=False))

        # Export Parquet optionnel
        if export_table and export_path:
            db.export_parquet(export_table, export_path)
            console.print(f"[green]✓[/green] Table [bold]{export_table}[/bold] exportée → {export_path}")


# ── report ────────────────────────────────────────────────────────────────────

@app.command("report")
def cmd_report(
    data_dir: Annotated[Path, typer.Argument(help="Répertoire de données du run")],
    query_str: Annotated[str, typer.Option("--query", "-q", help="Requête NCBI d'origine")] = "",
    out: Annotated[Optional[Path], typer.Option("-o", help="Chemin du HTML généré")] = None,
    top_kmers: Annotated[int, typer.Option()] = 20,
    top_clusters: Annotated[int, typer.Option()] = 15,
) -> None:
    """Génère un rapport HTML complet du run (dark mode, responsive)."""
    from .report import generate_report

    report_path = generate_report(
        data_dir=data_dir,
        query=query_str,
        out_path=out,
        top_kmers=top_kmers,
        top_clusters=top_clusters,
    )
    console.print(f"[green]✓[/green] Rapport généré → [bold]{report_path}[/bold]")


# ── codon_usage ───────────────────────────────────────────────────────────────

@app.command("codons")
def cmd_codons(
    fasta: Annotated[Path, typer.Argument(help="FASTA en entrée (DNA ou RNA)")],
    out: Annotated[Path, typer.Option("-o")] = Path("data/codon_usage"),
    reading_frame: Annotated[int, typer.Option(help="Cadre de lecture : 0, 1 ou 2")] = 0,
    min_codons: Annotated[int, typer.Option(help="Nombre min de codons pour traiter la séquence")] = 10,
) -> None:
    """Calcule les fréquences de codons, RSCU, GC1/2/3 et ENc par séquence."""
    from .codon_usage import CodonUsageRequest, compute_codon_usage

    req = CodonUsageRequest(
        input_fasta=fasta, out_dir=out,
        reading_frame=reading_frame, min_codons=min_codons,
    )
    res = compute_codon_usage(req)
    console.print(
        f"[green]✓[/green] {res.computed}/{res.total_records} séquences calculées "
        f"({res.skipped} ignorées) en {res.elapsed:.1f}s"
    )
    console.print(f"  Stats : {res.stats_csv}")
    console.print(f"  RSCU  : {res.rscu_csv}")


# ── embed ─────────────────────────────────────────────────────────────────────

@app.command("embed")
def cmd_embed(
    fasta: Annotated[Path, typer.Argument(help="FASTA en entrée")],
    out: Annotated[Path, typer.Option("-o")] = Path("data/embeddings"),
    backend: Annotated[str, typer.Option(
        help="Backend : auto | kmer_tfidf | nucleotide_transformer | esm2"
    )] = "auto",
    k: Annotated[int, typer.Option(help="Taille k-mers (TF-IDF uniquement)")] = 6,
    vocab_size: Annotated[int, typer.Option(help="Taille vocabulaire TF-IDF")] = 4096,
    batch_size: Annotated[int, typer.Option(help="Batch size (transformers)")] = 8,
    max_length: Annotated[int, typer.Option(help="Longueur max tokens (transformers)")] = 512,
    device: Annotated[str, typer.Option(help="Dispositif : cpu | cuda")] = "cpu",
    tsv: Annotated[bool, typer.Option("--tsv/--no-tsv", help="Export TSV en plus du Parquet")] = False,
) -> None:
    """Calcule les embeddings de séquences (k-mer TF-IDF ou transformers biologiques)."""
    from .embeddings import EmbeddingRequest, compute_embeddings

    req = EmbeddingRequest(
        input_fasta=fasta,
        out_dir=out,
        backend=backend,
        k=k,
        tfidf_vocab_size=vocab_size,
        batch_size=batch_size,
        max_length=max_length,
        device=device,
        save_tsv=tsv,
    )
    res = compute_embeddings(req)
    console.print(
        f"[green]✓[/green] {res.total_records} séquences, "
        f"dim={res.embedding_dim}, backend=[bold]{res.backend_used}[/bold] "
        f"({res.elapsed_s:.1f}s)"
    )
    console.print(f"  Parquet : {res.embeddings_path}")
    console.print(f"  Manifest: {res.manifest_path}")
    if res.tsv_path:
        console.print(f"  TSV     : {res.tsv_path}")


# ── api serve ─────────────────────────────────────────────────────────────────

@app.command("serve")
def cmd_serve(
    host: Annotated[str, typer.Option(help="Adresse d'écoute")] = "127.0.0.1",
    port: Annotated[int, typer.Option(help="Port")] = 8000,
    reload: Annotated[bool, typer.Option("--reload/--no-reload", help="Hot-reload (dev)")] = False,
) -> None:
    """Lance le serveur FastAPI (arnfinderv3 serve)."""
    try:
        import uvicorn
    except ImportError:
        console.print("[red]uvicorn non installé. Lance : pip install uvicorn[standard][/red]")
        raise typer.Exit(1)

    console.print(f"[bold cyan]ARN Finder V3 API[/bold cyan] → http://{host}:{port}/docs")
    uvicorn.run(
        "arn_finder_v3.api:create_app",
        factory=True,
        host=host, port=port, reload=reload,
    )
