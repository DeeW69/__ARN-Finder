![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Python](https://img.shields.io/badge/Python-3.11%2B-blue.svg)
![Version](https://img.shields.io/badge/version-3.2.0-brightgreen.svg)

# ARN Finder

**Open-source Python pipeline turning public DNA/RNA sequences into exploratory, ML-ready artifacts — with RNA secondary structure prediction, a REST API, a local DuckDB database, and an HTML report.**

---

## Table of contents

1. [Project overview](#project-overview)
2. [Changelog](#changelog)
3. [Architecture](#architecture)
4. [Installation](#installation)
5. [Quick start — one command](#quick-start--one-command)
6. [CLI reference](#cli-reference)
7. [REST API](#rest-api)
8. [DuckDB local database](#duckdb-local-database)
9. [Secondary structure prediction](#secondary-structure-prediction)
10. [HTML report](#html-report)
11. [Configuration](#configuration)
12. [Pipeline outputs](#pipeline-outputs)
13. [Scientific limits](#scientific-limits)
14. [NCBI compliance](#ncbi-compliance)
15. [Roadmap](#roadmap)
16. [License](#license)

---

## Project overview

Public repositories like GenBank are rich yet messy. ARN Finder bridges bioinformatics and data engineering:

- Downloads curated sequence subsets via NCBI E-utilities (async, with retry)
- Enforces quality filters (length, N%, GC%), deduplicates, enforces per-organism caps
- Extracts k-mer motifs and clusters similar fragments via Jaccard similarity
- Generates exploratory consensus sequences per cluster
- **Predicts RNA secondary structure** (dot-bracket + MFE) via ViennaRNA, RNAfold, or a built-in Nussinov fallback
- Exports everything as **Parquet ML-ready datasets** (validated with Pandera)
- Indexes all outputs in a **local DuckDB database** for SQL queries
- Generates a **dark-mode HTML report** at the end of every run
- Exposes a **FastAPI REST interface** with async job management and Server-Sent Events

Guiding principles: transparency, reproducibility, responsible use of public data.

---

## Changelog

### V3.2 — 2026-05-21 (current)
- **Codon usage** — fréquences de codons, RSCU (Sharp & Li 1986), GC1/GC2/GC3, ENc (Wright 1990) par séquence. Nouveau module `codon_usage/`, commande `arnfinderv3 codons`, tables DuckDB `codon_usage` + `codon_rscu`, section rapport HTML
- **Persistent job storage** — le registre de jobs FastAPI est maintenant backed par SQLite (`~/.arnfinder/jobs.db`). Les jobs survivent aux redémarrages de l'API
- **BLAST taxonomy enrichment** — les hits BLAST sont enrichis avec la lignée taxonomique complète (kingdom, phylum, class, order, family, genus) via NCBI efetch. Le CSV blast_summary inclut ces colonnes
- Pipeline étendu : `fetch → filter → motifs → cluster → consensus → structure → codon_usage → [blast] → export → report`
- Nouveau CLI : `arnfinderv3 codons`
- 46 fichiers, 14 modules de tests

### V3.1 — 2026-05-21
- **Secondary structure ARN** — dot-bracket + MFE via ViennaRNA bindings, RNAfold subprocess, or built-in Nussinov (always available, no extra install)
- **FastAPI** — async REST API (`arnfinderv3 serve`), SSE stream, job cancel, `/health` endpoint
- **DuckDB** — local SQL database (`arnfinderv3 db`), 6 tables, pre-built analytical queries
- **HTML report** — dark-mode responsive report auto-generated at end of every run
- Pipeline extended: `fetch → filter → motifs → cluster → consensus → structure → [blast] → export → report`
- New CLI commands: `structure`, `db`, `report`, `serve`
- 39 files, 12 test modules

### V3.0 — 2026-05-20
- **Async fetch** with `aiohttp` + `asyncio.gather` (3–5× faster on large datasets)
- **Retry with exponential backoff** via `tenacity` on all NCBI calls
- **Pydantic Settings** replace dataclasses — validation at startup, explicit errors
- **Typer CLI** replaces argparse — shell autocompletion, colored help
- **Pipeline orchestrator** (`arnfinderv3 run "..."`) — one command instead of six
- **GC filter** (`--min-gc` / `--max-gc`) added to the filter stage
- **Pandera schema validation** on Parquet exports
- **Rich** progress bars and colored logging throughout
- **Async BLAST polling** — non-blocking, tenacity retry
- `src/` layout, Python 3.11+ required
- 27 files, 8 test modules

### V2 — 2024
- ML-ready Parquet export (pandas + pyarrow)
- Extended test suite
- French-language documentation and logging

### V1 — 2024
- Complete pipeline: fetch → filter → motifs → cluster → consensus → BLAST
- Synchronous, argparse CLI, dataclasses config

---

## Architecture

```
arn_finderV3/
├── pyproject.toml
├── .env.example
└── src/arn_finder_v3/
    ├── cli.py                   # Typer — 10 commands
    ├── config.py                # Pydantic Settings
    ├── io/                      # FASTA/CSV/JSON helpers, RichHandler
    ├── fetch/                   # aiohttp async + tenacity retry
    ├── filters/                 # Quality filter (length, N%, GC%)
    ├── motifs/                  # K-mer extraction (canonical form)
    ├── clustering/              # Jaccard similarity + union-find
    ├── consensus/               # Overlap-based cluster consensus
    ├── secondary_structure/     # Dot-bracket + MFE (ViennaRNA / RNAfold / Nussinov)
    ├── codon_usage/             # RSCU, GC1/2/3, ENc (V3.2)
    ├── blast/                   # Async BLAST polling + parser + taxonomy enrichment
    ├── export/                  # Parquet + Pandera schema validation
    ├── db/                      # DuckDB — import, SQL queries, export
    ├── report/                  # HTML report (Jinja2, dark mode)
    ├── api/                     # FastAPI — jobs, SSE, cancel
    └── pipeline/                # Orchestrator — full run in one call
```

**Pipeline flow:**
```
fetch → filter → motifs → cluster → consensus → structure → codon_usage → [blast] → export → report
```

Each stage emits curated FASTA/CSV/JSON outputs plus a manifest capturing exact parameters, counts, and timestamps.

---

## Installation

### Requirements
- Python 3.11+
- (optional) ViennaRNA for accurate secondary structure: `conda install -c bioconda viennarna`

```powershell
# Clone
git clone https://github.com/DeeW69/__ARN-Finder.git
cd __ARN-Finder/arn_finderV3

# Virtual environment
python -m venv .venv
.\.venv\Scripts\activate        # Windows
# source .venv/bin/activate    # Linux/macOS

# Install
pip install -e ".[dev]"

# Configure
copy .env.example .env
# Edit .env : set NCBI_EMAIL at minimum
```

---

## Quick start — one command

```powershell
arnfinderv3 run "BRCA1[Gene Name] AND Homo sapiens[Organism]" `
  --out data `
  --limit 100 `
  --min-gc 0.35 `
  --max-gc 0.75

# Output:
#   data/raw/           sequences.fasta, metadata.jsonl
#   data/filtered/      filtered_sequences.fasta, filter_report.csv
#   data/motifs/        motifs.csv, motifs_by_record.csv
#   data/clusters/      clusters.jsonl, cluster_stats.csv
#   data/consensus/     consensus.fasta, consensus_stats.csv
#   data/structure/     structures.csv, structures.jsonl   ← NEW V3.1
#   data/export/        features.parquet, joined_dataset.parquet
#   data/report.html                                        ← NEW V3.1
```

---

## CLI reference

```
arnfinderv3 --help

Commands:
  run        Full pipeline (fetch → ... → report)
  fetch      Download sequences from NCBI Entrez
  filter     Quality filter (length, N%, GC%, dedupe)
  motifs     K-mer extraction
  cluster    Jaccard similarity clustering
  consensus  Cluster consensus sequences
  structure  RNA secondary structure prediction       ← V3.1
  codons     Codon usage : RSCU, GC1/2/3, ENc        ← V3.2
  blast      NCBI BLAST (async polling + taxonomy)   ← taxonomy V3.2
  export     Parquet ML-ready export (Pandera)
  db         Import outputs into DuckDB, run SQL      ← V3.1
  report     Generate HTML run report                 ← V3.1
  serve      Start FastAPI REST server                ← V3.1
```

### Key options for `run`

| Option | Default | Description |
|--------|---------|-------------|
| `--limit` | 200 | Max sequences to fetch |
| `--min-len` | 200 | Minimum sequence length |
| `--max-n-frac` | 0.05 | Max fraction of N bases |
| `--min-gc` | 0.0 | Min GC fraction (0–1) |
| `--max-gc` | 1.0 | Max GC fraction (0–1) |
| `--k` | 9 | K-mer size (motifs + clustering) |
| `--min-sim` | 0.35 | Jaccard threshold for clustering |
| `--no-structure` | — | Disable secondary structure |
| `--blast` | — | Enable BLAST validation |
| `--no-export` | — | Disable Parquet export |
| `--no-report` | — | Disable HTML report |

---

## REST API

```powershell
arnfinderv3 serve           # http://127.0.0.1:8000/docs
arnfinderv3 serve --port 9000 --reload   # dev mode
```

### Endpoints

| Method | Path | Description |
|--------|------|-------------|
| `GET` | `/health` | API health + available backends |
| `POST` | `/pipeline/run` | Launch full pipeline (async, returns `job_id`) |
| `GET` | `/pipeline/{id}/status` | Job status + per-stage elapsed time |
| `GET` | `/pipeline/{id}/stream` | **Server-Sent Events** real-time progress |
| `GET` | `/pipeline/{id}/results` | Full results (completed jobs only) |
| `DELETE` | `/pipeline/{id}` | Cancel a running job |
| `GET` | `/pipeline` | List all jobs |
| `POST` | `/structure/run` | Standalone structure prediction |

#### Example — launch pipeline via API

```python
import httpx, json, time

r = httpx.post("http://localhost:8000/pipeline/run", json={
    "query": "COI[Gene] AND Aves[Organism]",
    "limit": 50,
    "run_blast": False,
})
job_id = r.json()["job_id"]

# Poll status
while True:
    s = httpx.get(f"http://localhost:8000/pipeline/{job_id}/status").json()
    print(s["status"], [st["stage"] for st in s["stages"] if st["status"] == "completed"])
    if s["status"] in ("completed", "failed"):
        break
    time.sleep(5)
```

---

## DuckDB local database

```powershell
# Import all pipeline outputs and run a SQL query
arnfinderv3 db data/ -q "SELECT organism, COUNT(*) FROM sequences GROUP BY organism ORDER BY 2 DESC"

# Export a table to Parquet
arnfinderv3 db data/ --export-table structure --export-path data/export/structure.parquet
```

#### Available tables

| Table | Content |
|-------|---------|
| `sequences` | All filtered sequences (quality metrics, organism) |
| `motifs` | K-mers per sequence |
| `clusters` | Cluster membership for each sequence |
| `consensus` | Consensus stats per cluster |
| `structure` | Dot-bracket + MFE per sequence |
| `blast_hits` | BLAST hit summaries |

#### Python usage

```python
from arn_finder_v3.db import ARNFinderDB
from pathlib import Path

with ARNFinderDB("data/arnfinder.duckdb") as db:
    db.import_all(Path("data/"))

    # Pre-built analytical queries
    print(db.top_kmers(n=10))
    print(db.cluster_stats())
    print(db.sequence_overview())

    # Free SQL
    df = db.query("""
        SELECT s.organism, AVG(st.mfe) AS mean_mfe, COUNT(*) AS n
        FROM sequences s
        JOIN structure st ON s.record_id = st.record_id
        WHERE s.kept = TRUE AND st.skipped = FALSE
        GROUP BY s.organism
        ORDER BY mean_mfe
    """)
    print(df)
```

---

## Secondary structure prediction

ARN Finder automatically picks the best available backend:

| Priority | Backend | Install | Accuracy |
|----------|---------|---------|----------|
| 1 | **ViennaRNA** Python bindings | `conda install -c bioconda viennarna` | High (thermodynamic) |
| 2 | **RNAfold** subprocess | Included with ViennaRNA system package | High (thermodynamic) |
| 3 | **Nussinov** (built-in) | None — always available | Approximate (maximizes base pairs) |

```powershell
# Predict on consensus sequences
arnfinderv3 structure data/consensus/consensus.fasta -o data/structure

# Force a specific backend
arnfinderv3 structure data/consensus/consensus.fasta --backend nussinov_approx

# Skip sequences longer than 1000 nt
arnfinderv3 structure data/consensus/consensus.fasta --max-len 1000
```

Output columns: `record_id`, `length`, `gc_frac`, `dot_bracket`, `mfe_kcal_mol`, `backend`, `skipped`, `skip_reason`.

> **Note:** DNA sequences are automatically converted to RNA (T→U) before prediction.
> MFE from the Nussinov fallback is an approximation (−1 kcal/mol per pair) and should not be used for thermodynamic comparisons.

---

## HTML report

An HTML report is automatically generated at the end of every `run`.

```powershell
# Generate standalone report from existing run data
arnfinderv3 report data/ --query "COI[Gene] AND Aves[Organism]"
# → data/report.html
```

The report includes:
- Global KPIs (sequences fetched/kept, clusters, structures, backend)
- Per-stage pipeline status and elapsed times
- Filter rejection reasons (bar chart)
- Top k-mers table
- Top clusters with dominant organism and mean MFE
- Secondary structure preview (dot-bracket + MFE per sequence)

Dark-mode, responsive, no external dependencies — single self-contained HTML file.

---

## Configuration

Copy `.env.example` to `.env` in the `arn_finderV3/` directory:

```env
# Required
NCBI_EMAIL=you@example.com

# Optional
NCBI_TOOL=arnfinder_v3
NCBI_API_KEY=your_api_key       # Unlocks 10 req/s instead of 3 req/s

# BLAST (optional)
BLAST_PROGRAM=blastn
BLAST_DB=nt
BLAST_MAX_HITS=10
BLAST_TIMEOUT=600
```

All settings are validated by Pydantic at startup — a clear error is raised if `NCBI_EMAIL` is missing.

---

## Pipeline outputs

| Stage | Key files | Format |
|-------|-----------|--------|
| **fetch** | `sequences.fasta`, `metadata.jsonl`, `ids.txt`, `manifest.json` | FASTA / JSONL / JSON |
| **filter** | `filtered_sequences.fasta`, `filter_report.csv`, `filter_manifest.json` | FASTA / CSV / JSON |
| **motifs** | `motifs.csv`, `motifs_by_record.csv`, `motifs_summary.json` | CSV / JSON |
| **cluster** | `clusters.jsonl`, `cluster_edges.csv`, `cluster_stats.csv`, `cluster_manifest.json` | JSONL / CSV / JSON |
| **consensus** | `consensus.fasta`, `consensus_stats.csv`, `consensus_manifest.json` | FASTA / CSV / JSON |
| **structure** *(V3.1)* | `structures.csv`, `structures.jsonl`, `structure_manifest.json` | CSV / JSONL / JSON |
| **blast** | `blast_hits.jsonl`, `blast_summary.csv`, `blast_manifest.json` | JSONL / CSV / JSON |
| **export** | `features.parquet`, `joined_dataset.parquet`, `motifs_by_record.parquet` | Parquet |
| **report** *(V3.1)* | `report.html` | HTML |

---

## Running tests

```powershell
cd arn_finderV3
pytest tests/ -v --tb=short

# With coverage
pytest tests/ --cov=arn_finder_v3 --cov-report=term-missing
```

Test suite: 12 modules, covers filters, motifs, clustering, consensus, BLAST parser, export (Parquet + Pandera), secondary structure (all backends), FastAPI (TestClient), DuckDB, HTML report, and pipeline integration.

---

## Scientific limits

- Data is partial and biased towards well-studied organisms in GenBank.
- Clustering (Jaccard k-mer) is sequence-similarity based, not phylogenetic.
- Consensus sequences are heuristic representatives, not biological ground truth.
- MFE from the Nussinov fallback is an approximation — use ViennaRNA for thermodynamic accuracy.
- BLAST validation is exploratory, not authoritative annotation.
- Revalidation is mandatory before any biological claim.
- Intended for exploration, pedagogy, and ML prototyping.

---

## NCBI compliance

- Query-on-demand only; no bulk mirroring.
- Rate limits respected: 3 req/s without API key, 10 req/s with key.
- Local caching minimizes repeated requests.
- All requests identify the tool and a contact email.
- Cite NCBI GenBank/RefSeq per their guidelines when publishing results.

---

## Roadmap

- [x] Persistent job storage (SQLite) for the FastAPI layer ← V3.2
- [x] Richer BLAST summaries with full taxonomic lineage ← V3.2
- [x] Codon usage and advanced genomic descriptors (RSCU, ENc) ← V3.2
- [ ] ViennaRNA Windows wheel bundled in optional dependency
- [ ] MinHash / graph embeddings for large-scale clustering
- [ ] Sequence embeddings export (ESM-2, Nucleotide Transformer)

---

## Target audiences

- AI/bioinformatics research and internship projects
- Academic research groups working on public genomics data
- Medtech/biotech startups prototyping on omics sequences
- ML teams exploring public biological datasets

---

## License & citation

- **License:** MIT
- **Data source:** NCBI GenBank / RefSeq — cite according to [NCBI guidelines](https://www.ncbi.nlm.nih.gov/home/about/policies/)
- **Secondary structure:** If using ViennaRNA, cite: Lorenz et al. (2011) *ViennaRNA Package 2.0*, Algorithms for Molecular Biology, 6:26
