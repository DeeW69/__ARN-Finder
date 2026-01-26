![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

# ARN Finder 
**Open-source Python pipeline turning public DNA/RNA sequences into exploratory, ML-ready artifacts.**

## Project overview
Public repositories like GenBank are rich yet messy. ARN Finder acts as a bridge between bioinformatics and data engineering: it downloads curated subsets via NCBI E-utilities, enforces quality filters, deduplicates, extracts motifs, clusters similar fragments, and produces exploratory consensus sequences. Guiding principles: transparency, reproducibility, responsible use of public data.

## Key capabilities
- Targeted GenBank/RefSeq queries through Entrez E-utilities.
- Throttled downloads with local cache and manifest logging.
- Quality filters (length, N%, GC%), dedupe, per-organism caps.
- Global and per-record k-mer extraction.
- Similarity clustering using k-mer Jaccard sets.
- Experimental consensus generation per cluster.
- End-to-end CLI that can be scripted or integrated in ML workflows.

## Installation (Windows, Python 3.10+)
`powershell
python -m venv .venv
.\.venv\Scripts\activate
pip install -e .
`

## Pipeline overview
`
fetch ? filter ? motifs ? cluster ? consensus
`
Each step emits curated FASTA/CSV/JSON outputs plus a manifest capturing parameters and counts.

## Quick demo & Fetch & filter
`powershell
arnfinder fetch --query "Aves[Organism] AND COI[Gene] AND mitochondrion[Filter]" --limit 200 --out-dir data/raw --cache-dir data/cache --email you@example.com
arnfinder filter --in-fasta data/raw/sequences.fasta --in-metadata data/raw/metadata.jsonl --out-dir data/filtered --min-len 200 --max-n-frac 0.05 --alphabet AUTO --dedupe
`
Outputs: raw FASTA/JSONL/IDs, then filtered FASTA + filter_report + metadata subset.

## Motif analysis (k-mers)
Useful for feature engineering or motif discovery.
`powershell
arnfinder motifs --in-fasta data/filtered/filtered_sequences.fasta --out-dir data/motifs --k 9 --top 200 --min-count 2 --ignore-N --per-record
`
Produces: motifs.csv, motifs_by_record.csv, motifs_summary.json.

## Similarity clustering
Groups related fragments (e.g., mitochondrial COI in birds). k-mer sets ? Jaccard similarity ? union-find components. Includes cluster edges, per-cluster stats, manifest.

## Consensus (experimental)
Not a biological truth: a heuristic representative per cluster.
- Overlap-based merging with identity checks (fallback to longest sequence).
- Flagged low_quality if consensus > max_n_frac.
- Outputs: consensus.fasta, consensus_stats.csv, consensus_manifest.json.

## Validation (V2.1) : BLAST optionnel
Exploratory validation of consensus or filtered sequences against the NCBI BLAST service (Common URL API). Each query is throttled, polled until ready, and produces three artifacts:
- `blast_hits.jsonl`: 1 JSON object per sequence with RID, cached raw file path, and the parsed top hits (accession, title, organism, pident, coverage, e-value, bitscore).
- `blast_summary.csv`: compact table of the best hit per sequence for rapid inspection.
- `blast_manifest.json`: parameters (program, db, format, rate-limit, poll cadence), timestamps, totals (ok/failed) and cache hints.

### Configuration
Set the following environment variables (copy `.env.example` to `.env`, automatically loaded on CLI start):
```
NCBI_EMAIL=you@example.com          # required by NCBI
NCBI_TOOL=arnfinder                 # reported tool name, default: arnfinder
NCBI_API_KEY=your_api_key_optional  # optional but unlocks higher rate limits
BLAST_PROGRAM=blastn                # default program (env override)
BLAST_DB=nt                         # default database (env override)
```

Respect NCBI guidelines: identify yourself with email/tool, keep at most ~1 request/second without API key, enable caching via `--cache-dir` to avoid resubmitting identical sequences, and treat BLAST as on-demand validation (no mirroring).

### Usage examples
```powershell
# Validating consensus sequences (falls back to filtered_sequences.fasta if consensus is missing)
arnfinder blast --out-dir data/blast --cache-dir data/blast/cache --format JSON2

# Alternate database/format, smaller hit window, faster polling
arnfinder blast --program blastn --db nt --max-hits 5 --poll-seconds 5 --format Tabular --out-dir data/blast_tabular
```

BLAST outputs are exploratory sanity checks, not authoritative annotations. Always re-validate and interpret in context of the biological question.

## Scientific interpretation & limits
- Partial, biased public data; optional BLAST validation is exploratory, not definitive phylogenetics.
- Intended for exploration, pedagogy, ML prototyping; revalidation is mandatory for any biological claim.

## NCBI compliance
- Query-on-demand, no mirroring.
- Respect for rate limits; configurable sleep / API key.
- Local caching to minimize repeated requests.
- Explicit citation of NCBI as data source.

## Roadmap (V2)
- Enrich BLAST alignment summaries with more taxonomic context.
- Stronger clustering (MinHash / graph embeddings).
- Richer RNA-specific handling.
- Advanced genomic descriptors (codon usage, signatures).
- Enhanced ML-ready exports (Parquet, embeddings, metadata joins).

## Target audiences
- AI/bioinformatics internship proposals.
- Academic research groups.
- Medtech/biotech startups.
- ML teams prototyping on public omics data.


## License & citation
- License: MIT.
- Data source: NCBI GenBank / RefSeq (cite according to NCBI guidelines).
