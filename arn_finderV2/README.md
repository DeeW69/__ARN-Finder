# ARN Finder V2

Cette V2 prépare l'évolution d'ARN Finder sous la forme d'un package installable, avec une nouvelle CLI `arnfinderv2`.

## Installation locale

```bash
cd arn_finderV2
python -m venv .venv
.venv\Scripts\activate  # sous Windows
pip install --upgrade pip
pip install -e ".[dev]"
```

## Configuration

Copier `.env.example` vers `.env` (à la racine du dépôt) et remplir les valeurs :

```
NCBI_EMAIL=you@example.org
NCBI_TOOL=arnfinderv2
NCBI_API_KEY=changeme
BLAST_PROGRAM=blastn
BLAST_DB=nt
```

La CLI charge automatiquement ce fichier si présent.

## CLI

```bash
arnfinderv2 --help
arnfinderv2 blast --help
```

## BLAST (V2)

La commande `arnfinderv2 blast` soumet un FASTA à l'API BLAST (NCBI) avec une approche exploratoire (pas d'annotation experte).  
Le flux produit trois fichiers dans `data/blast/` par défaut :

- `blast_hits.jsonl` : une entrée par séquence avec ses hits principaux, chemin des résultats bruts.
- `blast_summary.csv` : table synthétique (accès, organisme, scores).
- `blast_manifest.json` : paramètres utilisés, timestamps, compteur succès/échecs, version CLI.

Respect NCBI :

- `--rate-limit-seconds` throttle toutes les requêtes (défaut 1.0s).
- `--poll-seconds` et `--timeout-seconds` contrôlent le polling (défaut 10s / 600s).
- `--cache-dir` (défaut `data/cache/blast`) conserve les RID déjà téléchargés pour éviter les re-fetch.

Exemple :

```bash
arnfinderv2 blast ^
  --in-fasta data/examples/sequences.fasta ^
  --program blastn ^
  --db nt ^
  --format JSON2 ^
  --max-hits 20
```

Notes :

- `--program` et `--db` retombent sur les variables `.env` (`BLAST_PROGRAM`, `BLAST_DB`) ou sur `blastn`/`nt`.
- Les séquences sont upper-case et `U` → `T` quand `blastn` est utilisé pour éviter les surprises.
- En cas d'échec réseau ou de parsing JSON2, le workflow continue et loggue un warning (top_hits = `[]` mais résultat brut conservé).

## Export ML-ready

Une fois les features filtrées et les metadata générées par la pipeline V1/V2, utilisez :

```bash
arnfinderv2 export ^
  --features-csv data/filtered/filtered_features.csv ^
  --metadata-jsonl data/raw/metadata.jsonl ^
  --motifs-by-record-csv data/motifs/motifs_by_record.csv ^
  --out-dir data/exports
```

Cette commande :

- charge les CSV/JSONL en vérifiant leur présence ;
- effectue une jointure `record_id`/`accession` quand possible ;
- exporte `features.parquet`, `motifs_by_record.parquet` (optionnel) et `joined_dataset.parquet` (features + metadata) ;
- écrit `export_manifest.json` avec version, paramètres et compteurs.

Les exports reposent sur `pandas` + `pyarrow` et sont pensés pour un Consommation ML (Spark, DuckDB, Polars…).
