# ARN Finder V3 — Plan d'implémentation

## Statut V3 : COMPLET ✓

Tous les modules sont implémentés et câblés.

---

## Ce qui change par rapport à V2

### 1. Orchestrateur de pipeline ✓ (killer-feature)
**Problème V1/V2** : 6 commandes manuelles à lancer dans le bon ordre.
**V3** : `arnfinderv3 run "..."` enchaîne tout automatiquement avec tableau de résultats Rich.
- Fichier : `pipeline/__init__.py`

### 2. Fetch asynchrone aiohttp ✓
**Problème V1/V2** : `requests` séquentiel, sleep 0.34s entre batches → lent.
**V3** : `aiohttp` + `asyncio.gather` → batches FASTA en parallèle. Gain estimé **3-5x**.
- Fichier : `fetch/__init__.py`

### 3. Retry tenacity sur tous les appels NCBI ✓
**Problème V1/V2** : crash immédiat sur erreur réseau temporaire.
**V3** : 5 tentatives, backoff exponentiel 2s→30s sur fetch, BLAST submit/poll/fetch.
- Fichiers : `fetch/__init__.py`, `blast/__init__.py`

### 4. Config Pydantic Settings ✓
**Problème V1/V2** : dataclasses sans validation, erreurs silencieuses.
**V3** : `pydantic-settings` valide chaque paramètre au démarrage.
- Fichier : `config.py`

### 5. CLI Typer + Rich ✓
**Problème V1/V2** : argparse verbeux, pas d'autocomplétion.
**V3** : `typer` génère autocomplétion shell native + help coloré Rich.
- Fichier : `cli.py`

### 6. Filtre GC ✓ (nouveau)
**Absent en V1/V2**. V3 ajoute `--min-gc` et `--max-gc` dans le filtre qualité.
- Fichier : `filters/__init__.py`

### 7. BLAST asynchrone ✓
**Problème V1/V2** : polling BLAST synchrone bloquant.
**V3** : polling asyncio avec `aiohttp`, retry tenacity, cache résultats raw/ par séquence.
- Fichier : `blast/__init__.py`

### 8. Validation Pandera sur l'export ✓
**Problème V1/V2** : export sans validation, erreurs silencieuses.
**V3** : `pandera` valide schémas avant écriture Parquet. Mode `lazy=True` = rapport complet.
- Fichier : `export/__init__.py`

### 9. Tests complets ✓
- `test_filters.py` — 6 tests (filtre, GC, dédup, métriques)
- `test_motifs.py` — 5 tests (k-mers, canonical, ignore_N)
- `test_clustering.py` — 6 tests (jaccard, union-find, outputs)
- `test_consensus.py` — 6 tests (merge, identité, seuils)
- `test_blast_parser.py` — 4 tests (parse JSON2, max_hits, fallback)
- `test_export.py` — 3 tests (Parquet, Pandera, manifest)
- `test_pipeline.py` — 2 tests intégration (mock fetch)
- `test_fetch.py` — squelette (respx à compléter)

---

## Structure finale V3

```
arn_finderV3/
├── pyproject.toml
├── .env.example
├── PLAN_V3.md
├── tests/
│   ├── test_filters.py       ✓
│   ├── test_motifs.py        ✓
│   ├── test_clustering.py    ✓
│   ├── test_consensus.py     ✓
│   ├── test_blast_parser.py  ✓
│   ├── test_export.py        ✓
│   ├── test_pipeline.py      ✓
│   └── test_fetch.py         (squelette respx)
└── src/arn_finder_v3/
    ├── __init__.py
    ├── cli.py                ✓  Typer, toutes commandes câblées
    ├── config.py             ✓  Pydantic Settings
    ├── io/__init__.py        ✓  RichHandler, helpers FASTA/CSV/JSON
    ├── fetch/__init__.py     ✓  aiohttp + tenacity
    ├── filters/__init__.py   ✓  + filtre GC, progress Rich
    ├── motifs/__init__.py    ✓  progress Rich
    ├── clustering/
    │   ├── __init__.py       ✓  progress Rich
    │   └── union_find.py     ✓
    ├── consensus/
    │   ├── __init__.py       ✓  progress Rich
    │   └── merge.py          ✓
    ├── blast/
    │   ├── __init__.py       ✓  async polling tenacity
    │   └── parser.py         ✓
    ├── export/__init__.py    ✓  Pandera + Parquet
    └── pipeline/__init__.py  ✓  orchestrateur complet
```

---

## V3.1 — Feuille de route (après stabilisation V3)

### A. Structure secondaire ARN
- Intégrer ViennaRNA via subprocess ou API RNAfold
- Ajouter colonne `mfe` (minimum free energy) dans l'export Parquet
- Nouveau sous-module : `arn_finder_v3/secondary_structure/`

### B. Rapport HTML par run
- `jinja2` + template Bootstrap minimaliste
- Génère `report.html` dans `data/` à la fin du pipeline
- Contient : stats par étape, top k-mers, distribution longueurs, top clusters
- Nouveau sous-module : `arn_finder_v3/report/`

### C. DuckDB local
- Remplacer les CSV intermédiaires par une base DuckDB `arnfinder.duckdb`
- Permet `SELECT * FROM clusters WHERE size > 5` directement
- Backward-compatible : export CSV/Parquet toujours disponible
- Nouveau sous-module : `arn_finder_v3/db/`

### D. FastAPI
- Wrapper HTTP du pipeline pour intégration CI/CD ou Jupyter
- `POST /run` → lance pipeline async, retourne job_id
- `GET /status/{job_id}` → statut en temps réel (SSE)
- Nouveau module : `arn_finder_v3/api/`

### E. Couverture de tests
- Compléter `test_fetch.py` avec `respx` pour mocker aiohttp
- Objectif : **>80% coverage** mesuré avec `pytest-cov`
- Ajouter `pytest.ini` ou `[tool.pytest]` dans `pyproject.toml`
