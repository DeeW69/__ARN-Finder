"""
DuckDB local — V3.1

Remplace les CSV intermédiaires par une base DuckDB requêtable.
Le fichier `arnfinder.duckdb` est créé dans le répertoire `data/` du run.

Tables créées :
  sequences   — toutes les séquences filtrées (métadonnées qualité + NCBI)
  motifs      — k-mers par séquence
  clusters    — appartenance cluster par séquence
  consensus   — séquences consensus par cluster
  structure   — prédictions structure secondaire ARN
  blast_hits  — résultats BLAST

Usage :
  db = ARNFinderDB("data/arnfinder.duckdb")
  db.import_all("data/")
  df = db.query("SELECT * FROM sequences WHERE gc_frac > 0.5")
  db.export_parquet("sequences", "data/export/sequences.parquet")
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Optional

from ..io import get_logger, now_iso

logger = get_logger("arnfinderv3.db")

# ── Schema DDL ────────────────────────────────────────────────────────────────

DDL_SEQUENCES = """
CREATE TABLE IF NOT EXISTS sequences (
    record_id     VARCHAR PRIMARY KEY,
    accession     VARCHAR,
    organism      VARCHAR,
    length        INTEGER,
    gc_frac       DOUBLE,
    n_frac        DOUBLE,
    kept          BOOLEAN,
    reason        VARCHAR,
    imported_at   VARCHAR
);
"""

DDL_MOTIFS = """
CREATE TABLE IF NOT EXISTS motifs (
    record_id        VARCHAR,
    kmer             VARCHAR,
    count            INTEGER,
    freq_in_record   DOUBLE,
    imported_at      VARCHAR
);
"""

DDL_CLUSTERS = """
CREATE TABLE IF NOT EXISTS clusters (
    cluster_id         VARCHAR,
    record_id          VARCHAR,
    accession          VARCHAR,
    length             INTEGER,
    organism           VARCHAR,
    size               INTEGER,
    dominant_organism  VARCHAR,
    mean_pairwise_sim  DOUBLE,
    imported_at        VARCHAR
);
"""

DDL_CONSENSUS = """
CREATE TABLE IF NOT EXISTS consensus (
    cluster_id     VARCHAR PRIMARY KEY,
    cluster_size   INTEGER,
    consensus_len  INTEGER,
    n_frac         DOUBLE,
    gc_percent     DOUBLE,
    method         VARCHAR,
    status         VARCHAR,
    imported_at    VARCHAR
);
"""

DDL_STRUCTURE = """
CREATE TABLE IF NOT EXISTS structure (
    record_id     VARCHAR PRIMARY KEY,
    length        INTEGER,
    gc_frac       DOUBLE,
    dot_bracket   VARCHAR,
    mfe           DOUBLE,
    backend       VARCHAR,
    skipped       BOOLEAN,
    skip_reason   VARCHAR,
    imported_at   VARCHAR
);
"""

DDL_BLAST = """
CREATE TABLE IF NOT EXISTS blast_hits (
    query_id    VARCHAR,
    accession   VARCHAR,
    organism    VARCHAR,
    taxid       VARCHAR,
    pident      DOUBLE,
    evalue      DOUBLE,
    bitscore    DOUBLE,
    coverage    DOUBLE,
    kingdom     VARCHAR,
    phylum      VARCHAR,
    class       VARCHAR,
    "order"     VARCHAR,
    family      VARCHAR,
    genus       VARCHAR,
    imported_at VARCHAR
);
"""

DDL_CODON_USAGE = """
CREATE TABLE IF NOT EXISTS codon_usage (
    record_id     VARCHAR,
    length        INTEGER,
    gc_frac       DOUBLE,
    gc1           DOUBLE,
    gc2           DOUBLE,
    gc3           DOUBLE,
    enc           DOUBLE,
    total_codons  INTEGER,
    stop_codons   INTEGER,
    imported_at   VARCHAR
);
"""

DDL_CODON_RSCU = """
CREATE TABLE IF NOT EXISTS codon_rscu (
    record_id   VARCHAR,
    codon       VARCHAR,
    amino_acid  VARCHAR,
    count       INTEGER,
    freq        DOUBLE,
    rscu        DOUBLE,
    imported_at VARCHAR
);
"""

ALL_DDL = [
    DDL_SEQUENCES, DDL_MOTIFS, DDL_CLUSTERS,
    DDL_CONSENSUS, DDL_STRUCTURE, DDL_BLAST,
    DDL_CODON_USAGE, DDL_CODON_RSCU,
]


# ── Classe principale ─────────────────────────────────────────────────────────

class ARNFinderDB:
    """Wrapper autour d'une connexion DuckDB locale."""

    def __init__(self, db_path: Path | str) -> None:
        try:
            import duckdb
        except ImportError:
            raise ImportError(
                "duckdb non installé. Lance : pip install duckdb>=0.10"
            ) from None

        self._path = Path(db_path)
        self._path.parent.mkdir(parents=True, exist_ok=True)
        self._con = duckdb.connect(str(self._path))
        self._init_schema()
        logger.info("DuckDB ouvert : %s", self._path)

    # ── Initialisation ────────────────────────────────────────────────────────

    def _init_schema(self) -> None:
        for ddl in ALL_DDL:
            self._con.execute(ddl)

    # ── Import depuis les sorties pipeline ─────────────────────────────────────

    def import_all(self, data_dir: Path | str) -> dict[str, int]:
        """
        Importe toutes les sorties du pipeline depuis `data_dir`.
        Retourne le nombre de lignes importées par table.
        """
        data_dir = Path(data_dir)
        counts: dict[str, int] = {}
        ts = now_iso()

        counts["sequences"]   = self._import_sequences(data_dir, ts)
        counts["motifs"]      = self._import_motifs(data_dir, ts)
        counts["clusters"]    = self._import_clusters(data_dir, ts)
        counts["consensus"]   = self._import_consensus(data_dir, ts)
        counts["structure"]   = self._import_structure(data_dir, ts)
        counts["blast_hits"]  = self._import_blast(data_dir, ts)
        counts["codon_usage"] = self._import_codon_usage(data_dir, ts)
        counts["codon_rscu"]  = self._import_codon_rscu(data_dir, ts)

        logger.info("Import DuckDB terminé : %s", counts)
        return counts

    def _import_sequences(self, data_dir: Path, ts: str) -> int:
        report = data_dir / "filtered" / "filter_report.csv"
        if not report.exists():
            return 0
        self._con.execute("DELETE FROM sequences")
        self._con.execute(f"""
            INSERT INTO sequences
            SELECT
                record_id,
                accession,
                NULL AS organism,
                CAST(length AS INTEGER),
                CAST(gc_percent AS DOUBLE) / 100.0,
                CAST(n_frac AS DOUBLE),
                (kept = 'true') AS kept,
                reason,
                '{ts}' AS imported_at
            FROM read_csv_auto('{report}')
        """)
        return self._con.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]

    def _import_motifs(self, data_dir: Path, ts: str) -> int:
        csv = data_dir / "motifs" / "motifs_by_record.csv"
        if not csv.exists():
            return 0
        self._con.execute("DELETE FROM motifs")
        self._con.execute(f"""
            INSERT INTO motifs
            SELECT record_id, kmer,
                   CAST(count AS INTEGER),
                   CAST(freq_in_record AS DOUBLE),
                   '{ts}'
            FROM read_csv_auto('{csv}')
        """)
        return self._con.execute("SELECT COUNT(*) FROM motifs").fetchone()[0]

    def _import_clusters(self, data_dir: Path, ts: str) -> int:
        jsonl = data_dir / "clusters" / "clusters.jsonl"
        if not jsonl.exists():
            return 0
        self._con.execute("DELETE FROM clusters")
        rows_inserted = 0
        with jsonl.open(encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                cluster = json.loads(line)
                cid = cluster["cluster_id"]
                size = cluster["size"]
                dom = cluster.get("dominant_organism", "")
                mean_sim = cluster.get("mean_pairwise_sim", 0.0)
                for member in cluster.get("members", []):
                    self._con.execute(
                        "INSERT INTO clusters VALUES (?,?,?,?,?,?,?,?,?)",
                        [cid, member["record_id"], member.get("accession"),
                         member.get("length"), member.get("organism"),
                         size, dom, mean_sim, ts]
                    )
                    rows_inserted += 1
        return rows_inserted

    def _import_consensus(self, data_dir: Path, ts: str) -> int:
        csv = data_dir / "consensus" / "consensus_stats.csv"
        if not csv.exists():
            return 0
        self._con.execute("DELETE FROM consensus")
        self._con.execute(f"""
            INSERT INTO consensus
            SELECT
                cluster_id,
                CAST(cluster_size AS INTEGER),
                CAST(consensus_len AS INTEGER),
                CAST(n_frac AS DOUBLE),
                CAST(gc_percent AS DOUBLE),
                method,
                status,
                '{ts}'
            FROM read_csv_auto('{csv}')
        """)
        return self._con.execute("SELECT COUNT(*) FROM consensus").fetchone()[0]

    def _import_structure(self, data_dir: Path, ts: str) -> int:
        csv = data_dir / "structure" / "structures.csv"
        if not csv.exists():
            return 0
        self._con.execute("DELETE FROM structure")
        self._con.execute(f"""
            INSERT INTO structure
            SELECT
                record_id,
                CAST(length AS INTEGER),
                CAST(gc_frac AS DOUBLE),
                dot_bracket,
                CAST(mfe_kcal_mol AS DOUBLE),
                backend,
                (skipped = 'True'),
                skip_reason,
                '{ts}'
            FROM read_csv_auto('{csv}')
        """)
        return self._con.execute("SELECT COUNT(*) FROM structure").fetchone()[0]

    def _import_blast(self, data_dir: Path, ts: str) -> int:
        csv_path = data_dir / "blast" / "blast_summary.csv"
        if not csv_path.exists():
            return 0
        self._con.execute("DELETE FROM blast_hits")
        # Charge toutes les colonnes présentes — taxonomy colonnes optionnelles
        tmp = self._con.execute(
            f"SELECT * FROM read_csv_auto('{csv_path}') LIMIT 0"
        ).description
        cols = [d[0] for d in tmp] if tmp else []
        has_taxonomy = "kingdom" in cols

        if has_taxonomy:
            self._con.execute(f"""
                INSERT INTO blast_hits
                SELECT query_id, accession, organism,
                       COALESCE(taxid, ''),
                       CAST(pident AS DOUBLE),
                       CAST(evalue AS DOUBLE),
                       CAST(bitscore AS DOUBLE),
                       CAST(coverage AS DOUBLE),
                       COALESCE(kingdom,''), COALESCE(phylum,''),
                       COALESCE("class",''), COALESCE("order",''),
                       COALESCE(family,''), COALESCE(genus,''),
                       '{ts}'
                FROM read_csv_auto('{csv_path}')
            """)
        else:
            self._con.execute(f"""
                INSERT INTO blast_hits
                SELECT query_id, accession, organism,
                       '' AS taxid,
                       CAST(pident AS DOUBLE),
                       CAST(evalue AS DOUBLE),
                       CAST(bitscore AS DOUBLE),
                       CAST(coverage AS DOUBLE),
                       '','','','','','',
                       '{ts}'
                FROM read_csv_auto('{csv_path}')
            """)
        return self._con.execute("SELECT COUNT(*) FROM blast_hits").fetchone()[0]

    def _import_codon_usage(self, data_dir: Path, ts: str) -> int:
        csv_path = data_dir / "codon_usage" / "codon_stats.csv"
        if not csv_path.exists():
            return 0
        self._con.execute("DELETE FROM codon_usage")
        self._con.execute(f"""
            INSERT INTO codon_usage
            SELECT
                record_id,
                CAST(length AS INTEGER),
                CAST(gc_frac AS DOUBLE),
                CAST(gc1 AS DOUBLE),
                CAST(gc2 AS DOUBLE),
                CAST(gc3 AS DOUBLE),
                CAST(enc AS DOUBLE),
                CAST(total_codons AS INTEGER),
                CAST(stop_codons AS INTEGER),
                '{ts}'
            FROM read_csv_auto('{csv_path}')
        """)
        return self._con.execute("SELECT COUNT(*) FROM codon_usage").fetchone()[0]

    def _import_codon_rscu(self, data_dir: Path, ts: str) -> int:
        csv_path = data_dir / "codon_usage" / "codon_rscu.csv"
        if not csv_path.exists():
            return 0
        self._con.execute("DELETE FROM codon_rscu")
        self._con.execute(f"""
            INSERT INTO codon_rscu
            SELECT
                record_id, codon, amino_acid,
                CAST(count AS INTEGER),
                CAST(freq AS DOUBLE),
                CAST(rscu AS DOUBLE),
                '{ts}'
            FROM read_csv_auto('{csv_path}')
        """)
        return self._con.execute("SELECT COUNT(*) FROM codon_rscu").fetchone()[0]

    # ── Requêtes ──────────────────────────────────────────────────────────────

    def query(self, sql: str) -> Any:
        """Exécute une requête SQL et retourne un DataFrame pandas."""
        return self._con.execute(sql).df()

    def tables(self) -> list[str]:
        """Liste les tables disponibles."""
        result = self._con.execute(
            "SELECT table_name FROM information_schema.tables WHERE table_schema = 'main'"
        ).fetchall()
        return [r[0] for r in result]

    def count(self, table: str) -> int:
        return self._con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]

    # ── Export ────────────────────────────────────────────────────────────────

    def export_parquet(self, table: str, out_path: Path | str) -> Path:
        """Exporte une table en Parquet."""
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        self._con.execute(f"COPY {table} TO '{out_path}' (FORMAT PARQUET)")
        logger.info("Table %s exportée → %s", table, out_path)
        return out_path

    def export_csv(self, table: str, out_path: Path | str) -> Path:
        """Exporte une table en CSV."""
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        self._con.execute(f"COPY {table} TO '{out_path}' (HEADER, DELIMITER ',')")
        return out_path

    # ── Requêtes prédéfinies utiles ───────────────────────────────────────────

    def top_kmers(self, n: int = 20) -> Any:
        """Top N k-mers les plus fréquents toutes séquences confondues."""
        return self.query(f"""
            SELECT kmer, SUM(count) AS total_count, COUNT(*) AS n_sequences
            FROM motifs
            GROUP BY kmer
            ORDER BY total_count DESC
            LIMIT {n}
        """)

    def cluster_stats(self) -> Any:
        """Statistiques par cluster (taille, GC moyen, MFE moyen)."""
        return self.query("""
            SELECT
                c.cluster_id,
                c.size,
                c.dominant_organism,
                AVG(s.gc_frac)  AS mean_gc,
                AVG(st.mfe)     AS mean_mfe,
                COUNT(st.record_id) AS n_with_structure
            FROM clusters c
            LEFT JOIN sequences s  ON c.record_id = s.record_id
            LEFT JOIN structure st ON c.record_id = st.record_id
            GROUP BY c.cluster_id, c.size, c.dominant_organism
            ORDER BY c.size DESC
        """)

    def sequence_overview(self) -> Any:
        """Vue jointe séquences + structure + cluster."""
        return self.query("""
            SELECT
                s.record_id,
                s.organism,
                s.length,
                s.gc_frac,
                s.n_frac,
                st.mfe,
                st.dot_bracket,
                st.backend,
                c.cluster_id
            FROM sequences s
            LEFT JOIN structure st ON s.record_id = st.record_id
            LEFT JOIN clusters  c  ON s.record_id = c.record_id
            WHERE s.kept = TRUE
            ORDER BY s.length DESC
        """)

    def codon_bias_overview(self) -> Any:
        """Vue du biais codon : ENc moyen, GC3 moyen, par organisme si disponible."""
        return self.query("""
            SELECT
                s.organism,
                AVG(cu.enc)  AS mean_enc,
                AVG(cu.gc3)  AS mean_gc3,
                AVG(cu.gc1)  AS mean_gc1,
                COUNT(*)     AS n_sequences
            FROM codon_usage cu
            LEFT JOIN sequences s ON cu.record_id = s.record_id
            GROUP BY s.organism
            ORDER BY mean_enc ASC
        """)

    def top_preferred_codons(self, n: int = 20) -> Any:
        """Top N codons avec RSCU moyen le plus élevé (codons préférés)."""
        return self.query(f"""
            SELECT codon, amino_acid,
                   AVG(rscu) AS mean_rscu,
                   SUM(count) AS total_count,
                   COUNT(DISTINCT record_id) AS n_sequences
            FROM codon_rscu
            WHERE amino_acid NOT IN ('Stop', '?')
            GROUP BY codon, amino_acid
            ORDER BY mean_rscu DESC
            LIMIT {n}
        """)

    def close(self) -> None:
        self._con.close()
        logger.debug("DuckDB fermé : %s", self._path)

    def __enter__(self) -> "ARNFinderDB":
        return self

    def __exit__(self, *_: Any) -> None:
        self.close()
