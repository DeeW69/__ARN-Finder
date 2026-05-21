"""Configuration V3 — Pydantic Settings (remplace les dataclasses V1/V2)."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from pydantic import Field, field_validator
from pydantic_settings import BaseSettings, SettingsConfigDict


class NCBISettings(BaseSettings):
    """Paramètres NCBI lus depuis .env ou variables d'environnement."""

    model_config = SettingsConfigDict(env_prefix="NCBI_", env_file=".env", extra="ignore")

    email: str = Field(..., description="Email pour NCBI (obligatoire)")
    tool: str = Field("arnfinder_v3", description="Nom de l'outil déclaré à NCBI")
    api_key: str | None = Field(None, description="Clé API NCBI (débloque 10 req/s)")

    @property
    def requests_per_second(self) -> float:
        return 10.0 if self.api_key else 3.0


class BlastSettings(BaseSettings):
    model_config = SettingsConfigDict(env_prefix="BLAST_", env_file=".env", extra="ignore")

    program: Literal["blastn", "blastx", "blastp", "tblastn", "tblastx"] = "blastn"
    db: str = "nt"
    max_hits: int = Field(10, ge=1, le=500)
    poll_interval: int = Field(10, ge=5, description="Secondes entre chaque polling")
    timeout: int = Field(600, ge=60, description="Timeout total BLAST en secondes")


class FetchDefaults(BaseSettings):
    model_config = SettingsConfigDict(env_prefix="FETCH_", env_file=".env", extra="ignore")

    limit: int = Field(200, ge=1, le=10_000)
    batch_size: int = Field(50, ge=10, le=500)
    db: str = "nucleotide"
    rettype: str = "fasta"
    retmode: str = "text"


class FilterDefaults(BaseSettings):
    model_config = SettingsConfigDict(env_prefix="FILTER_", env_file=".env", extra="ignore")

    min_len: int = Field(200, ge=1)
    max_len: int | None = None
    max_n_frac: float = Field(0.05, ge=0.0, le=1.0)
    min_gc: float = Field(0.2, ge=0.0, le=1.0)
    max_gc: float = Field(0.8, ge=0.0, le=1.0)
    max_per_organism: int | None = None


class MotifDefaults(BaseSettings):
    model_config = SettingsConfigDict(env_prefix="MOTIF_", env_file=".env", extra="ignore")

    k: int = Field(9, ge=3, le=21)
    top: int = Field(200, ge=1)
    min_count: int = Field(1, ge=1)
    canonical: bool = True


class ClusterDefaults(BaseSettings):
    model_config = SettingsConfigDict(env_prefix="CLUSTER_", env_file=".env", extra="ignore")

    k: int = Field(9, ge=3, le=21)
    min_similarity: float = Field(0.35, ge=0.0, le=1.0)
    max_pairs: int = Field(200_000, ge=1000)


class DataDirs(BaseSettings):
    """Répertoires de données — configurables via env ou code."""

    model_config = SettingsConfigDict(env_prefix="DATA_", env_file=".env", extra="ignore")

    base: Path = Path("data")

    @property
    def raw(self) -> Path:
        return self.base / "raw"

    @property
    def cache(self) -> Path:
        return self.base / "cache"

    @property
    def filtered(self) -> Path:
        return self.base / "filtered"

    @property
    def motifs(self) -> Path:
        return self.base / "motifs"

    @property
    def clusters(self) -> Path:
        return self.base / "clusters"

    @property
    def consensus(self) -> Path:
        return self.base / "consensus"

    @property
    def blast(self) -> Path:
        return self.base / "blast"

    @property
    def export(self) -> Path:
        return self.base / "export"


# Instances globales (chargées une seule fois)
ncbi = NCBISettings()  # type: ignore[call-arg]  # email requis, doit venir du .env
blast = BlastSettings()
fetch_defaults = FetchDefaults()
filter_defaults = FilterDefaults()
motif_defaults = MotifDefaults()
cluster_defaults = ClusterDefaults()
dirs = DataDirs()
