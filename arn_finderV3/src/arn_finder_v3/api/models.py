"""Modèles Pydantic pour l'API FastAPI."""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Any, Optional

from pydantic import BaseModel, Field


class JobStatus(str, Enum):
    PENDING   = "pending"
    RUNNING   = "running"
    COMPLETED = "completed"
    FAILED    = "failed"
    CANCELLED = "cancelled"


# ── Requêtes ──────────────────────────────────────────────────────────────────

class PipelineRunRequest(BaseModel):
    query: str = Field(..., description="Requête NCBI Entrez")
    out_dir: str = Field("data", description="Répertoire de sortie")
    limit: int = Field(200, ge=1, le=10_000)
    min_len: int = Field(200, ge=1)
    max_n_frac: float = Field(0.05, ge=0.0, le=1.0)
    min_gc: float = Field(0.0, ge=0.0, le=1.0)
    max_gc: float = Field(1.0, ge=0.0, le=1.0)
    dedupe: bool = False
    k: int = Field(9, ge=3, le=21)
    min_sim: float = Field(0.35, ge=0.0, le=1.0)
    run_blast: bool = False
    run_structure: bool = True
    run_export: bool = True
    stop_on_error: bool = True


class StructureRequest(BaseModel):
    fasta_path: str = Field(..., description="Chemin vers le FASTA")
    out_dir: str = Field("data/structure")
    max_len: int = Field(2000, ge=1)
    temperature: float = Field(37.0, ge=0.0, le=100.0)


# ── Réponses ──────────────────────────────────────────────────────────────────

class JobCreatedResponse(BaseModel):
    job_id: str
    status: JobStatus
    message: str


class StageProgress(BaseModel):
    stage: str
    status: JobStatus
    elapsed_s: Optional[float] = None
    error: Optional[str] = None


class JobStatusResponse(BaseModel):
    job_id: str
    status: JobStatus
    stages: list[StageProgress]
    elapsed_s: Optional[float] = None
    error: Optional[str] = None
    created_at: str
    updated_at: str


class JobResultResponse(BaseModel):
    job_id: str
    status: JobStatus
    stage_summaries: dict[str, Any]
    output_dir: str


class HealthResponse(BaseModel):
    status: str = "ok"
    version: str
    backends: dict[str, bool]
