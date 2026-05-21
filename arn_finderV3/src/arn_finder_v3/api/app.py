"""Application FastAPI — ARN Finder V3 API."""

from __future__ import annotations

import asyncio
import json
import time
from pathlib import Path
from typing import AsyncGenerator

from fastapi import FastAPI, HTTPException, BackgroundTasks
from fastapi.responses import StreamingResponse

from .. import __version__
from ..io import get_logger
from .jobs import Job, JobRegistry, registry
from .models import (
    HealthResponse, JobCreatedResponse, JobResultResponse,
    JobStatus, JobStatusResponse, PipelineRunRequest, StructureRequest,
)

logger = get_logger("arnfinderv3.api")


def create_app() -> FastAPI:
    app = FastAPI(
        title="ARN Finder V3 API",
        description="Pipeline bioinformatique ML-ready — interface HTTP",
        version=__version__,
        docs_url="/docs",
        redoc_url="/redoc",
    )

    # ── Health ────────────────────────────────────────────────────────────────

    @app.get("/health", response_model=HealthResponse, tags=["system"])
    async def health() -> HealthResponse:
        """Vérifie l'état de l'API et des backends disponibles."""
        backends: dict[str, bool] = {}
        try:
            import RNA  # noqa: F401
            backends["viennarna"] = True
        except ImportError:
            backends["viennarna"] = False

        import subprocess
        try:
            r = subprocess.run(["RNAfold", "--version"], capture_output=True, timeout=3)
            backends["rnafold_subprocess"] = r.returncode == 0
        except Exception:
            backends["rnafold_subprocess"] = False

        backends["nussinov_approx"] = True  # toujours disponible
        return HealthResponse(version=__version__, backends=backends)

    # ── Pipeline run ──────────────────────────────────────────────────────────

    @app.post("/pipeline/run", response_model=JobCreatedResponse, status_code=202, tags=["pipeline"])
    async def pipeline_run(req: PipelineRunRequest) -> JobCreatedResponse:
        """
        Lance le pipeline complet en tâche de fond.
        Retourne immédiatement un `job_id` pour suivre la progression.
        """
        job = registry.create(output_dir=req.out_dir)
        job.status = JobStatus.RUNNING
        job.touch()

        task = asyncio.create_task(_run_pipeline_task(job, req))
        job._task = task

        logger.info("Job %s créé pour query=%r", job.job_id, req.query)
        return JobCreatedResponse(
            job_id=job.job_id,
            status=JobStatus.PENDING,
            message=f"Pipeline démarré. Suivez la progression sur /pipeline/{job.job_id}/status",
        )

    # ── Status ────────────────────────────────────────────────────────────────

    @app.get("/pipeline/{job_id}/status", response_model=JobStatusResponse, tags=["pipeline"])
    async def pipeline_status(job_id: str) -> JobStatusResponse:
        """Retourne l'état courant d'un job."""
        job = _get_job_or_404(job_id)
        return JobStatusResponse(
            job_id=job.job_id,
            status=job.status,
            stages=job.stages,
            elapsed_s=job.elapsed_s,
            error=job.error,
            created_at=job.created_at,
            updated_at=job.updated_at,
        )

    # ── Résultats ─────────────────────────────────────────────────────────────

    @app.get("/pipeline/{job_id}/results", response_model=JobResultResponse, tags=["pipeline"])
    async def pipeline_results(job_id: str) -> JobResultResponse:
        """Retourne les résultats complets d'un job terminé."""
        job = _get_job_or_404(job_id)
        if job.status not in (JobStatus.COMPLETED, JobStatus.FAILED):
            raise HTTPException(status_code=409, detail=f"Job {job_id} non terminé (statut : {job.status})")
        return JobResultResponse(
            job_id=job.job_id,
            status=job.status,
            stage_summaries=_serialize_results(job),
            output_dir=job.output_dir,
        )

    # ── SSE stream ────────────────────────────────────────────────────────────

    @app.get("/pipeline/{job_id}/stream", tags=["pipeline"])
    async def pipeline_stream(job_id: str) -> StreamingResponse:
        """
        Server-Sent Events — envoie un événement JSON à chaque mise à jour du job.
        Compatible avec EventSource côté client / notebook.
        """
        job = _get_job_or_404(job_id)

        async def event_generator() -> AsyncGenerator[str, None]:
            last_update = ""
            while True:
                if job.updated_at != last_update:
                    last_update = job.updated_at
                    data = json.dumps({
                        "job_id": job.job_id,
                        "status": job.status,
                        "stages": [s.model_dump() for s in job.stages],
                        "elapsed_s": job.elapsed_s,
                    })
                    yield f"data: {data}\n\n"

                if job.status in (JobStatus.COMPLETED, JobStatus.FAILED, JobStatus.CANCELLED):
                    yield "event: done\ndata: {}\n\n"
                    break

                await asyncio.sleep(1)

        return StreamingResponse(event_generator(), media_type="text/event-stream")

    # ── Annulation ────────────────────────────────────────────────────────────

    @app.delete("/pipeline/{job_id}", tags=["pipeline"])
    async def pipeline_cancel(job_id: str) -> dict:
        """Annule un job en cours."""
        if not registry.cancel(job_id):
            raise HTTPException(status_code=404, detail=f"Job {job_id} introuvable.")
        return {"job_id": job_id, "status": "cancelled"}

    # ── Liste des jobs ────────────────────────────────────────────────────────

    @app.get("/pipeline", tags=["pipeline"])
    async def pipeline_list() -> list[dict]:
        """Liste tous les jobs (limité à la session courante)."""
        return [
            {"job_id": j.job_id, "status": j.status, "created_at": j.created_at}
            for j in registry.all()
        ]

    # ── Structure secondaire (endpoint dédié) ─────────────────────────────────

    @app.post("/structure/run", response_model=JobCreatedResponse, status_code=202, tags=["structure"])
    async def structure_run(req: StructureRequest) -> JobCreatedResponse:
        """Lance la prédiction de structure secondaire sur un FASTA existant."""
        fasta = Path(req.fasta_path)
        if not fasta.exists():
            raise HTTPException(status_code=422, detail=f"FASTA introuvable : {req.fasta_path}")

        job = registry.create(output_dir=req.out_dir)
        job.status = JobStatus.RUNNING
        task = asyncio.create_task(_run_structure_task(job, req))
        job._task = task

        return JobCreatedResponse(
            job_id=job.job_id,
            status=JobStatus.PENDING,
            message=f"Structure job démarré. Suivez sur /pipeline/{job.job_id}/status",
        )

    return app


# ── Tâches de fond ────────────────────────────────────────────────────────────

async def _run_pipeline_task(job: Job, req: PipelineRunRequest) -> None:
    """Exécute le pipeline complet dans un thread séparé (non-bloquant)."""
    from ..pipeline import PipelineConfig, run_pipeline

    t0 = time.perf_counter()
    try:
        cfg = PipelineConfig(
            query=req.query,
            out_base=Path(req.out_dir),
            fetch_limit=req.limit,
            filter_min_len=req.min_len,
            filter_max_n_frac=req.max_n_frac,
            filter_min_gc=req.min_gc,
            filter_max_gc=req.max_gc,
            filter_dedupe=req.dedupe,
            motif_k=req.k,
            cluster_k=req.k,
            cluster_min_sim=req.min_sim,
            run_blast=req.run_blast,
            run_structure=req.run_structure,
            run_export=req.run_export,
            stop_on_error=req.stop_on_error,
        )

        # Lance dans un executor pour ne pas bloquer la boucle asyncio
        loop = asyncio.get_running_loop()
        result = await loop.run_in_executor(None, run_pipeline, cfg)

        for stage, elapsed in result.elapsed.items():
            if stage in result.errors:
                job.set_stage(stage, JobStatus.FAILED, elapsed, str(result.errors[stage]))
            else:
                job.set_stage(stage, JobStatus.COMPLETED, elapsed)

        job.stage_results = _serialize_results_from_pipeline(result)
        job.status = JobStatus.FAILED if result.errors else JobStatus.COMPLETED

    except asyncio.CancelledError:
        job.status = JobStatus.CANCELLED
    except Exception as exc:
        job.status = JobStatus.FAILED
        job.error = str(exc)
        logger.exception("Erreur pipeline job %s", job.job_id)
    finally:
        job.elapsed_s = round(time.perf_counter() - t0, 2)
        job.touch()


async def _run_structure_task(job: Job, req: StructureRequest) -> None:
    from ..secondary_structure import StructureRequest as SReq, compute_secondary_structure

    t0 = time.perf_counter()
    try:
        s_req = SReq(
            input_fasta=Path(req.fasta_path),
            out_dir=Path(req.out_dir),
            max_len=req.max_len,
            temperature=req.temperature,
        )
        loop = asyncio.get_running_loop()
        result = await loop.run_in_executor(None, compute_secondary_structure, s_req)
        job.set_stage("structure", JobStatus.COMPLETED, result.elapsed)
        job.stage_results["structure"] = {
            "total": result.total, "computed": result.computed,
            "skipped": result.skipped, "backend": result.backend_used,
        }
        job.status = JobStatus.COMPLETED
    except asyncio.CancelledError:
        job.status = JobStatus.CANCELLED
    except Exception as exc:
        job.set_stage("structure", JobStatus.FAILED, error=str(exc))
        job.status = JobStatus.FAILED
        job.error = str(exc)
    finally:
        job.elapsed_s = round(time.perf_counter() - t0, 2)
        job.touch()


# ── Helpers ───────────────────────────────────────────────────────────────────

def _get_job_or_404(job_id: str) -> Job:
    job = registry.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job {job_id} introuvable.")
    return job


def _serialize_results(job: Job) -> dict:
    return job.stage_results


def _serialize_results_from_pipeline(result) -> dict:
    summaries = {}
    for stage, res in result.stage_results.items():
        try:
            if hasattr(res, "__dataclass_fields__"):
                summaries[stage] = {
                    k: str(v) if isinstance(v, Path) else v
                    for k, v in res.__dict__.items()
                    if not k.startswith("_")
                }
            elif isinstance(res, dict):
                summaries[stage] = res
            else:
                summaries[stage] = str(res)
        except Exception:
            summaries[stage] = "non sérialisable"
    return summaries
