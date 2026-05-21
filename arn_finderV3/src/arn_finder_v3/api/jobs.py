"""Gestionnaire de jobs en mémoire — état, progress et annulation."""

from __future__ import annotations

import asyncio
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Optional

from .models import JobStatus, StageProgress


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass
class Job:
    job_id: str
    status: JobStatus = JobStatus.PENDING
    stages: list[StageProgress] = field(default_factory=list)
    stage_results: dict[str, Any] = field(default_factory=dict)
    output_dir: str = "data"
    error: Optional[str] = None
    created_at: str = field(default_factory=_now)
    updated_at: str = field(default_factory=_now)
    elapsed_s: Optional[float] = None
    _task: Optional[asyncio.Task] = field(default=None, repr=False, compare=False)

    def touch(self) -> None:
        self.updated_at = _now()

    def set_stage(self, stage: str, status: JobStatus, elapsed: Optional[float] = None, error: Optional[str] = None) -> None:
        for s in self.stages:
            if s.stage == stage:
                s.status = status
                s.elapsed_s = elapsed
                s.error = error
                self.touch()
                return
        self.stages.append(StageProgress(stage=stage, status=status, elapsed_s=elapsed, error=error))
        self.touch()


class JobRegistry:
    """Registre en mémoire des jobs (remplacé par Redis/DB en production)."""

    def __init__(self) -> None:
        self._jobs: dict[str, Job] = {}

    def create(self, output_dir: str = "data") -> Job:
        job_id = str(uuid.uuid4())[:8]
        job = Job(job_id=job_id, output_dir=output_dir)
        self._jobs[job_id] = job
        return job

    def get(self, job_id: str) -> Optional[Job]:
        return self._jobs.get(job_id)

    def all(self) -> list[Job]:
        return list(self._jobs.values())

    def cancel(self, job_id: str) -> bool:
        job = self._jobs.get(job_id)
        if not job:
            return False
        if job._task and not job._task.done():
            job._task.cancel()
        job.status = JobStatus.CANCELLED
        job.touch()
        return True


# Instance globale partagée par l'app
registry = JobRegistry()
