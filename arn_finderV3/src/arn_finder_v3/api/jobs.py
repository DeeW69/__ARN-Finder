"""
Gestionnaire de jobs FastAPI avec persistence SQLite — V3.2

Les jobs survivent aux redémarrages du serveur grâce à un store SQLite
(~/.arnfinder/jobs.db par défaut). Les asyncio.Task sont conservés en
mémoire uniquement (non sérialisables).

Compatibilité V3.1 : l'interface publique de JobRegistry est inchangée.
"""

from __future__ import annotations

import asyncio
import json
import sqlite3
import threading
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Optional

from .models import JobStatus, StageProgress

_DEFAULT_DB: Path = Path.home() / ".arnfinder" / "jobs.db"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


# ── Dataclass Job ─────────────────────────────────────────────────────────────

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
    # Champs non persistés (en mémoire uniquement)
    _task: Optional[asyncio.Task] = field(default=None, repr=False, compare=False)
    _persist_cb: Optional[Callable[["Job"], None]] = field(
        default=None, repr=False, compare=False
    )

    def touch(self) -> None:
        """Met à jour updated_at et persiste si un callback est défini."""
        self.updated_at = _now()
        if self._persist_cb is not None:
            try:
                self._persist_cb(self)
            except Exception:
                pass  # La persistence ne doit jamais casser le pipeline

    def set_stage(
        self,
        stage: str,
        status: JobStatus,
        elapsed: Optional[float] = None,
        error: Optional[str] = None,
    ) -> None:
        for s in self.stages:
            if s.stage == stage:
                s.status = status
                s.elapsed_s = elapsed
                s.error = error
                self.touch()
                return
        self.stages.append(
            StageProgress(stage=stage, status=status, elapsed_s=elapsed, error=error)
        )
        self.touch()


# ── SQLiteJobStore ─────────────────────────────────────────────────────────────

class SQLiteJobStore:
    """
    Store SQLite pour les jobs API.

    Chaque job est sérialisé avec ses étapes (JSON) et ses résultats (JSON).
    Thread-safe grâce à threading.Lock.
    Supporte ":memory:" pour les tests.
    """

    DDL = """
    CREATE TABLE IF NOT EXISTS jobs (
        job_id        TEXT PRIMARY KEY,
        status        TEXT NOT NULL,
        output_dir    TEXT NOT NULL DEFAULT 'data',
        error         TEXT,
        elapsed_s     REAL,
        stages        TEXT NOT NULL DEFAULT '[]',
        stage_results TEXT NOT NULL DEFAULT '{}',
        created_at    TEXT NOT NULL,
        updated_at    TEXT NOT NULL
    );
    """

    def __init__(self, db_path: Path | str = _DEFAULT_DB) -> None:
        self._path = str(db_path)
        if self._path != ":memory:":
            Path(self._path).parent.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()
        # Pour :memory:, on garde une connexion unique persistante
        self._mem_con: Optional[sqlite3.Connection] = None
        if self._path == ":memory:":
            self._mem_con = sqlite3.connect(":memory:", check_same_thread=False)
            self._mem_con.row_factory = sqlite3.Row
            self._mem_con.execute(self.DDL)
            self._mem_con.commit()
        else:
            con = sqlite3.connect(self._path, check_same_thread=False)
            con.row_factory = sqlite3.Row
            con.execute(self.DDL)
            con.commit()
            con.close()

    def _connect(self) -> sqlite3.Connection:
        if self._mem_con is not None:
            return self._mem_con
        con = sqlite3.connect(self._path, check_same_thread=False)
        con.row_factory = sqlite3.Row
        return con

    def _commit_close(self, con: sqlite3.Connection) -> None:
        if self._mem_con is None:
            con.commit()
            con.close()

    # ── CRUD ──────────────────────────────────────────────────────────────────

    def create(self, output_dir: str = "data") -> Job:
        job = Job(job_id=str(uuid.uuid4())[:8], output_dir=output_dir)
        with self._lock:
            con = self._connect()
            con.execute(
                "INSERT INTO jobs VALUES (?,?,?,?,?,?,?,?,?)",
                [
                    job.job_id, job.status.value, job.output_dir,
                    job.error, job.elapsed_s,
                    "[]", "{}",
                    job.created_at, job.updated_at,
                ],
            )
            self._commit_close(con)
        return job

    def save(self, job: Job) -> None:
        """Persiste l'état courant du job (update)."""
        stages_json = json.dumps([s.model_dump() for s in job.stages])
        try:
            results_json = json.dumps(_jsonify(job.stage_results))
        except Exception:
            results_json = "{}"

        with self._lock:
            con = self._connect()
            con.execute(
                """UPDATE jobs
                   SET status=?, error=?, elapsed_s=?, stages=?,
                       stage_results=?, updated_at=?, output_dir=?
                   WHERE job_id=?""",
                [
                    job.status.value, job.error, job.elapsed_s,
                    stages_json, results_json, job.updated_at,
                    job.output_dir, job.job_id,
                ],
            )
            self._commit_close(con)

    def get(self, job_id: str) -> Optional[Job]:
        with self._lock:
            con = self._connect()
            row = con.execute(
                "SELECT * FROM jobs WHERE job_id=?", [job_id]
            ).fetchone()
            self._commit_close(con)
        return self._row_to_job(dict(row)) if row else None

    def all(self) -> list[Job]:
        with self._lock:
            con = self._connect()
            rows = con.execute(
                "SELECT * FROM jobs ORDER BY created_at DESC"
            ).fetchall()
            self._commit_close(con)
        return [self._row_to_job(dict(r)) for r in rows]

    def clear(self) -> None:
        """Supprime tous les jobs — pour les tests uniquement."""
        with self._lock:
            con = self._connect()
            con.execute("DELETE FROM jobs")
            self._commit_close(con)

    # ── Désérialisation ───────────────────────────────────────────────────────

    def _row_to_job(self, row: dict) -> Job:
        job = Job(
            job_id=row["job_id"],
            status=JobStatus(row["status"]),
            output_dir=row.get("output_dir") or "data",
            error=row.get("error"),
            elapsed_s=row.get("elapsed_s"),
            created_at=row["created_at"],
            updated_at=row["updated_at"],
        )
        stages_raw = json.loads(row.get("stages") or "[]")
        for s in stages_raw:
            job.stages.append(StageProgress(**s))
        try:
            job.stage_results = json.loads(row.get("stage_results") or "{}")
        except Exception:
            job.stage_results = {}
        return job


# ── JobRegistry ───────────────────────────────────────────────────────────────

class JobRegistry:
    """
    Registre de jobs — backed par SQLiteJobStore pour la persistence.

    Interface publique identique à V3.1 pour compatibilité avec app.py.
    Les jobs actifs sont également mis en cache en mémoire pour les accès
    SSE rapides (polling updated_at).
    """

    def __init__(self, db_path: Path | str = _DEFAULT_DB) -> None:
        self._store = SQLiteJobStore(db_path)
        # Cache mémoire pour les jobs actifs (évite des SELECT à chaque poll SSE)
        self._cache: dict[str, Job] = {}

    def create(self, output_dir: str = "data") -> Job:
        job = self._store.create(output_dir)
        job._persist_cb = self._store.save  # auto-persist à chaque touch()
        self._cache[job.job_id] = job
        return job

    def get(self, job_id: str) -> Optional[Job]:
        # Cache d'abord (jobs en cours d'exécution)
        if job_id in self._cache:
            return self._cache[job_id]
        # Sinon depuis SQLite (jobs terminés ou rechargés après restart)
        job = self._store.get(job_id)
        if job:
            self._cache[job_id] = job
        return job

    def all(self) -> list[Job]:
        return self._store.all()

    def cancel(self, job_id: str) -> bool:
        job = self.get(job_id)
        if not job:
            return False
        if job._task and not job._task.done():
            job._task.cancel()
        job.status = JobStatus.CANCELLED
        job.touch()  # déclenche _persist_cb → SQLite
        return True

    def _clear(self) -> None:
        """Réinitialise le registre — à utiliser dans les tests uniquement."""
        self._cache.clear()
        self._store.clear()


# ── Instance globale ──────────────────────────────────────────────────────────
# Utilise ~/.arnfinder/jobs.db — crée le répertoire si besoin.
registry = JobRegistry()


# ── Helpers ───────────────────────────────────────────────────────────────────

def _jsonify(obj: Any) -> Any:
    """Rend récursivement un objet JSON-sérialisable."""
    if isinstance(obj, dict):
        return {k: _jsonify(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonify(v) for v in obj]
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, (str, int, float, bool, type(None))):
        return obj
    return str(obj)
