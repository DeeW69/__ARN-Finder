"""Tests du job store SQLite V3.2."""

from __future__ import annotations

import pytest

from arn_finder_v3.api.jobs import Job, JobRegistry, SQLiteJobStore
from arn_finder_v3.api.models import JobStatus


@pytest.fixture()
def store() -> SQLiteJobStore:
    """Store en mémoire pour les tests."""
    return SQLiteJobStore(":memory:")


@pytest.fixture()
def registry_mem() -> JobRegistry:
    """Registre avec SQLite in-memory pour les tests."""
    return JobRegistry(":memory:")


# ── SQLiteJobStore ─────────────────────────────────────────────────────────────

def test_store_create_and_get(store: SQLiteJobStore) -> None:
    job = store.create("data/test")
    assert len(job.job_id) == 8
    loaded = store.get(job.job_id)
    assert loaded is not None
    assert loaded.job_id == job.job_id
    assert loaded.status == JobStatus.PENDING


def test_store_save_status(store: SQLiteJobStore) -> None:
    job = store.create()
    job.status = JobStatus.RUNNING
    job.updated_at = "2026-01-01T00:00:00"
    store.save(job)

    loaded = store.get(job.job_id)
    assert loaded.status == JobStatus.RUNNING


def test_store_save_stages(store: SQLiteJobStore) -> None:
    from arn_finder_v3.api.models import StageProgress

    job = store.create()
    job.stages.append(StageProgress(stage="fetch", status=JobStatus.COMPLETED, elapsed_s=1.5))
    store.save(job)

    loaded = store.get(job.job_id)
    assert len(loaded.stages) == 1
    assert loaded.stages[0].stage == "fetch"
    assert loaded.stages[0].status == JobStatus.COMPLETED


def test_store_all(store: SQLiteJobStore) -> None:
    store.create()
    store.create()
    store.create()
    jobs = store.all()
    assert len(jobs) == 3


def test_store_clear(store: SQLiteJobStore) -> None:
    store.create()
    store.create()
    store.clear()
    assert store.all() == []


def test_store_get_nonexistent(store: SQLiteJobStore) -> None:
    assert store.get("nonexistent") is None


def test_store_persist_elapsed(store: SQLiteJobStore) -> None:
    job = store.create()
    job.elapsed_s = 42.5
    store.save(job)
    loaded = store.get(job.job_id)
    assert loaded.elapsed_s == pytest.approx(42.5)


def test_store_persist_error(store: SQLiteJobStore) -> None:
    job = store.create()
    job.error = "Test error message"
    store.save(job)
    loaded = store.get(job.job_id)
    assert loaded.error == "Test error message"


# ── JobRegistry ────────────────────────────────────────────────────────────────

def test_registry_create_and_get(registry_mem: JobRegistry) -> None:
    job = registry_mem.create("data/out")
    fetched = registry_mem.get(job.job_id)
    assert fetched is not None
    assert fetched.job_id == job.job_id


def test_registry_auto_persist_on_touch(registry_mem: JobRegistry) -> None:
    """touch() doit déclencher la persistence SQLite automatiquement."""
    job = registry_mem.create()
    job.status = JobStatus.RUNNING
    job.touch()

    # Vide le cache pour forcer la lecture depuis SQLite
    registry_mem._cache.clear()
    loaded = registry_mem.get(job.job_id)
    assert loaded is not None
    assert loaded.status == JobStatus.RUNNING


def test_registry_cancel(registry_mem: JobRegistry) -> None:
    job = registry_mem.create()
    job.status = JobStatus.RUNNING
    job.touch()
    result = registry_mem.cancel(job.job_id)
    assert result is True

    # Vide le cache
    registry_mem._cache.clear()
    loaded = registry_mem.get(job.job_id)
    assert loaded.status == JobStatus.CANCELLED


def test_registry_cancel_nonexistent(registry_mem: JobRegistry) -> None:
    assert registry_mem.cancel("unknown") is False


def test_registry_all(registry_mem: JobRegistry) -> None:
    registry_mem.create()
    registry_mem.create()
    jobs = registry_mem.all()
    assert len(jobs) == 2


def test_registry_clear(registry_mem: JobRegistry) -> None:
    registry_mem.create()
    registry_mem._clear()
    assert registry_mem.all() == []
    assert registry_mem._cache == {}


def test_registry_persist_across_instances() -> None:
    """Simule un redémarrage : le second registre recharge depuis SQLite."""
    import tempfile, os
    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as f:
        db_path = f.name
    try:
        reg1 = JobRegistry(db_path)
        job = reg1.create()
        job.status = JobStatus.COMPLETED
        job.elapsed_s = 10.5
        job.touch()

        # "Redémarrage" : nouvelle instance, même fichier
        reg2 = JobRegistry(db_path)
        loaded = reg2.get(job.job_id)
        assert loaded is not None
        assert loaded.status == JobStatus.COMPLETED
        assert loaded.elapsed_s == pytest.approx(10.5)
    finally:
        os.unlink(db_path)
