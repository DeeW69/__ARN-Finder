"""Tests de l'API FastAPI V3.1."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import AsyncMock, MagicMock, patch

import pytest
from fastapi.testclient import TestClient

from arn_finder_v3.api.app import create_app
from arn_finder_v3.api.jobs import registry, Job
from arn_finder_v3.api.models import JobStatus


@pytest.fixture()
def client():
    app = create_app()
    with TestClient(app) as c:
        yield c


@pytest.fixture(autouse=True)
def reset_registry():
    """Remet le registre à zéro entre les tests."""
    registry._clear()
    yield
    registry._clear()


def test_health(client: TestClient) -> None:
    r = client.get("/health")
    assert r.status_code == 200
    data = r.json()
    assert data["status"] == "ok"
    assert "nussinov_approx" in data["backends"]
    assert data["backends"]["nussinov_approx"] is True


def test_pipeline_run_returns_job_id(client: TestClient) -> None:
    with patch("arn_finder_v3.api.app._run_pipeline_task", new_callable=AsyncMock):
        r = client.post("/pipeline/run", json={
            "query": "BRCA1[Gene Name]",
            "limit": 10,
            "run_blast": False,
        })
    assert r.status_code == 202
    data = r.json()
    assert "job_id" in data
    assert len(data["job_id"]) == 8


def test_pipeline_status_not_found(client: TestClient) -> None:
    r = client.get("/pipeline/nonexistent/status")
    assert r.status_code == 404


def test_pipeline_status_found(client: TestClient) -> None:
    job = registry.create()
    job.status = JobStatus.RUNNING
    r = client.get(f"/pipeline/{job.job_id}/status")
    assert r.status_code == 200
    assert r.json()["status"] == "running"


def test_pipeline_results_not_ready(client: TestClient) -> None:
    job = registry.create()
    job.status = JobStatus.RUNNING
    r = client.get(f"/pipeline/{job.job_id}/results")
    assert r.status_code == 409


def test_pipeline_results_ready(client: TestClient) -> None:
    job = registry.create()
    job.status = JobStatus.COMPLETED
    job.stage_results = {"filter": {"kept": 5}}
    r = client.get(f"/pipeline/{job.job_id}/results")
    assert r.status_code == 200
    data = r.json()
    assert data["status"] == "completed"
    assert "filter" in data["stage_summaries"]


def test_pipeline_cancel(client: TestClient) -> None:
    job = registry.create()
    job.status = JobStatus.RUNNING
    r = client.delete(f"/pipeline/{job.job_id}")
    assert r.status_code == 200
    assert registry.get(job.job_id).status == JobStatus.CANCELLED


def test_pipeline_list(client: TestClient) -> None:
    registry.create()
    registry.create()
    r = client.get("/pipeline")
    assert r.status_code == 200
    assert len(r.json()) == 2


def test_structure_endpoint_invalid_fasta(client: TestClient) -> None:
    r = client.post("/structure/run", json={
        "fasta_path": "/nonexistent/path.fasta",
        "out_dir": "/tmp/out",
    })
    assert r.status_code == 422
