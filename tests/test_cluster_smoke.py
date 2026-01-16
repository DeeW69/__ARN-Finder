"""Smoke test for clustering command."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PKG_ROOT = ROOT


def _run_cli(args: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = f"{PKG_ROOT}{os.pathsep}{env.get('PYTHONPATH', '')}"
    return subprocess.run(
        [sys.executable, "-m", "arn_finder.cli", *args],
        cwd=cwd or ROOT,
        env=env,
        check=False,
        text=True,
        capture_output=True,
    )


def _write_fasta(path: Path) -> None:
    path.write_text(
        """>A
ATGCGATGCA
>B
ATGCGATGCC
>C
CCCCGGGGTT
>D
CCCGGGGTTT
""",
        encoding="utf-8",
    )


def test_cluster_jaccard(tmp_path: Path) -> None:
    fasta = tmp_path / "seqs.fasta"
    _write_fasta(fasta)
    out_dir = tmp_path / "clusters"
    proc = _run_cli(
        [
            "cluster",
            "--in-fasta",
            str(fasta),
            "--out-dir",
            str(out_dir),
            "--k",
            "3",
            "--min-sim",
            "0.3",
            "--max-pairs",
            "100",
        ],
        cwd=tmp_path,
    )
    assert proc.returncode == 0, proc.stderr
    clusters_jsonl = out_dir / "clusters.jsonl"
    edges_csv = out_dir / "cluster_edges.csv"
    assert clusters_jsonl.exists()
    assert edges_csv.exists()

    clusters = []
    with clusters_jsonl.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                clusters.append(json.loads(line))
    assert len(clusters) == 2
    cluster_sizes = sorted(cluster["size"] for cluster in clusters)
    assert cluster_sizes == [2, 2]
