"""Smoke tests for consensus command."""

from __future__ import annotations

import csv
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
        """>seq1
AACCGGTTAA
>seq2
GGTTAAAACC
>seq3
AACCGGTAAA
""",
        encoding="utf-8",
    )


def _write_clusters(path: Path) -> None:
    clusters = [
        {
            "cluster_id": "C0001",
            "members": [
                {"record_id": "seq1"},
                {"record_id": "seq2"},
                {"record_id": "seq3"},
            ],
            "dominant_organism": "Testus example",
        },
        {
            "cluster_id": "C0002",
            "members": [
                {"record_id": "seq4"},
                {"record_id": "seq5"},
            ],
        },
    ]
    with path.open("w", encoding="utf-8") as handle:
        for cluster in clusters:
            handle.write(json.dumps(cluster) + "\n")


def test_consensus_kmer(tmp_path: Path) -> None:
    fasta = tmp_path / "seqs.fasta"
    clusters = tmp_path / "clusters.jsonl"
    _write_fasta(fasta)
    _write_clusters(clusters)
    out_dir = tmp_path / "consensus"

    proc = _run_cli(
        [
            "consensus",
            "--clusters",
            str(clusters),
            "--in-fasta",
            str(fasta),
            "--out-dir",
            str(out_dir),
            "--min-cluster-size",
            "2",
            "--method",
            "kmer",
            "--min-overlap",
            "4",
            "--min-identity",
            "0.5",
        ],
        cwd=tmp_path,
    )
    assert proc.returncode == 0, proc.stderr
    fasta_out = out_dir / "consensus.fasta"
    stats_csv = out_dir / "consensus_stats.csv"
    assert fasta_out.exists()
    assert stats_csv.exists()
    contents = fasta_out.read_text(encoding="utf-8")
    assert "CONS_C0001" in contents
    rows = list(csv.DictReader(stats_csv.open("r", encoding="utf-8", newline="")))
    assert rows[0]["status"] in {"ok", "low_quality"}
