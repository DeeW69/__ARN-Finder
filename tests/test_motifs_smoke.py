"""Smoke tests for motifs extraction."""

from __future__ import annotations

import csv
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
        """>SEQ1 sample
ATGCATGC
>SEQ2 sample
GCATGCAT
""",
        encoding="utf-8",
    )


def test_motifs_canonical(tmp_path: Path) -> None:
    fasta = tmp_path / "seqs.fasta"
    _write_fasta(fasta)
    out_dir = tmp_path / "motifs"
    proc = _run_cli(
        [
            "motifs",
            "--in-fasta",
            str(fasta),
            "--out-dir",
            str(out_dir),
            "--k",
            "3",
            "--top",
            "10",
            "--canonical",
        ],
        cwd=tmp_path,
    )
    assert proc.returncode == 0, proc.stderr
    motifs_csv = out_dir / "motifs.csv"
    assert motifs_csv.exists()
    rows = list(csv.DictReader(motifs_csv.open("r", encoding="utf-8", newline="")))
    assert rows, "motifs.csv should not be empty"
    kmers = {row["kmer"] for row in rows}
    assert "ATG" in kmers
    assert "CAT" not in kmers  # canonical merges reverse complement
