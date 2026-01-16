"""Network-assisted smoke test for the fetch command."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

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


@pytest.mark.skipif(os.environ.get("ARNFINDER_SKIP_NETWORK") == "1", reason="network access disabled")
def test_fetch_small_query(tmp_path: Path) -> None:
    out_dir = tmp_path / "data" / "raw"
    cache_dir = tmp_path / "cache"
    proc = _run_cli(
        [
            "fetch",
            "--query",
            "Mammuthus primigenius[Organism] AND mitochondrion[Filter] AND complete genome",
            "--limit",
            "2",
            "--out-dir",
            str(out_dir),
            "--cache-dir",
            str(cache_dir),
            "--sleep",
            "0.34",
        ],
        cwd=tmp_path,
    )
    assert proc.returncode == 0, proc.stderr
    sequences = out_dir / "sequences.fasta"
    metadata = out_dir / "metadata.jsonl"
    assert sequences.exists() and sequences.stat().st_size > 0
    assert metadata.exists() and metadata.stat().st_size > 0
