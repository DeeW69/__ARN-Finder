"""Offline smoke test for the filter command."""

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


def _write_test_fasta(path: Path) -> None:
    content = """>PX802414.1 valid record
ACGTACGTACGTACGTACGT
>SHORT_1 short seq
ACG
>MANYN awful
NNNNNNNNNA
>AMBIG1 sample
ACGTUACGTRYSW
>DUPLICATE same as first
acgtacgtacgtacgtacgt
"""
    path.write_text(content, encoding="utf-8")


def _write_test_metadata(path: Path) -> None:
    entries = [
        {"uid": "PX802414.1", "accession": "PX802414.1", "organism": "Troglodytes musculus"},
        {"uid": "AMBIG1", "accession": "AMBIG1", "organism": "Troglodytes musculus"},
        {"uid": "SHORT_1", "accession": "SHORT_1", "organism": "Troglodytes musculus"},
    ]
    with path.open("w", encoding="utf-8") as handle:
        for entry in entries:
            handle.write(json.dumps(entry) + "\n")


def test_filter_command(tmp_path: Path) -> None:
    fasta = tmp_path / "seqs.fasta"
    metadata = tmp_path / "metadata.jsonl"
    _write_test_fasta(fasta)
    _write_test_metadata(metadata)
    out_dir = tmp_path / "filtered"

    proc = _run_cli(
        [
            "filter",
            "--in-fasta",
            str(fasta),
            "--in-metadata",
            str(metadata),
            "--out-dir",
            str(out_dir),
            "--min-len",
            "5",
            "--max-n-frac",
            "0.4",
            "--alphabet",
            "AUTO",
            "--dedupe",
        ],
        cwd=tmp_path,
    )
    assert proc.returncode == 0, proc.stderr
    fasta_out = out_dir / "filtered_sequences.fasta"
    report_out = out_dir / "filter_report.csv"
    metadata_out = out_dir / "filtered_metadata.jsonl"
    manifest = out_dir / "filter_manifest.json"
    assert fasta_out.exists() and fasta_out.read_text().count(">") == 2
    assert metadata_out.exists()

    with report_out.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    rows_by_id = {row["record_id"]: row for row in rows}
    reasons = {rid: row["reason"] for rid, row in rows_by_id.items()}
    assert reasons["PX802414.1"] == "ok"
    assert reasons["SHORT_1"] == "too_short"
    assert "duplicate" in reasons.values()
    assert rows_by_id["DUPLICATE"]["duplicate_of"] == "PX802414.1"

    manifest_data = json.loads(manifest.read_text(encoding="utf-8"))
    assert manifest_data["kept_records"] == 2


def test_filter_max_per_organism_limit(tmp_path: Path) -> None:
    fasta = tmp_path / "orgs.fasta"
    metadata = tmp_path / "orgs.jsonl"
    fasta.write_text(
        """>SEQ1 sample Homo sapiens
ACGTACGTAC
>SEQ2 sample Homo sapiens
ACGTACGTAA
>SEQ3 sample Homo sapiens
ACGTACGTCC
""",
        encoding="utf-8",
    )
    with metadata.open("w", encoding="utf-8") as handle:
        for seq_id in ["SEQ1", "SEQ2", "SEQ3"]:
            handle.write(json.dumps({"uid": seq_id, "accession": seq_id, "organism": "Homo sapiens"}) + "\n")
    out_dir = tmp_path / "filtered_orgs"

    proc = _run_cli(
        [
            "filter",
            "--in-fasta",
            str(fasta),
            "--in-metadata",
            str(metadata),
            "--out-dir",
            str(out_dir),
            "--min-len",
            "5",
            "--max-n-frac",
            "0.5",
            "--max-per-organism",
            "2",
            "--dedupe",
        ],
        cwd=tmp_path,
    )
    assert proc.returncode == 0, proc.stderr
    report_out = out_dir / "filter_report.csv"
    manifest = out_dir / "filter_manifest.json"
    rows = list(csv.DictReader(report_out.open("r", encoding="utf-8")))
    reasons = {row["record_id"]: row["reason"] for row in rows}
    assert list(reasons.values()).count("ok") == 2
    assert "max_per_organism" in reasons.values()
    manifest_data = json.loads(manifest.read_text(encoding="utf-8"))
    assert manifest_data["kept_records"] == 2
    assert manifest_data["dropped_max_per_organism"] == 1
