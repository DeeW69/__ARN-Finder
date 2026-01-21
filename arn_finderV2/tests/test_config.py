from pathlib import Path

from arn_finder_v2 import config


def test_collect_runtime_config_reads_env(tmp_path, monkeypatch):
    env_file = tmp_path / ".env"
    env_file.write_text(
        "\n".join(
            [
                "NCBI_EMAIL=test@example.org",
                "NCBI_TOOL=arnfinderv2",
                "NCBI_API_KEY=abc123",
                "BLAST_PROGRAM=blastn",
                "BLAST_DB=nt",
            ]
        ),
        encoding="utf-8",
    )

    nested_dir = tmp_path / "nested" / "deeper"
    nested_dir.mkdir(parents=True)

    cfg = config.collect_runtime_config(start_path=nested_dir)

    assert cfg.ncbi_email == "test@example.org"
    assert cfg.ncbi_tool == "arnfinderv2"
    assert list(cfg.missing_keys()) == []


def test_missing_keys_reports_unset_variables(monkeypatch):
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    cfg = config.BlastConfig(None, "tool", None, "blastn", "nt")
    assert list(cfg.missing_keys()) == ["NCBI_EMAIL", "NCBI_API_KEY"]
