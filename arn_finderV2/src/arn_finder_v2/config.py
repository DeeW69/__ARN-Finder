"""Environment and configuration helpers for the ARN Finder V2 CLI."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, Optional, Union

ENVIRONMENT_VARIABLES = (
    "NCBI_EMAIL",
    "NCBI_TOOL",
    "NCBI_API_KEY",
    "BLAST_PROGRAM",
    "BLAST_DB",
)


@dataclass
class BlastConfig:
    """Structured access to BLAST/NCBI credentials."""

    ncbi_email: Optional[str]
    ncbi_tool: Optional[str]
    ncbi_api_key: Optional[str]
    blast_program: Optional[str]
    blast_db: Optional[str]

    def as_dict(self) -> Dict[str, Optional[str]]:
        return {
            "NCBI_EMAIL": self.ncbi_email,
            "NCBI_TOOL": self.ncbi_tool,
            "NCBI_API_KEY": self.ncbi_api_key,
            "BLAST_PROGRAM": self.blast_program,
            "BLAST_DB": self.blast_db,
        }

    def missing_keys(self) -> Iterator[str]:
        for key, value in self.as_dict().items():
            if not value:
                yield key


PathLike = Union[str, Path]


def find_env_file(start_path: Optional[PathLike] = None) -> Optional[Path]:
    """Search for the closest .env file starting from start_path or CWD."""
    search_root = Path(start_path).resolve() if start_path else Path.cwd().resolve()
    for candidate_dir in _walk_upwards(search_root):
        candidate = candidate_dir / ".env"
        if candidate.exists():
            return candidate

    # Fallback: repository root relative to the installed package (useful in editable installs).
    repo_candidate = Path(__file__).resolve().parents[3] / ".env"
    if repo_candidate.exists():
        return repo_candidate

    return None


def load_dotenv(start_path: Optional[PathLike] = None) -> Optional[Path]:
    """Load .env values into os.environ without overriding pre-existing values."""
    env_path = find_env_file(start_path)
    if not env_path:
        return None

    for key, value in parse_env_file(env_path).items():
        os.environ.setdefault(key, value)

    return env_path


def parse_env_file(env_path: Path) -> Dict[str, str]:
    """Parse a .env file manually to avoid bringing an extra dependency."""
    values: Dict[str, str] = {}
    for raw_line in env_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            continue
        key, raw_value = line.split("=", 1)
        values[key.strip()] = _strip_quotes(raw_value.strip())
    return values


def collect_runtime_config(start_path: Optional[PathLike] = None) -> BlastConfig:
    """Ensure the .env file is sourced and expose the values as a dataclass."""
    load_dotenv(start_path)
    env = os.environ
    return BlastConfig(
        ncbi_email=env.get("NCBI_EMAIL"),
        ncbi_tool=env.get("NCBI_TOOL"),
        ncbi_api_key=env.get("NCBI_API_KEY"),
        blast_program=env.get("BLAST_PROGRAM"),
        blast_db=env.get("BLAST_DB"),
    )


def _walk_upwards(start: Path) -> Iterable[Path]:
    current = start
    last = None
    while last != current:
        yield current
        last = current
        current = current.parent


def _strip_quotes(value: str) -> str:
    if (value.startswith('"') and value.endswith('"')) or (value.startswith("'") and value.endswith("'")):
        return value[1:-1]
    return value
