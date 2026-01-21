"""Central location for default paths and settings."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

# Default directories for each stage of the pipeline
DEFAULT_RAW_DIR = Path("data") / "raw"
DEFAULT_CACHE_DIR = Path("data") / "cache"
DEFAULT_FILTERED_DIR = Path("data") / "filtered"
DEFAULT_MOTIF_DIR = Path("data") / "motifs"
DEFAULT_CLUSTER_DIR = Path("data") / "clusters"
DEFAULT_CONSENSUS_DIR = Path("data") / "consensus"
DEFAULT_BLAST_DIR = Path("data") / "blast"

# Placeholder filenames used by the remaining stubs
FILTER_PLACEHOLDER = "filter_placeholder.json"
MOTIF_PLACEHOLDER = "motifs_placeholder.json"
CLUSTER_PLACEHOLDER = "clusters_placeholder.json"
CONSENSUS_PLACEHOLDER = "consensus_placeholder.json"


@dataclass(slots=True)
class FetchDefaults:
    """Options surfaced on the CLI fetch command."""

    limit: int = 200
    out_dir: Path = DEFAULT_RAW_DIR
    cache_dir: Path = DEFAULT_CACHE_DIR
    db: str = "nuccore"
    rettype: str = "fasta"
    retmode: str = "text"
    tool: str = "arnfinder"
    sleep: float = 0.34
    timeout: int = 30
    retries: int = 3
    batch_size: int = 100


@dataclass(slots=True)
class FilterDefaults:
    """Default values for filter command arguments."""

    min_len: int = 200
    max_n_fraction: float = 0.05
    out_dir: Path = DEFAULT_FILTERED_DIR


@dataclass(slots=True)
class MotifDefaults:
    """Default values for motif scanning stubs."""

    k: int = 9
    out_dir: Path = DEFAULT_MOTIF_DIR


@dataclass(slots=True)
class ClusterDefaults:
    """Defaults for clustering command."""

    k: int = 9
    min_similarity: float = 0.35
    out_dir: Path = DEFAULT_CLUSTER_DIR


@dataclass(slots=True)
class ConsensusDefaults:
    """Defaults for consensus command."""

    out_dir: Path = DEFAULT_CONSENSUS_DIR
    min_cluster_size: int = 2
    method: str = "kmer"
    k: int = 9
    min_overlap: int = 80
    min_identity: float = 0.97
    max_n_frac: float = 0.05


@dataclass(slots=True)
class BlastDefaults:
    """Defaults for BLAST command arguments."""

    out_dir: Path = DEFAULT_BLAST_DIR
    program: str = "blastn"
    db: str = "nt"
    format_type: str = "JSON2"
    max_hits: int = 10
    poll_seconds: int = 10
    timeout_seconds: int = 600
    rate_limit_seconds: float = 1.0


FETCH_DEFAULTS = FetchDefaults()
FILTER_DEFAULTS = FilterDefaults()
MOTIF_DEFAULTS = MotifDefaults()
CLUSTER_DEFAULTS = ClusterDefaults()
CONSENSUS_DEFAULTS = ConsensusDefaults()
BLAST_DEFAULTS = BlastDefaults()
