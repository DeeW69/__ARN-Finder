"""Tests embeddings de séquences — V3.4."""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from arn_finder_v3.embeddings import (
    BACKEND_ESM2,
    BACKEND_KMER_TFIDF,
    BACKEND_NT,
    EmbeddingRequest,
    EmbeddingResult,
    compute_embeddings,
    _embed_kmer_tfidf,
    _kmers_from_seq,
)

# ── Helpers ────────────────────────────────────────────────────────────────────

_SEQS = [
    ("SEQ1", "ATGCATGCATGCATGCATGCATGCATGCATGC"),
    ("SEQ2", "ATGCATGCATGCATGCATGCATGCGGGGGGGG"),
    ("SEQ3", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"),
]


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w") as f:
        for sid, seq in records:
            f.write(f">{sid}\n{seq}\n")


# ── Tests k-mers helpers ──────────────────────────────────────────────────────

def test_kmers_from_seq_length() -> None:
    seq = "ATGCATGC"
    kmers = _kmers_from_seq(seq, k=3)
    assert len(kmers) == len(seq) - 3 + 1


def test_kmers_from_seq_content() -> None:
    kmers = _kmers_from_seq("ATGC", k=2)
    assert kmers == ["AT", "TG", "GC"]


def test_kmers_rna_conversion() -> None:
    """U doit être converti en T avant extraction."""
    kmers = _kmers_from_seq("AUGC", k=2)
    assert "AU" not in kmers
    assert "AT" in kmers


def test_kmers_empty_seq() -> None:
    assert _kmers_from_seq("", k=3) == []


def test_kmers_shorter_than_k() -> None:
    assert _kmers_from_seq("AT", k=5) == []


# ── Tests TF-IDF backend ──────────────────────────────────────────────────────

def _make_tfidf_req(**kwargs):  # type: ignore[no-untyped-def]
    """Crée un EmbeddingRequest minimal pour kmer_tfidf."""
    defaults = {"k": 3, "tfidf_vocab_size": 64}
    defaults.update(kwargs)
    return type("Req", (), defaults)()


def test_tfidf_output_length() -> None:
    req = _make_tfidf_req(k=3, tfidf_vocab_size=64)
    seqs = ["ATGCATGC", "TTTTTTTT"]
    vecs = _embed_kmer_tfidf(req, seqs)
    assert len(vecs) == 2
    # La dimension réelle = min(vocab_size, nb_kmers_distincts_dans_le_corpus)
    assert all(len(v) == len(vecs[0]) for v in vecs)
    assert len(vecs[0]) > 0


def test_tfidf_l2_normalized() -> None:
    req = _make_tfidf_req(k=3, tfidf_vocab_size=64)
    seqs = ["ATGCATGCATGC", "GCGCGCGCGCGC"]
    vecs = _embed_kmer_tfidf(req, seqs)
    for vec in vecs:
        norm = math.sqrt(sum(v * v for v in vec))
        # Norme ~1.0 (ou 0 si vecteur nul)
        assert norm == pytest.approx(1.0, abs=1e-5) or norm == pytest.approx(0.0, abs=1e-5)


def test_tfidf_identical_seqs_same_vector() -> None:
    req = _make_tfidf_req(k=3, tfidf_vocab_size=32)
    seq = "ATGCATGCATGC"
    vecs = _embed_kmer_tfidf(req, [seq, seq])
    assert vecs[0] == vecs[1]


def test_tfidf_different_seqs_different_vectors() -> None:
    req = _make_tfidf_req(k=3, tfidf_vocab_size=64)
    seqs = ["ATGCATGCATGCATGC", "TTTTTTTTTTTTTTTT"]
    vecs = _embed_kmer_tfidf(req, seqs)
    assert vecs[0] != vecs[1]


def test_tfidf_single_seq() -> None:
    req = _make_tfidf_req(k=3, tfidf_vocab_size=32)
    vecs = _embed_kmer_tfidf(req, ["ATGCATGCATGC"])
    assert len(vecs) == 1
    # Dim = min(32, nb_kmers_distincts) — pour une séquence courte, peut être < 32
    assert len(vecs[0]) > 0


def test_tfidf_very_short_seq_no_crash() -> None:
    req = _make_tfidf_req(k=9, tfidf_vocab_size=32)
    # Séquence plus courte que k → aucun k-mer → vecteur nul
    vecs = _embed_kmer_tfidf(req, ["ATG"])
    assert len(vecs) == 1
    assert all(v == 0.0 for v in vecs[0])


# ── Tests compute_embeddings (intégration) ────────────────────────────────────

def test_compute_embeddings_kmer_tfidf(tmp_path: Path) -> None:
    fasta = tmp_path / "seqs.fasta"
    _write_fasta(fasta, _SEQS)

    req = EmbeddingRequest(
        input_fasta=fasta,
        out_dir=tmp_path / "emb",
        backend=BACKEND_KMER_TFIDF,
        k=3,
        tfidf_vocab_size=64,
    )
    result = compute_embeddings(req)

    assert isinstance(result, EmbeddingResult)
    assert result.total_records == 3
    # dim = min(vocab_size, distinct_kmers_in_corpus) — toujours > 0
    assert result.embedding_dim > 0
    assert result.backend_used == BACKEND_KMER_TFIDF
    assert result.embeddings_path.exists()
    assert result.manifest_path.exists()
    assert result.elapsed_s > 0


def test_compute_embeddings_parquet_content(tmp_path: Path) -> None:
    """Le Parquet doit contenir record_id et embedding."""
    import pyarrow.parquet as pq

    fasta = tmp_path / "seqs.fasta"
    _write_fasta(fasta, _SEQS)

    req = EmbeddingRequest(
        input_fasta=fasta,
        out_dir=tmp_path / "emb",
        backend=BACKEND_KMER_TFIDF,
        k=3,
        tfidf_vocab_size=32,
    )
    result = compute_embeddings(req)

    table = pq.read_table(str(result.embeddings_path))
    assert "record_id" in table.column_names
    assert "embedding" in table.column_names
    assert table.num_rows == 3
    ids = table.column("record_id").to_pylist()
    assert "SEQ1" in ids
    assert "SEQ3" in ids


def test_compute_embeddings_manifest_content(tmp_path: Path) -> None:
    import json

    fasta = tmp_path / "seqs.fasta"
    _write_fasta(fasta, _SEQS)

    req = EmbeddingRequest(
        input_fasta=fasta,
        out_dir=tmp_path / "emb",
        backend=BACKEND_KMER_TFIDF,
        k=4,
        tfidf_vocab_size=128,
    )
    result = compute_embeddings(req)

    manifest = json.loads(result.manifest_path.read_text())
    assert manifest["backend"] == BACKEND_KMER_TFIDF
    assert manifest["total_records"] == 3
    assert manifest["k"] == 4
    assert manifest["tfidf_vocab_size"] == 128
    # dim = min(128, distinct_kmers) — cohérent avec embedding_dim du résultat
    assert manifest["embedding_dim"] == result.embedding_dim
    assert manifest["embedding_dim"] > 0


def test_compute_embeddings_tsv_option(tmp_path: Path) -> None:
    fasta = tmp_path / "seqs.fasta"
    _write_fasta(fasta, _SEQS)

    req = EmbeddingRequest(
        input_fasta=fasta,
        out_dir=tmp_path / "emb",
        backend=BACKEND_KMER_TFIDF,
        k=3,
        tfidf_vocab_size=16,
        save_tsv=True,
    )
    result = compute_embeddings(req)
    assert result.tsv_path is not None
    assert result.tsv_path.exists()
    lines = result.tsv_path.read_text().strip().splitlines()
    assert len(lines) == 3
    # Chaque ligne : record_id\tv1\tv2\t...
    first_cols = lines[0].split("\t")
    assert first_cols[0] in {"SEQ1", "SEQ2", "SEQ3"}
    # dim = min(16, distinct_kmers) ; on vérifie cohérence avec embedding_dim
    assert len(first_cols) == 1 + result.embedding_dim


def test_compute_embeddings_auto_fallback(tmp_path: Path) -> None:
    """Sans torch/transformers, 'auto' doit retomber sur kmer_tfidf."""
    fasta = tmp_path / "seqs.fasta"
    _write_fasta(fasta, _SEQS)

    req = EmbeddingRequest(
        input_fasta=fasta,
        out_dir=tmp_path / "emb",
        backend="auto",
        k=3,
        tfidf_vocab_size=32,
    )
    result = compute_embeddings(req)
    # kmer_tfidf toujours disponible → doit réussir
    assert result.total_records == 3
    assert result.embeddings_path.exists()


def test_compute_embeddings_invalid_backend(tmp_path: Path) -> None:
    fasta = tmp_path / "seqs.fasta"
    _write_fasta(fasta, _SEQS)

    req = EmbeddingRequest(
        input_fasta=fasta,
        out_dir=tmp_path / "emb",
        backend="inexistant",
    )
    with pytest.raises(ValueError, match="Backend inconnu"):
        compute_embeddings(req)


def test_compute_embeddings_empty_fasta(tmp_path: Path) -> None:
    fasta = tmp_path / "empty.fasta"
    fasta.write_text("")

    req = EmbeddingRequest(
        input_fasta=fasta,
        out_dir=tmp_path / "emb",
        backend=BACKEND_KMER_TFIDF,
    )
    with pytest.raises(ValueError, match="Aucune séquence"):
        compute_embeddings(req)


# ── Tests avec torch (skip si non disponible) ─────────────────────────────────

torch_available = pytest.mark.skipif(
    condition=True,  # remplacer par : not importlib.util.find_spec("torch")
    reason="torch non installé",
)

# Ces tests sont gardés en placeholder — ils seraient activés avec une installation
# complète de torch + transformers.

@pytest.mark.skip(reason="Nécessite torch + transformers installés")
def test_nucleotide_transformer_backend(tmp_path: Path) -> None:
    fasta = tmp_path / "seqs.fasta"
    _write_fasta(fasta, _SEQS)

    req = EmbeddingRequest(
        input_fasta=fasta,
        out_dir=tmp_path / "emb",
        backend=BACKEND_NT,
        batch_size=2,
        max_length=64,
    )
    result = compute_embeddings(req)
    assert result.backend_used == BACKEND_NT
    assert result.embedding_dim > 0


@pytest.mark.skip(reason="Nécessite torch + transformers installés")
def test_esm2_backend(tmp_path: Path) -> None:
    fasta = tmp_path / "seqs.fasta"
    _write_fasta(fasta, _SEQS)

    req = EmbeddingRequest(
        input_fasta=fasta,
        out_dir=tmp_path / "emb",
        backend=BACKEND_ESM2,
        batch_size=2,
        max_length=64,
    )
    result = compute_embeddings(req)
    assert result.backend_used == BACKEND_ESM2
    assert result.embedding_dim > 0
