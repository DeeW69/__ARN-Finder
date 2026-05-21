"""
Embeddings de séquences — V3.4

Trois backends disponibles :
  kmer_tfidf           — TF-IDF sur k-mers (numpy, toujours disponible)
  nucleotide_transformer — InstaDeepAI/nucleotide-transformer-v2-50m-multi-species
  esm2                 — facebook/esm2_t6_8M_UR50D

Usage :
    from arn_finder_v3.embeddings import EmbeddingRequest, compute_embeddings
    req = EmbeddingRequest(input_fasta=Path("seqs.fasta"), out_dir=Path("embeddings"))
    result = compute_embeddings(req)
    # → embeddings/embeddings.parquet  (record_id + embedding list[float32])
    # → embeddings/embedding_manifest.json
"""

from __future__ import annotations

import hashlib
import json
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Optional

from ..io import get_logger, now_iso

if TYPE_CHECKING:
    pass

logger = get_logger("arnfinderv3.embeddings")

# ── Constantes ────────────────────────────────────────────────────────────────

BACKEND_KMER_TFIDF = "kmer_tfidf"
BACKEND_NT = "nucleotide_transformer"
BACKEND_ESM2 = "esm2"

_KNOWN_BACKENDS = {BACKEND_KMER_TFIDF, BACKEND_NT, BACKEND_ESM2}

# Modèles HuggingFace
_NT_MODEL_ID = "InstaDeepAI/nucleotide-transformer-v2-50m-multi-species"
_ESM2_MODEL_ID = "facebook/esm2_t6_8M_UR50D"


# ── Dataclasses ───────────────────────────────────────────────────────────────

@dataclass(slots=True)
class EmbeddingRequest:
    input_fasta: Path
    out_dir: Path

    # Choix du backend : "kmer_tfidf" | "nucleotide_transformer" | "esm2" | "auto"
    backend: str = "auto"

    # Paramètres TF-IDF
    k: int = 6                    # taille des k-mers
    tfidf_vocab_size: int = 4096  # taille max du vocabulaire

    # Paramètres transformers
    batch_size: int = 8           # nombre de séquences par batch
    max_length: int = 512         # max tokens (troncature)
    device: str = "cpu"           # "cpu" | "cuda"

    # Sorties
    save_parquet: bool = True
    save_tsv: bool = False        # TSV optionnel (plus gros)


@dataclass
class EmbeddingResult:
    embeddings_path: Path         # Parquet principal
    manifest_path: Path
    total_records: int
    embedding_dim: int
    backend_used: str
    elapsed_s: float
    tsv_path: Optional[Path] = None


# ── Entrée principale ─────────────────────────────────────────────────────────

def compute_embeddings(request: EmbeddingRequest) -> EmbeddingResult:
    """Calcule les embeddings pour toutes les séquences du FASTA."""
    t0 = time.perf_counter()

    if request.backend not in _KNOWN_BACKENDS and request.backend != "auto":
        raise ValueError(
            f"Backend inconnu : {request.backend!r}. "
            f"Valeurs acceptées : {sorted(_KNOWN_BACKENDS | {'auto'})}"
        )

    request.out_dir.mkdir(parents=True, exist_ok=True)

    # ── Lecture FASTA ────────────────────────────────────────────────────────
    records = _read_fasta(request.input_fasta)
    if not records:
        raise ValueError(f"Aucune séquence trouvée dans {request.input_fasta}")

    ids = [r[0] for r in records]
    seqs = [r[1] for r in records]
    logger.info("Embeddings : %d séquences, backend=%s", len(ids), request.backend)

    # ── Sélection du backend ─────────────────────────────────────────────────
    backend_used, embeddings = _dispatch(request, ids, seqs)

    dim = len(embeddings[0]) if embeddings else 0

    # ── Export Parquet ────────────────────────────────────────────────────────
    parquet_path = request.out_dir / "embeddings.parquet"
    if request.save_parquet:
        _save_parquet(ids, embeddings, parquet_path)

    # ── Export TSV optionnel ──────────────────────────────────────────────────
    tsv_path: Optional[Path] = None
    if request.save_tsv:
        tsv_path = request.out_dir / "embeddings.tsv"
        _save_tsv(ids, embeddings, tsv_path)

    elapsed = time.perf_counter() - t0

    # ── Manifest ──────────────────────────────────────────────────────────────
    manifest_path = request.out_dir / "embedding_manifest.json"
    manifest = {
        "input_fasta": str(request.input_fasta),
        "backend": backend_used,
        "total_records": len(ids),
        "embedding_dim": dim,
        "k": request.k if backend_used == BACKEND_KMER_TFIDF else None,
        "tfidf_vocab_size": request.tfidf_vocab_size if backend_used == BACKEND_KMER_TFIDF else None,
        "max_length": request.max_length if backend_used != BACKEND_KMER_TFIDF else None,
        "elapsed_s": round(elapsed, 3),
        "generated_at": now_iso(),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    logger.info(
        "Embeddings terminés : %d séquences, dim=%d, backend=%s, %.1fs",
        len(ids), dim, backend_used, elapsed,
    )
    return EmbeddingResult(
        embeddings_path=parquet_path,
        manifest_path=manifest_path,
        total_records=len(ids),
        embedding_dim=dim,
        backend_used=backend_used,
        elapsed_s=elapsed,
        tsv_path=tsv_path,
    )


# ── Dispatch backend ──────────────────────────────────────────────────────────

def _dispatch(
    req: EmbeddingRequest,
    ids: list[str],
    seqs: list[str],
) -> tuple[str, list[list[float]]]:
    """Résout le backend effectif et calcule les embeddings."""
    backend = req.backend

    if backend == "auto":
        backend = _detect_best_backend()
        logger.info("Backend auto-détecté : %s", backend)

    if backend == BACKEND_NT:
        try:
            return BACKEND_NT, _embed_nucleotide_transformer(req, ids, seqs)
        except Exception as exc:
            logger.warning(
                "Nucleotide Transformer indisponible (%s) — fallback kmer_tfidf", exc
            )
            backend = BACKEND_KMER_TFIDF

    if backend == BACKEND_ESM2:
        try:
            return BACKEND_ESM2, _embed_esm2(req, ids, seqs)
        except Exception as exc:
            logger.warning("ESM-2 indisponible (%s) — fallback kmer_tfidf", exc)
            backend = BACKEND_KMER_TFIDF

    return BACKEND_KMER_TFIDF, _embed_kmer_tfidf(req, seqs)


def _detect_best_backend() -> str:
    """Sélectionne le meilleur backend disponible."""
    try:
        import transformers  # noqa: F401
        import torch  # noqa: F401
        return BACKEND_NT
    except ImportError:
        pass
    return BACKEND_KMER_TFIDF


# ── Backend k-mer TF-IDF ──────────────────────────────────────────────────────

def _kmers_from_seq(seq: str, k: int) -> list[str]:
    seq = seq.upper().replace("U", "T")
    return [seq[i:i+k] for i in range(len(seq) - k + 1)]


def _embed_kmer_tfidf(req: EmbeddingRequest, seqs: list[str]) -> list[list[float]]:
    """
    Embedding TF-IDF sur k-mers.

    1. Construit le vocabulaire à partir des top-N k-mers les plus fréquents
    2. Calcule la matrice TF (term frequency normalisée par séquence)
    3. Applique IDF (log(N / df + 1))
    4. L2-normalise chaque vecteur
    """
    import math

    k = req.k
    vocab_size = req.tfidf_vocab_size
    n = len(seqs)

    # ── Vocabulaire ───────────────────────────────────────────────────────────
    global_counts: Counter[str] = Counter()
    per_seq_kmers: list[list[str]] = []
    for seq in seqs:
        kmers = _kmers_from_seq(seq, k)
        per_seq_kmers.append(kmers)
        global_counts.update(set(kmers))   # document frequency

    # Top-N k-mers par document frequency
    vocab_tokens = [t for t, _ in global_counts.most_common(vocab_size)]
    vocab = {t: i for i, t in enumerate(vocab_tokens)}
    V = len(vocab)

    if V == 0:
        # Séquences trop courtes → vecteur nul
        return [[0.0] * max(vocab_size, 1) for _ in seqs]

    # IDF
    df = [0.0] * V
    for kmers in per_seq_kmers:
        seen = set(kmers)
        for t in seen:
            if t in vocab:
                df[vocab[t]] += 1.0
    idf = [math.log((n + 1.0) / (d + 1.0)) + 1.0 for d in df]

    # TF-IDF + L2 norm
    result: list[list[float]] = []
    for kmers in per_seq_kmers:
        if not kmers:
            result.append([0.0] * V)
            continue
        # TF (raw count normalisé par longueur)
        counts: Counter[str] = Counter(kmers)
        total = len(kmers)
        vec = [0.0] * V
        for t, cnt in counts.items():
            if t in vocab:
                vec[vocab[t]] = (cnt / total) * idf[vocab[t]]

        # L2 norm
        norm = math.sqrt(sum(v * v for v in vec))
        if norm > 0:
            vec = [v / norm for v in vec]
        result.append(vec)

    logger.debug("kmer_tfidf : k=%d vocab=%d dim=%d", k, vocab_size, V)
    return result


# ── Backend Nucleotide Transformer ────────────────────────────────────────────

def _embed_nucleotide_transformer(
    req: EmbeddingRequest,
    ids: list[str],
    seqs: list[str],
) -> list[list[float]]:
    """
    Embedding via InstaDeepAI/nucleotide-transformer-v2-50m-multi-species.
    Nécessite : torch + transformers
    """
    import torch
    from transformers import AutoTokenizer, AutoModel

    logger.info("Chargement du modèle Nucleotide Transformer : %s", _NT_MODEL_ID)
    tokenizer = AutoTokenizer.from_pretrained(_NT_MODEL_ID, trust_remote_code=True)
    model = AutoModel.from_pretrained(_NT_MODEL_ID, trust_remote_code=True)
    model.eval()
    device = torch.device(req.device)
    model.to(device)

    return _run_transformer_batched(model, tokenizer, seqs, req, device)


# ── Backend ESM-2 ─────────────────────────────────────────────────────────────

def _embed_esm2(
    req: EmbeddingRequest,
    ids: list[str],
    seqs: list[str],
) -> list[list[float]]:
    """
    Embedding via facebook/esm2_t6_8M_UR50D.
    Nécessite : torch + transformers
    ESM-2 est pensé pour les protéines — utile si les séquences sont des CDS traductibles.
    """
    import torch
    from transformers import EsmTokenizer, EsmModel

    logger.info("Chargement du modèle ESM-2 : %s", _ESM2_MODEL_ID)
    tokenizer = EsmTokenizer.from_pretrained(_ESM2_MODEL_ID)
    model = EsmModel.from_pretrained(_ESM2_MODEL_ID)
    model.eval()
    device = torch.device(req.device)
    model.to(device)

    return _run_transformer_batched(model, tokenizer, seqs, req, device)


def _run_transformer_batched(
    model: object,
    tokenizer: object,
    seqs: list[str],
    req: EmbeddingRequest,
    device: object,
) -> list[list[float]]:
    """Mean-pooling sur la dernière couche cachée, en batches."""
    import torch

    all_embeddings: list[list[float]] = []
    batch_size = req.batch_size

    for start in range(0, len(seqs), batch_size):
        batch_seqs = seqs[start:start + batch_size]
        inputs = tokenizer(  # type: ignore[operator]
            batch_seqs,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=req.max_length,
        )
        inputs = {k: v.to(device) for k, v in inputs.items()}  # type: ignore[union-attr]

        with torch.no_grad():
            outputs = model(**inputs)  # type: ignore[operator]

        # Mean pooling (masque d'attention pour ignorer le padding)
        hidden = outputs.last_hidden_state            # (B, L, D)
        mask = inputs["attention_mask"].unsqueeze(-1)  # (B, L, 1)
        summed = (hidden * mask).sum(dim=1)            # (B, D)
        counts = mask.sum(dim=1).clamp(min=1e-9)       # (B, 1)
        mean_pooled = (summed / counts).cpu().float()  # (B, D)

        for vec in mean_pooled:
            all_embeddings.append(vec.tolist())

        logger.debug(
            "Batch %d/%d traité (%d séquences)",
            start // batch_size + 1, (len(seqs) - 1) // batch_size + 1, len(batch_seqs),
        )

    return all_embeddings


# ── I/O ───────────────────────────────────────────────────────────────────────

def _read_fasta(fasta_path: Path) -> list[tuple[str, str]]:
    """Lit un fichier FASTA et retourne [(record_id, sequence), ...]."""
    records: list[tuple[str, str]] = []
    current_id: Optional[str] = None
    current_seq: list[str] = []

    with fasta_path.open(encoding="utf-8") as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(current_seq)))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())

    if current_id is not None:
        records.append((current_id, "".join(current_seq)))

    return records


def _save_parquet(
    ids: list[str],
    embeddings: list[list[float]],
    out_path: Path,
) -> None:
    """Sauvegarde les embeddings en Parquet via pyarrow."""
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError:
        raise ImportError("pyarrow non installé. Lance : pip install pyarrow>=15.0") from None

    table = pa.table({
        "record_id": pa.array(ids, type=pa.string()),
        "embedding": pa.array(
            [emb for emb in embeddings],
            type=pa.list_(pa.float32()),
        ),
    })
    pq.write_table(table, str(out_path), compression="snappy")
    logger.info("Embeddings Parquet sauvegardés → %s (%d lignes)", out_path, len(ids))


def _save_tsv(
    ids: list[str],
    embeddings: list[list[float]],
    out_path: Path,
) -> None:
    """Sauvegarde les embeddings en TSV (record_id + valeurs séparées par tabulation)."""
    with out_path.open("w", encoding="utf-8") as f:
        for rid, emb in zip(ids, embeddings):
            row = rid + "\t" + "\t".join(f"{v:.6f}" for v in emb)
            f.write(row + "\n")
    logger.info("Embeddings TSV sauvegardés → %s", out_path)
