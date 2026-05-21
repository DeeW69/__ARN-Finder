"""
Structure secondaire ARN — V3.1

Stratégie de backend (ordre de priorité) :
  1. ViennaRNA Python bindings  (`import RNA`)      → conda install -c bioconda viennarna
  2. RNAfold subprocess          (`RNAfold --version`) → inclus dans ViennaRNA system package
  3. Approximation Python pure   (nussinov simplifié)  → toujours disponible, moins précis

Le backend utilisé est tracé dans le manifest de sortie.
Les séquences DNA sont converties en RNA (T→U) avant prédiction.
"""

from __future__ import annotations

import json
import re
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Optional

from rich.progress import Progress, BarColumn, TextColumn, TaskProgressColumn, TimeElapsedColumn

from ..io import ensure_dir, now_iso, read_fasta_records, write_csv, write_json, get_logger

logger = get_logger("arnfinderv3.secondary_structure")

Backend = Literal["viennarna", "rnafold_subprocess", "nussinov_approx"]

# ── Détection du backend disponible ──────────────────────────────────────────

def _detect_backend() -> Backend:
    """Choisit le meilleur backend disponible au runtime."""
    try:
        import RNA  # noqa: F401
        logger.debug("Backend : ViennaRNA Python bindings")
        return "viennarna"
    except ImportError:
        pass

    try:
        result = subprocess.run(
            ["RNAfold", "--version"],
            capture_output=True, text=True, timeout=5,
        )
        if result.returncode == 0:
            logger.debug("Backend : RNAfold subprocess")
            return "rnafold_subprocess"
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass

    logger.warning(
        "ViennaRNA non disponible — utilisation de l'approximation Nussinov (moins précise). "
        "Pour installer ViennaRNA : conda install -c bioconda viennarna"
    )
    return "nussinov_approx"


# ── Dataclasses ───────────────────────────────────────────────────────────────

@dataclass(slots=True)
class StructureRequest:
    input_fasta: Path
    out_dir: Path
    max_len: int = 2000          # Séquences > max_len ignorées (trop lentes)
    temperature: float = 37.0    # °C pour RNAfold
    backend: Optional[Backend] = None   # None = auto-détection


@dataclass(slots=True)
class StructureRecord:
    record_id: str
    sequence_rna: str            # Séquence après T→U
    dot_bracket: str             # Structure en notation point-parenthèse
    mfe: float                   # Énergie libre minimale (kcal/mol)
    length: int
    gc_frac: float
    backend: Backend
    skipped: bool = False
    skip_reason: str = ""


@dataclass
class StructureResult:
    records: list[StructureRecord] = field(default_factory=list)
    structures_csv: Optional[Path] = None
    structures_jsonl: Optional[Path] = None
    manifest_path: Optional[Path] = None
    total: int = 0
    computed: int = 0
    skipped: int = 0
    elapsed: float = 0.0
    backend_used: str = ""


# ── Prédiction par backend ────────────────────────────────────────────────────

def _predict_viennarna(seq_rna: str, temperature: float) -> tuple[str, float]:
    """ViennaRNA Python bindings — le plus précis."""
    import RNA
    md = RNA.md()
    md.temperature = temperature
    fc = RNA.fold_compound(seq_rna, md)
    structure, mfe = fc.mfe()
    return structure, mfe


def _predict_rnafold_subprocess(seq_rna: str, temperature: float) -> tuple[str, float]:
    """Appelle RNAfold en ligne de commande."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False, encoding="utf-8") as f:
        f.write(f">seq\n{seq_rna}\n")
        tmp_path = f.name

    try:
        result = subprocess.run(
            ["RNAfold", "--noPS", f"-T{temperature}", tmp_path],
            capture_output=True, text=True, timeout=120,
        )
        if result.returncode != 0:
            raise RuntimeError(f"RNAfold error: {result.stderr.strip()}")
        lines = [l for l in result.stdout.strip().splitlines() if l]
        # Format : ">seq\n<sequence>\n<structure> (<mfe>)"
        for line in reversed(lines):
            m = re.match(r"^([.()]+)\s+\(\s*(-?\d+\.\d+)\s*\)", line)
            if m:
                return m.group(1), float(m.group(2))
        raise ValueError(f"Format de sortie RNAfold inattendu : {result.stdout[:200]}")
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def _predict_nussinov(seq_rna: str) -> tuple[str, float]:
    """
    Algorithme de Nussinov simplifié — Python pur, sans dépendance.
    Maximise les paires de bases (A-U, G-C, G-U) par programmation dynamique.
    MFE approximée : -1 kcal/mol par paire (valeur symbolique, non physique).
    """
    n = len(seq_rna)
    PAIRS = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}
    MIN_LOOP = 3  # distance minimale entre deux bases appariées

    # Matrice DP
    dp = [[0] * n for _ in range(n)]
    for length in range(MIN_LOOP + 1, n):
        for i in range(n - length):
            j = i + length
            # Option 1 : j non apparié
            dp[i][j] = dp[i][j - 1]
            # Option 2 : j apparié avec k
            for k in range(i, j - MIN_LOOP):
                if (seq_rna[k], seq_rna[j]) in PAIRS:
                    val = (dp[i][k - 1] if k > i else 0) + 1 + dp[k + 1][j - 1]
                    dp[i][j] = max(dp[i][j], val)

    # Traceback
    structure = ["."] * n

    def traceback(i: int, j: int) -> None:
        if i >= j:
            return
        if dp[i][j] == dp[i][j - 1]:
            traceback(i, j - 1)
        else:
            for k in range(i, j - MIN_LOOP):
                if (seq_rna[k], seq_rna[j]) in PAIRS:
                    val = (dp[i][k - 1] if k > i else 0) + 1 + dp[k + 1][j - 1]
                    if dp[i][j] == val:
                        structure[k] = "("
                        structure[j] = ")"
                        if k > i:
                            traceback(i, k - 1)
                        traceback(k + 1, j - 1)
                        return

    if n > 0:
        traceback(0, n - 1)

    dot_bracket = "".join(structure)
    pairs_count = dot_bracket.count("(")
    mfe_approx = -1.0 * pairs_count   # approximation symbolique
    return dot_bracket, mfe_approx


def _predict(seq_rna: str, backend: Backend, temperature: float) -> tuple[str, float]:
    if backend == "viennarna":
        return _predict_viennarna(seq_rna, temperature)
    if backend == "rnafold_subprocess":
        return _predict_rnafold_subprocess(seq_rna, temperature)
    return _predict_nussinov(seq_rna)


# ── Calcul GC ────────────────────────────────────────────────────────────────

def _gc_frac(seq: str) -> float:
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq)


# ── Point d'entrée ────────────────────────────────────────────────────────────

def compute_secondary_structure(request: StructureRequest) -> StructureResult:
    """
    Prédit la structure secondaire de chaque séquence ARN du FASTA.

    Pipeline : DNA/RNA → normalise en RNA (T→U) → ViennaRNA / RNAfold / Nussinov
    Sorties : structures.csv, structures.jsonl, manifest.json
    """
    out_dir = ensure_dir(request.out_dir)
    csv_path = out_dir / "structures.csv"
    jsonl_path = out_dir / "structures.jsonl"
    manifest_path = out_dir / "structure_manifest.json"

    backend = request.backend or _detect_backend()
    t0 = time.perf_counter()
    result = StructureResult(backend_used=backend)

    rows: list[list[str]] = []

    with jsonl_path.open("w", encoding="utf-8") as jf, \
         Progress(
             TextColumn("[cyan]structure[/cyan]"),
             BarColumn(), TaskProgressColumn(), TimeElapsedColumn(), transient=True,
         ) as progress:

        task = progress.add_task("prédiction structure...", total=None)

        for header, sequence in read_fasta_records(request.input_fasta):
            result.total += 1
            progress.advance(task)
            record_id = header.split(" ", 1)[0]

            # Normalise en RNA
            seq_rna = sequence.strip().upper().replace("T", "U")
            # Supprime les caractères non-RNA
            seq_rna = "".join(c for c in seq_rna if c in "ACGUN")

            # Séquence trop longue → skip
            if len(seq_rna) > request.max_len:
                rec = StructureRecord(
                    record_id=record_id, sequence_rna=seq_rna,
                    dot_bracket="", mfe=0.0,
                    length=len(seq_rna), gc_frac=_gc_frac(seq_rna),
                    backend=backend, skipped=True, skip_reason="too_long",
                )
                result.skipped += 1
                jf.write(json.dumps({
                    "record_id": record_id, "skipped": True, "skip_reason": "too_long",
                    "length": len(seq_rna),
                }) + "\n")
                result.records.append(rec)
                continue

            try:
                dot_bracket, mfe = _predict(seq_rna, backend, request.temperature)
                rec = StructureRecord(
                    record_id=record_id, sequence_rna=seq_rna,
                    dot_bracket=dot_bracket, mfe=mfe,
                    length=len(seq_rna), gc_frac=_gc_frac(seq_rna),
                    backend=backend,
                )
                result.computed += 1
            except Exception as exc:
                logger.warning("Prédiction échouée pour %s : %s", record_id, exc)
                rec = StructureRecord(
                    record_id=record_id, sequence_rna=seq_rna,
                    dot_bracket="", mfe=0.0,
                    length=len(seq_rna), gc_frac=_gc_frac(seq_rna),
                    backend=backend, skipped=True, skip_reason=str(exc),
                )
                result.skipped += 1

            result.records.append(rec)
            jf.write(json.dumps({
                "record_id": rec.record_id,
                "length": rec.length,
                "gc_frac": round(rec.gc_frac, 4),
                "dot_bracket": rec.dot_bracket,
                "mfe": rec.mfe,
                "backend": rec.backend,
                "skipped": rec.skipped,
                "skip_reason": rec.skip_reason,
            }) + "\n")

            rows.append([
                rec.record_id,
                str(rec.length),
                f"{rec.gc_frac:.4f}",
                rec.dot_bracket,
                f"{rec.mfe:.2f}",
                rec.backend,
                str(rec.skipped),
                rec.skip_reason,
            ])

    write_csv(
        csv_path,
        ["record_id", "length", "gc_frac", "dot_bracket", "mfe_kcal_mol",
         "backend", "skipped", "skip_reason"],
        rows,
    )

    result.elapsed = time.perf_counter() - t0
    result.structures_csv = csv_path
    result.structures_jsonl = jsonl_path
    result.manifest_path = manifest_path

    manifest = {
        "timestamp": now_iso(),
        "backend": backend,
        "params": {
            "max_len": request.max_len,
            "temperature_celsius": request.temperature,
        },
        "inputs": {"fasta": str(request.input_fasta)},
        "total": result.total,
        "computed": result.computed,
        "skipped": result.skipped,
        "elapsed_s": round(result.elapsed, 2),
        "outputs": {
            "structures_csv": str(csv_path),
            "structures_jsonl": str(jsonl_path),
        },
    }
    write_json(manifest_path, manifest)

    logger.info(
        "Structure secondaire : %s séquences, %s calculées, %s ignorées, backend=%s, %.1fs",
        result.total, result.computed, result.skipped, backend, result.elapsed,
    )
    return result
