"""
Module codon_usage — V3.2

Calcule par séquence :
  - Fréquences de codons (compte brut + fréquence relative)
  - RSCU (Relative Synonymous Codon Usage, Sharp & Li 1986)
  - GC par position codon : GC1, GC2, GC3
  - ENc (Effective Number of Codons, Wright 1990) — indicateur de biais codon

Les séquences ADN sont automatiquement converties en ARN (T→U).
Seul le cadre de lecture +1 est utilisé (à partir de la position 0).
Les codons stop sont comptés mais exclus du RSCU et de l'ENc.

Sorties :
  codon_stats.csv   — par séquence : record_id, length, gc_frac, gc1, gc2, gc3,
                      enc, total_codons, stop_codons
  codon_rscu.csv    — long format  : record_id, codon, amino_acid, count, freq, rscu
  codon_manifest.json

Usage :
    from arn_finder_v3.codon_usage import CodonUsageRequest, compute_codon_usage
    result = compute_codon_usage(CodonUsageRequest(
        input_fasta=Path("data/filtered/filtered_sequences.fasta"),
        out_dir=Path("data/codon_usage"),
    ))
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from ..io import ensure_dir, get_logger, now_iso, read_fasta_records, write_csv, write_json

logger = get_logger("arnfinderv3.codon_usage")

# ── Table du code génétique standard ─────────────────────────────────────────

GENETIC_CODE: dict[str, str] = {
    "UUU": "Phe", "UUC": "Phe",
    "UUA": "Leu", "UUG": "Leu",
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile",
    "AUG": "Met",
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "UAU": "Tyr", "UAC": "Tyr",
    "UAA": "Stop", "UAG": "Stop",
    "CAU": "His", "CAC": "His",
    "CAA": "Gln", "CAG": "Gln",
    "AAU": "Asn", "AAC": "Asn",
    "AAA": "Lys", "AAG": "Lys",
    "GAU": "Asp", "GAC": "Asp",
    "GAA": "Glu", "GAG": "Glu",
    "UGU": "Cys", "UGC": "Cys",
    "UGA": "Stop",
    "UGG": "Trp",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGU": "Ser", "AGC": "Ser",
    "AGA": "Arg", "AGG": "Arg",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

# Codons regroupés par acide aminé (pour RSCU et ENc)
AA_CODON_GROUPS: dict[str, list[str]] = {
    "Phe": ["UUU", "UUC"],
    "Leu": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
    "Ile": ["AUU", "AUC", "AUA"],
    "Met": ["AUG"],
    "Val": ["GUU", "GUC", "GUA", "GUG"],
    "Ser": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
    "Pro": ["CCU", "CCC", "CCA", "CCG"],
    "Thr": ["ACU", "ACC", "ACA", "ACG"],
    "Ala": ["GCU", "GCC", "GCA", "GCG"],
    "Tyr": ["UAU", "UAC"],
    "Stop": ["UAA", "UAG", "UGA"],
    "His": ["CAU", "CAC"],
    "Gln": ["CAA", "CAG"],
    "Asn": ["AAU", "AAC"],
    "Lys": ["AAA", "AAG"],
    "Asp": ["GAU", "GAC"],
    "Glu": ["GAA", "GAG"],
    "Cys": ["UGU", "UGC"],
    "Trp": ["UGG"],
    "Arg": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "Gly": ["GGU", "GGC", "GGA", "GGG"],
}

# Classes de dégénérescence pour l'ENc (Wright 1990)
# 2 × 1-fold (Met, Trp) contribuent 2.0 constants à l'ENc
_ENc_2FOLD = ["Phe", "Tyr", "His", "Gln", "Asn", "Lys", "Asp", "Glu", "Cys"]  # 9 familles
_ENc_3FOLD = ["Ile"]                                                              # 1 famille
_ENc_4FOLD = ["Val", "Pro", "Thr", "Ala", "Gly"]                                # 5 familles
_ENc_6FOLD = ["Leu", "Ser", "Arg"]                                               # 3 familles


# ── Dataclasses ───────────────────────────────────────────────────────────────

@dataclass(slots=True)
class CodonUsageRequest:
    input_fasta: Path
    out_dir: Path
    reading_frame: int = 0            # 0, 1 ou 2 (offset de départ)
    min_codons: int = 10              # séquences avec moins de codons → ignorées


@dataclass
class CodonUsageResult:
    stats_csv: Path
    rscu_csv: Path
    manifest_path: Path
    total_records: int
    computed: int
    skipped: int
    elapsed: float


# ── Calculs ───────────────────────────────────────────────────────────────────

def _dna_to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def _codon_counts(rna: str, offset: int) -> dict[str, int]:
    """Compte les triplets valides dans la séquence ARN."""
    counts: dict[str, int] = {codon: 0 for codon in GENETIC_CODE}
    for i in range(offset, len(rna) - 2, 3):
        triplet = rna[i:i+3]
        if triplet in counts:
            counts[triplet] += 1
    return counts


def _gc_frac(seq: str) -> float:
    if not seq:
        return 0.0
    return round(sum(1 for b in seq.upper() if b in "GC") / len(seq), 4)


def _gc_by_position(codon_list: list[str]) -> tuple[float, float, float]:
    """GC aux positions 1, 2, 3 du codon (indices 0, 1, 2)."""
    if not codon_list:
        return 0.0, 0.0, 0.0
    n = len(codon_list)
    gc1 = round(sum(1 for c in codon_list if c[0] in "GC") / n, 4)
    gc2 = round(sum(1 for c in codon_list if c[1] in "GC") / n, 4)
    gc3 = round(sum(1 for c in codon_list if c[2] in "GC") / n, 4)
    return gc1, gc2, gc3


def _rscu(counts: dict[str, int]) -> dict[str, float]:
    """
    RSCU par codon (Sharp & Li 1986).
    RSCU = 1.0 → usage neutre ; >1.0 → codon préféré ; <1.0 → codon évité.
    Les codons stop sont exclus (RSCU = 0).
    """
    rscu: dict[str, float] = {}
    for aa, codons in AA_CODON_GROUPS.items():
        if aa == "Stop":
            for c in codons:
                rscu[c] = 0.0
            continue
        n_syn = len(codons)
        total = sum(counts.get(c, 0) for c in codons)
        expected = total / n_syn if n_syn > 0 else 0.0
        for c in codons:
            obs = counts.get(c, 0)
            rscu[c] = round(obs / expected, 4) if expected > 0 else 0.0
    return rscu


def _f_stat(codon_counts_for_aa: list[int]) -> float:
    """Statistique F de Wright pour une famille d'acides aminés."""
    ni = sum(codon_counts_for_aa)
    n_syn = len(codon_counts_for_aa)
    if ni <= 1 or n_syn <= 1:
        return 1.0
    sum_pi_sq = sum((c / ni) ** 2 for c in codon_counts_for_aa)
    fc = (ni * sum_pi_sq - 1.0) / (ni - 1.0)
    # Clamp dans l'intervalle valide [1/n_syn, 1.0]
    return max(1.0 / n_syn, min(1.0, fc))


def _f_mean_for_class(aa_list: list[str], counts: dict[str, int]) -> float:
    """Moyenne pondérée de F pour une classe de dégénérescence."""
    total_n = 0
    weighted_sum = 0.0
    for aa in aa_list:
        codons = AA_CODON_GROUPS.get(aa, [])
        aa_counts = [counts.get(c, 0) for c in codons]
        ni = sum(aa_counts)
        if ni == 0:
            continue
        fc = _f_stat(aa_counts)
        weighted_sum += ni * fc
        total_n += ni
    return weighted_sum / total_n if total_n > 0 else 1.0


def _enc(counts: dict[str, int]) -> float:
    """
    ENc (Effective Number of Codons) — Wright (1990).
    Range : 20 (biais maximal, un seul codon par AA) → 61 (usage non biaisé).
    Formule : ENc = 2 + 9/F̂2 + 1/F̂3 + 5/F̂4 + 3/F̂6
    """
    f2 = _f_mean_for_class(_ENc_2FOLD, counts)
    f3 = _f_mean_for_class(_ENc_3FOLD, counts)
    f4 = _f_mean_for_class(_ENc_4FOLD, counts)
    f6 = _f_mean_for_class(_ENc_6FOLD, counts)

    enc = 2.0  # Met + Trp : 1-fold, contribuent toujours 1.0 chacun
    enc += 9.0 / f2
    enc += 1.0 / f3
    enc += 5.0 / f4
    enc += 3.0 / f6

    return round(max(20.0, min(61.0, enc)), 2)


# ── Orchestrateur ─────────────────────────────────────────────────────────────

def compute_codon_usage(request: CodonUsageRequest) -> CodonUsageResult:
    """
    Calcule les fréquences de codons, RSCU, GC1/2/3 et ENc pour chaque séquence.
    """
    out_dir = ensure_dir(request.out_dir)
    stats_path = out_dir / "codon_stats.csv"
    rscu_path = out_dir / "codon_rscu.csv"
    manifest_path = out_dir / "codon_manifest.json"

    t0 = time.perf_counter()
    stats_rows: list[list] = []
    rscu_rows: list[list] = []
    total = computed = skipped = 0

    try:
        from rich.progress import Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
        _rich = True
    except ImportError:
        _rich = False

    def _process() -> None:
        nonlocal total, computed, skipped
        for header, sequence in read_fasta_records(request.input_fasta):
            total += 1
            record_id = header.split(" ", 1)[0]
            rna = _dna_to_rna(sequence)
            counts = _codon_counts(rna, request.reading_frame)

            total_codons = sum(counts.values())
            stop_codons = sum(counts.get(c, 0) for c in AA_CODON_GROUPS["Stop"])
            coding_codons = total_codons - stop_codons

            if coding_codons < request.min_codons:
                skipped += 1
                logger.debug("Ignoré %s : trop peu de codons (%d)", record_id, coding_codons)
                continue

            # GC global
            gc = _gc_frac(rna)

            # GC par position (sur codons coding uniquement)
            coding_list = [
                rna[i:i+3]
                for i in range(request.reading_frame, len(rna) - 2, 3)
                if rna[i:i+3] in GENETIC_CODE and GENETIC_CODE[rna[i:i+3]] != "Stop"
            ]
            gc1, gc2, gc3 = _gc_by_position(coding_list) if coding_list else (0.0, 0.0, 0.0)

            # ENc
            enc_val = _enc(counts)

            stats_rows.append([
                record_id, len(sequence), gc, gc1, gc2, gc3,
                enc_val, total_codons, stop_codons,
            ])

            # RSCU
            rscu_vals = _rscu(counts)
            for codon, r_val in rscu_vals.items():
                aa = GENETIC_CODE.get(codon, "?")
                count = counts.get(codon, 0)
                freq = round(count / total_codons, 6) if total_codons else 0.0
                rscu_rows.append([record_id, codon, aa, count, freq, r_val])

            computed += 1

    if _rich:
        with Progress(
            SpinnerColumn(), TextColumn("{task.description}"), TimeElapsedColumn()
        ) as prog:
            task = prog.add_task("codon_usage", total=None)
            _process()
            prog.update(task, description=f"codon_usage : {computed} séquences")
    else:
        _process()

    stats_headers = [
        "record_id", "length", "gc_frac", "gc1", "gc2", "gc3",
        "enc", "total_codons", "stop_codons",
    ]
    write_csv(stats_path, stats_headers, stats_rows)

    rscu_headers = ["record_id", "codon", "amino_acid", "count", "freq", "rscu"]
    write_csv(rscu_path, rscu_headers, rscu_rows)

    elapsed = time.perf_counter() - t0
    manifest = {
        "timestamp": now_iso(),
        "params": {
            "reading_frame": request.reading_frame,
            "min_codons": request.min_codons,
        },
        "inputs": {"fasta": str(request.input_fasta)},
        "total_records": total,
        "computed": computed,
        "skipped": skipped,
        "elapsed_s": round(elapsed, 1),
        "outputs": {
            "stats_csv": str(stats_path),
            "rscu_csv": str(rscu_path),
        },
    }
    write_json(manifest_path, manifest)
    logger.info(
        "Codon usage terminé : %d/%d séquences calculées en %.1fs",
        computed, total, elapsed,
    )

    return CodonUsageResult(
        stats_csv=stats_path,
        rscu_csv=rscu_path,
        manifest_path=manifest_path,
        total_records=total,
        computed=computed,
        skipped=skipped,
        elapsed=elapsed,
    )
