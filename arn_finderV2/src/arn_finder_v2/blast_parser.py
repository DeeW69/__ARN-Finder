"""Helpers to interpret BLAST JSON responses."""

from __future__ import annotations

import json
from typing import Any, Dict, List, Optional


class BlastParseError(RuntimeError):
    """Raised when BLAST results cannot be parsed."""


def parse_json2_top_hits(raw_json: str, max_hits: int) -> List[Dict[str, Optional[float]]]:
    """Extract top hits from a BLAST JSON2 payload."""
    try:
        payload = json.loads(raw_json)
    except json.JSONDecodeError as exc:
        raise BlastParseError("Invalid JSON2 payload") from exc

    try:
        search = payload["BlastOutput2"][0]["report"]["results"]["search"]
        hits = search.get("hits", [])
        query_len = search.get("query_len")
    except (KeyError, IndexError, TypeError) as exc:
        raise BlastParseError("Unexpected JSON2 schema") from exc

    parsed: List[Dict[str, Any]] = []
    for hit in hits:
        description = (hit.get("description") or [{}])[0]
        organism = description.get("organism") or description.get("sciname") or description.get("taxname")
        taxid = description.get("taxid")
        accession = hit.get("accession") or description.get("accession")
        title = description.get("title")

        hsp = (hit.get("hsps") or [{}])[0]
        align_len = hsp.get("align_len")
        identity = hsp.get("identity")
        pident = _compute_percentage(identity, align_len)
        evalue = hsp.get("evalue")
        bitscore = hsp.get("bit_score")
        qcov = _compute_query_coverage(hsp, query_len)

        parsed.append(
            {
                "accession": accession,
                "title": title,
                "taxid": taxid,
                "organism": organism,
                "pident": pident,
                "align_len": align_len,
                "evalue": evalue,
                "bitscore": bitscore,
                "qcov": qcov,
            }
        )
        if len(parsed) >= max_hits:
            break

    return parsed


def _compute_percentage(identity: Optional[int], align_len: Optional[int]) -> Optional[float]:
    if identity is None or not align_len:
        return None
    return round(100 * identity / align_len, 2)


def _compute_query_coverage(hsp: Dict[str, Any], query_len: Optional[int]) -> Optional[float]:
    if not query_len:
        return None
    q_from = hsp.get("query_from")
    q_to = hsp.get("query_to")
    if q_from is None or q_to is None:
        return None
    covered = abs(q_to - q_from) + 1
    return round(100 * covered / query_len, 2)
