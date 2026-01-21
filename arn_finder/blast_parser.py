"""Parsing helpers for BLAST JSON2 payloads."""

from __future__ import annotations

import json
from typing import Any, Dict, List


def _extract_number(container: Dict[str, Any], keys: List[str], default: float = 0.0) -> float:
    for key in keys:
        if key in container and container[key] is not None:
            try:
                return float(container[key])
            except (TypeError, ValueError):
                continue
    return default


def _as_int(value: Any) -> int:
    try:
        return int(round(float(value)))
    except (TypeError, ValueError):
        return 0


def _extract_organism(description: Dict[str, Any], title: str) -> str:
    organism = description.get("sciname")
    if organism:
        return str(organism)
    if "[" in title and "]" in title:
        # Use the text in the last brackets as organism shorthand.
        try:
            return title.split("[")[-1].rstrip("]")
        except Exception:
            pass
    taxid = description.get("taxid")
    return str(taxid) if taxid is not None else ""


def parse_blast_json(raw_json: str, max_hits: int, query_len: int) -> List[Dict[str, Any]]:
    """Extract simplified hit summaries from a BLAST JSON2 payload."""

    if not raw_json.strip():
        return []
    payload = json.loads(raw_json)
    reports = payload.get("BlastOutput2", [])
    hits_out: List[Dict[str, Any]] = []
    for report in reports:
        search = (
            report.get("report", {})
            .get("results", {})
            .get("search", {})
        )
        hits = search.get("hits", [])
        for hit in hits:
            descriptions = hit.get("description") or [{}]
            desc = descriptions[0]
            title = desc.get("title") or desc.get("id") or ""
            hsp_list = hit.get("hsps") or []
            if not hsp_list:
                continue
            hsp = hsp_list[0]
            align_len = _extract_number(hsp, ["align_len", "align-len", "alignLen"], default=0.0)
            identity = _extract_number(
                hsp,
                ["identity", "identity_num", "identity-val", "num_identical"],
                default=0.0,
            )
            length = hit.get("len") or align_len
            pident = (identity / align_len * 100.0) if align_len else 0.0
            coverage = (align_len / query_len) if query_len and align_len else 0.0
            hits_out.append(
                {
                    "accession": desc.get("accession") or desc.get("id") or "",
                    "title": title,
                    "organism": _extract_organism(desc, title),
                    "pident": round(pident, 3),
                    "length": _as_int(length),
                    "evalue": _extract_number(hsp, ["evalue", "e-value", "expect"], default=0.0),
                    "bitscore": _extract_number(hsp, ["bit_score", "bitScore", "bit-score"], default=0.0),
                    "coverage": round(coverage, 3),
                }
            )
            if len(hits_out) >= max_hits:
                return hits_out
    return hits_out


__all__ = ["parse_blast_json"]
