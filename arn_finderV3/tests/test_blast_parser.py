"""Tests du parser BLAST JSON2 V3."""

from __future__ import annotations

import json

import pytest

from arn_finder_v3.blast.parser import parse_blast_json


def _make_blast_payload(accession: str, sciname: str, pident_frac: float, align_len: int, evalue: float) -> str:
    return json.dumps({
        "BlastOutput2": [{
            "report": {
                "results": {
                    "search": {
                        "hits": [{
                            "description": [{"accession": accession, "sciname": sciname, "title": f"Test [{sciname}]"}],
                            "len": align_len,
                            "hsps": [{
                                "align_len": align_len,
                                "identity": int(pident_frac * align_len),
                                "evalue": evalue,
                                "bit_score": 100.0,
                            }]
                        }]
                    }
                }
            }
        }]
    })


def test_parse_single_hit() -> None:
    raw = _make_blast_payload("NC_001", "Homo sapiens", 0.99, 300, 1e-50)
    hits = parse_blast_json(raw, max_hits=10, query_len=300)
    assert len(hits) == 1
    assert hits[0]["accession"] == "NC_001"
    assert hits[0]["organism"] == "Homo sapiens"
    assert abs(hits[0]["pident"] - 99.0) < 0.1
    assert hits[0]["coverage"] == 1.0


def test_max_hits_respected() -> None:
    payload = {
        "BlastOutput2": [{
            "report": {"results": {"search": {
                "hits": [
                    {
                        "description": [{"accession": f"AC{i}", "sciname": "Org", "title": "T"}],
                        "len": 100,
                        "hsps": [{"align_len": 100, "identity": 99, "evalue": 1e-10, "bit_score": 50.0}]
                    }
                    for i in range(20)
                ]
            }}}
        ]}
    }
    hits = parse_blast_json(json.dumps(payload), max_hits=5, query_len=100)
    assert len(hits) == 5


def test_empty_input_returns_empty() -> None:
    assert parse_blast_json("", max_hits=10, query_len=100) == []


def test_organism_fallback_from_brackets() -> None:
    payload = {
        "BlastOutput2": [{
            "report": {"results": {"search": {
                "hits": [{
                    "description": [{"accession": "X1", "title": "Sequence 1 [Mus musculus]"}],
                    "len": 100,
                    "hsps": [{"align_len": 100, "identity": 100, "evalue": 0.0, "bit_score": 200.0}]
                }]
            }}}
        ]}
    }
    hits = parse_blast_json(json.dumps(payload), max_hits=1, query_len=100)
    assert hits[0]["organism"] == "Mus musculus"
