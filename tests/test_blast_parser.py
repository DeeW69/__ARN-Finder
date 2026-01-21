"""Smoke tests for the BLAST JSON parser."""

from __future__ import annotations

import json

from arn_finder.blast_parser import parse_blast_json


def test_parse_blast_json_extracts_hits() -> None:
    payload = {
        "BlastOutput2": [
            {
                "report": {
                    "results": {
                        "search": {
                            "hits": [
                                {
                                    "description": [
                                        {"accession": "ABC123", "title": "Example hit [Test organism]"},
                                    ],
                                    "len": 250,
                                    "hsps": [
                                        {
                                            "align_len": 250,
                                            "identity": 240,
                                            "bit_score": 310.5,
                                            "evalue": 5e-20,
                                        }
                                    ],
                                },
                                {
                                    "description": [
                                        {"accession": "XYZ789", "title": "Second hit"},
                                    ],
                                    "len": 220,
                                    "hsps": [
                                        {
                                            "alignLen": 200,
                                            "identity": 180,
                                            "bit_score": 250,
                                            "evalue": 1e-5,
                                        }
                                    ],
                                },
                            ]
                        }
                    }
                }
            }
        ]
    }
    hits = parse_blast_json(json.dumps(payload), max_hits=2, query_len=500)
    assert len(hits) == 2
    first = hits[0]
    assert first["accession"] == "ABC123"
    assert first["organism"] == "Test organism"
    assert first["length"] == 250
    assert first["bitscore"] == 310.5
    assert first["pident"] == 96.0
    assert first["coverage"] == 0.5


def test_parse_blast_json_limits_hits() -> None:
    payload = {
        "BlastOutput2": [
            {
                "report": {
                    "results": {
                        "search": {
                            "hits": [
                                {
                                    "description": [{"accession": "A1", "title": "hit1"}],
                                    "hsps": [{"align_len": 100, "identity": 80, "bit_score": 200, "evalue": 1e-6}],
                                },
                                {
                                    "description": [{"accession": "A2", "title": "hit2"}],
                                    "hsps": [{"align_len": 100, "identity": 70, "bit_score": 190, "evalue": 1e-5}],
                                },
                            ]
                        }
                    }
                }
            }
        ]
    }
    hits = parse_blast_json(json.dumps(payload), max_hits=1, query_len=100)
    assert len(hits) == 1
    assert hits[0]["accession"] == "A1"
