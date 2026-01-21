import json

import pytest

from arn_finder_v2.blast_parser import BlastParseError, parse_json2_top_hits


def test_parse_json2_top_hits_extracts_values():
    payload = {
        "BlastOutput2": [
            {
                "report": {
                    "results": {
                        "search": {
                            "query_len": 100,
                            "hits": [
                                {
                                    "description": [
                                        {
                                            "title": "Example virus",
                                            "organism": "Virus example",
                                            "taxid": 1234,
                                            "accession": "ABC123",
                                        }
                                    ],
                                    "accession": "ABC123",
                                    "hsps": [
                                        {
                                            "identity": 90,
                                            "align_len": 100,
                                            "evalue": 1e-10,
                                            "bit_score": 50.0,
                                            "query_from": 1,
                                            "query_to": 90,
                                        }
                                    ],
                                }
                            ],
                        }
                    }
                }
            }
        ]
    }
    hits = parse_json2_top_hits(json.dumps(payload), max_hits=5)
    assert hits[0]["accession"] == "ABC123"
    assert hits[0]["organism"] == "Virus example"
    assert hits[0]["pident"] == 90.0
    assert hits[0]["qcov"] == 90.0


def test_parse_json2_invalid_payload():
    with pytest.raises(BlastParseError):
        parse_json2_top_hits("{}", max_hits=1)
