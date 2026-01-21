from pathlib import Path

from arn_finder_v2.fasta import parse_fasta_lines, read_fasta


def test_parse_fasta_lines_merges_sequences():
    fasta = [
        ">seq1 description here",
        "acgt",
        "tgca",
        ">seq2",
        "UUUU",
    ]
    records = list(parse_fasta_lines(fasta, uppercase=True, convert_u_to_t=True))
    assert records[0].identifier == "seq1"
    assert records[0].sequence == "ACGTTGCA"
    assert records[1].sequence == "TTTT"


def test_read_fasta_from_file(tmp_path: Path):
    content = ">a\nAA\n> b c\ncc"
    in_path = tmp_path / "in.fasta"
    in_path.write_text(content, encoding="utf-8")
    records = read_fasta(in_path)
    assert len(records) == 2
    assert records[1].identifier == "b"
    assert records[1].description == "c"
