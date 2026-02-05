import io
import re
from Bio import SeqIO


def validate_sequence(dna_seq: str, allow_n: bool = False) -> tuple[bool, str]:
    if dna_seq is None:
        return False, "No DNA sequence provided."

    seq = dna_seq.upper().replace("\n", "").replace(" ", "")
    if len(seq) < 23:
        return False, "Sequence too short. Minimum 23bp required."

    pattern = r"[ATCGN]+" if allow_n else r"[ATCG]+"
    if not re.fullmatch(pattern, seq):
        allowed = "A/T/C/G/N" if allow_n else "A/T/C/G"
        return False, f"Invalid characters in DNA sequence. Only {allowed} allowed."

    return True, seq


def load_fasta_text(contents: str) -> tuple[str, str]:
    fasta_io = io.StringIO(contents)
    records = list(SeqIO.parse(fasta_io, "fasta"))
    if records:
        if len(records) > 1:
            return str(records[0].seq), "Multiple FASTA records found. Using first sequence."
        return str(records[0].seq), ""

    seq = contents.strip().upper()
    if seq and re.fullmatch(r"[ATCGN\n\r ]+", seq):
        return seq.replace("\n", "").replace("\r", "").replace(" ", ""), ""

    return "", "No FASTA records found and content is not a valid plain DNA sequence."
