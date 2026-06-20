"""Genome-wide off-target search.

Extends the off-target engine from "scan a pasted background" to "scan a whole
genome FASTA" — multi-record (chromosomes/contigs), both strands, with global
coordinates and the same CFD + MIT/Hsu + aggregate-specificity scoring.

Memory safety: each record is scanned in overlapping chunks, so peak memory is
bounded by ``chunk_size`` rather than the genome size. The NumPy
sliding-window mismatch count makes each chunk fast. This handles bacterial,
viral, plasmid, organellar and small-eukaryotic genomes (and any user FASTA)
comfortably; mammalian whole-genome scans run but are slower (a seed/FM-index
backend is the documented next step).

Off-target semantics match the in-memory engine: perfect (0-mismatch) matches
are treated as the intended target and excluded; reported hits have 1..max
mismatches with a valid PAM.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd

try:
    from offtarget import aggregate_specificity, calculate_cfd_score, mit_hit_score
except ImportError:  # pragma: no cover
    from .offtarget import aggregate_specificity, calculate_cfd_score, mit_hit_score

_PAM_LEN = {"NGG": 3, "NAG": 3, "NG": 2, "TTTV": 4}


def _check_pam(pam_seq: str, pam: str) -> bool:
    if pam == "NGG":
        return len(pam_seq) == 3 and pam_seq[1:] == "GG"
    if pam == "NAG":
        return len(pam_seq) == 3 and pam_seq[1:] == "AG"
    if pam == "NG":
        return len(pam_seq) == 2 and pam_seq[1] == "G"
    if pam == "TTTV":
        return len(pam_seq) == 4 and pam_seq[:3] == "TTT" and pam_seq[3] in "ACG"
    return False


_COMP = str.maketrans("ACGTN", "TGCAN")


def _revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def iter_fasta(path: str) -> Iterator[tuple[str, str]]:
    """Stream (record_id, sequence) pairs from a FASTA file, one record at a time."""
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0] if len(line) > 1 else "seq"
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks).upper()


def _scan_strand(seq: str, strand: str, rec_len: int, guide: str, guide_arr: np.ndarray,
                 guide_len: int, pam: str, pam_len: int, max_mismatches: int,
                 chunk_size: int) -> list[dict]:
    """Scan one strand of one record in overlapping chunks (bounded memory)."""
    rows: list[dict] = []
    n = len(seq)
    overlap = guide_len + pam_len
    start = 0
    while start <= n - guide_len - pam_len:
        end = min(n, start + chunk_size + overlap)
        sub = seq[start:end]
        arr = np.frombuffer(sub.encode("ascii", "replace"), dtype=np.int8)
        usable = len(arr) - pam_len
        if usable > guide_len:
            windows = np.lib.stride_tricks.sliding_window_view(arr[:usable], guide_len)
            mm = np.sum(windows != guide_arr, axis=1)
            # Only accept window starts owned by this chunk to avoid double counting.
            limit = min(len(mm), chunk_size)
            cand = np.where((mm[:limit] >= 1) & (mm[:limit] <= max_mismatches))[0]
            for off in cand:
                i = int(off)
                pam_seq = sub[i + guide_len:i + guide_len + pam_len]
                if _check_pam(pam_seq, pam):
                    window_seq = sub[i:i + guide_len]
                    gpos = start + i  # position on the scanned strand
                    pos = gpos if strand == "+" else rec_len - gpos - guide_len - pam_len
                    rows.append({
                        "gRNA": guide,
                        "Strand": strand,
                        "Pos": int(pos),
                        "Mismatches": int(mm[i]),
                        "TargetSeq": window_seq,
                        "PAM": pam_seq,
                        "CFD_Score": calculate_cfd_score(guide, window_seq, pam_seq),
                        "MIT_Score": mit_hit_score(guide, window_seq),
                    })
        start += chunk_size
    return rows


def scan_records(records: Iterator[tuple[str, str]], guides: list[str],
                 max_mismatches: int = 3, pam: str = "NGG",
                 chunk_size: int = 2_000_000) -> pd.DataFrame:
    """Scan an iterator of (record_id, sequence) for off-targets of each guide."""
    pam_len = _PAM_LEN.get(pam, 3)
    guide_arrs = {g: np.frombuffer(g.encode(), dtype=np.int8) for g in guides}
    rows: list[dict] = []
    for rec_id, seq in records:
        rc = _revcomp(seq)
        rec_len = len(seq)
        for g, garr in guide_arrs.items():
            gl = len(g)
            for strand, s in (("+", seq), ("-", rc)):
                hits = _scan_strand(s, strand, rec_len, g, garr, gl, pam, pam_len,
                                    max_mismatches, chunk_size)
                for h in hits:
                    h["Chrom"] = rec_id
                rows.extend(hits)
    cols = ["gRNA", "Chrom", "Strand", "Pos", "Mismatches", "TargetSeq", "PAM", "CFD_Score", "MIT_Score"]
    return pd.DataFrame(rows, columns=cols)


def scan_fasta_file(path: str, guides: list[str], max_mismatches: int = 3,
                    pam: str = "NGG", chunk_size: int = 2_000_000) -> pd.DataFrame:
    return scan_records(iter_fasta(path), guides, max_mismatches, pam, chunk_size)


def scan_fasta_text(text: str, guides: list[str], max_mismatches: int = 3,
                    pam: str = "NGG", chunk_size: int = 2_000_000) -> pd.DataFrame:
    def _records():
        if text.lstrip().startswith(">"):
            name, chunks = None, []
            for line in text.splitlines():
                if line.startswith(">"):
                    if name is not None:
                        yield name, "".join(chunks).upper()
                    name = line[1:].split()[0] if len(line) > 1 else "seq"
                    chunks = []
                else:
                    chunks.append(line.strip())
            if name is not None:
                yield name, "".join(chunks).upper()
        else:
            yield "input", "".join(text.split()).upper()
    return scan_records(_records(), guides, max_mismatches, pam, chunk_size)


def genome_specificity(hits: pd.DataFrame, guides: list[str]) -> pd.DataFrame:
    """Per-guide genome-wide aggregate specificity (CRISPOR convention)."""
    rows = []
    for g in guides:
        sub = hits[hits["gRNA"] == g] if not hits.empty else hits
        cfd = list(sub["CFD_Score"]) if not sub.empty else []
        mit = list(sub["MIT_Score"]) if not sub.empty else []
        rows.append({
            "gRNA": g,
            "OffTargetCount": int(len(sub)),
            "CFD_Specificity": aggregate_specificity(cfd, scale=100.0),
            "MIT_Specificity": aggregate_specificity(mit, scale=1.0),
        })
    return pd.DataFrame(rows)


def _cli() -> None:
    import argparse
    ap = argparse.ArgumentParser(description="Genome-wide CRISPR off-target search")
    ap.add_argument("fasta", help="genome FASTA file")
    ap.add_argument("guides", nargs="+", help="20-nt guide sequences")
    ap.add_argument("--max-mismatches", type=int, default=3)
    ap.add_argument("--pam", default="NGG")
    args = ap.parse_args()

    hits = scan_fasta_file(args.fasta, [g.upper() for g in args.guides],
                           args.max_mismatches, args.pam)
    spec = genome_specificity(hits, [g.upper() for g in args.guides])
    print(spec.to_string(index=False))
    print(f"\n{len(hits)} off-target site(s) found. Top by CFD:")
    if not hits.empty:
        print(hits.sort_values("CFD_Score", ascending=False).head(20).to_string(index=False))


if __name__ == "__main__":
    sys.exit(_cli())
