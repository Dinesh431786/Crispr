from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from Bio.Seq import Seq
import pandas as pd


@dataclass(frozen=True)
class ScoringConfig:
    """Research-informed sequence constraints for CRISPR guide ranking."""

    gc_low: float = 0.40
    gc_high: float = 0.65
    seed_window: int = 8
    homopolymer_penalty: float = 0.15
    off_target_penalty: float = 0.06


def _clean_seq(seq: str) -> str:
    return seq.upper().replace("\n", "").replace(" ", "")


def _gc_fraction(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / max(len(seq), 1)


def hybrid_score(guide: str, off_target_count: int = 0, cfg: ScoringConfig = ScoringConfig()) -> float:
    guide = guide.upper()
    gc = _gc_fraction(guide)
    score = 0.55

    if cfg.gc_low <= gc <= cfg.gc_high:
        score += 0.20
    else:
        score -= 0.12

    seed = guide[-cfg.seed_window:]
    if seed.count("G") + seed.count("C") >= cfg.seed_window // 2:
        score += 0.10

    if any(base * 4 in guide for base in "ATCG"):
        score -= cfg.homopolymer_penalty

    if guide.startswith("G"):
        score += 0.04
    if guide.endswith("GG"):
        score += 0.06

    score -= cfg.off_target_penalty * max(0, off_target_count)
    return round(max(0.0, min(score, 1.0)), 3)


def ml_gRNA_score(guide: str) -> float:
    guide = guide.upper()
    gc = _gc_fraction(guide)
    score = 0.50

    if 0.45 <= gc <= 0.60:
        score += 0.18
    elif 0.35 <= gc < 0.45 or 0.60 < gc <= 0.70:
        score += 0.08

    seed = guide[-10:]
    if seed.count("G") + seed.count("C") >= 5:
        score += 0.10
    if guide[0] in {"G", "C"}:
        score += 0.04
    if len(guide) > 16 and guide[16] in {"A", "T"}:
        score += 0.03

    if "TTTT" in guide:
        score -= 0.15
    if any(motif in guide for motif in ("AAAA", "CCCC", "GGGG")):
        score -= 0.10

    return round(max(0.0, min(score, 1.0)), 3)


def check_pam(pam_seq: str, pam: str) -> bool:
    pam_seq = pam_seq.upper()
    if pam == "NGG":
        return len(pam_seq) == 3 and pam_seq[1:] == "GG"
    if pam == "NAG":
        return len(pam_seq) == 3 and pam_seq[1:] == "AG"
    if pam == "NG":
        return len(pam_seq) == 2 and pam_seq[1] == "G"
    if pam == "TTTV":
        return len(pam_seq) == 4 and pam_seq[:3] == "TTT" and pam_seq[3] in "ACG"
    return False


def find_gRNAs(
    dna_seq: str,
    pam: str = "NGG",
    guide_length: int = 20,
    min_gc: int = 40,
    max_gc: int = 70,
    add_5prime_g: bool = False,
) -> pd.DataFrame:
    sequence = _clean_seq(dna_seq)
    pam_len = 2 if pam == "NG" else (4 if pam == "TTTV" else 3)
    guides: list[dict[str, str | float | int]] = []

    def _append_guide(strand: str, start: int, guide: str, pam_seq: str) -> None:
        gc = _gc_fraction(guide) * 100
        if min_gc <= gc <= max_gc and "TTTT" not in guide:
            g_out = guide if not add_5prime_g or guide.startswith("G") else f"G{guide[:-1]}"
            ml = ml_gRNA_score(g_out)
            hy = hybrid_score(g_out)
            guides.append(
                {
                    "Strand": strand,
                    "Start": start,
                    "gRNA": g_out,
                    "PAM": pam_seq,
                    "GC%": round(gc, 2),
                    "HybridScore": hy,
                    "MLScore": ml,
                    "ConsensusScore": round((hy + ml) / 2, 3),
                }
            )

    for i in range(len(sequence) - guide_length - pam_len + 1):
        guide = sequence[i : i + guide_length]
        pam_seq = sequence[i + guide_length : i + guide_length + pam_len]
        if check_pam(pam_seq, pam):
            _append_guide("+", i, guide, pam_seq)

    rc_sequence = str(Seq(sequence).reverse_complement())
    for i in range(len(rc_sequence) - guide_length - pam_len + 1):
        guide = rc_sequence[i : i + guide_length]
        pam_seq = rc_sequence[i + guide_length : i + guide_length + pam_len]
        if check_pam(pam_seq, pam):
            _append_guide("-", len(sequence) - i - guide_length - pam_len, guide, pam_seq)

    if not guides:
        return pd.DataFrame(columns=["Strand", "Start", "gRNA", "PAM", "GC%", "HybridScore", "MLScore", "ConsensusScore"])

    return pd.DataFrame(guides).sort_values("ConsensusScore", ascending=False).reset_index(drop=True)


def count_mismatches(a: str, b: str) -> int:
    a = a.upper()
    b = b.upper()
    if len(a) != len(b):
        return 10**9
    return sum(1 for x, y in zip(a, b) if x != y)


CFD_SCORES = {
    "rA:dA": 1.0,
    "rA:dC": 1.0,
    "rA:dG": 0.857,
    "rA:dT": 1.0,
    "rC:dA": 1.0,
    "rC:dC": 0.913,
    "rC:dG": 1.0,
    "rC:dT": 1.0,
    "rG:dA": 1.0,
    "rG:dC": 1.0,
    "rG:dG": 0.714,
    "rG:dT": 0.9,
    "rU:dA": 1.0,
    "rU:dC": 0.957,
    "rU:dG": 0.857,
    "rU:dT": 1.0,
}
PAM_SCORES = {
    "AA": 0.0,
    "AC": 0.0,
    "AG": 0.259,
    "AT": 0.0,
    "CA": 0.0,
    "CC": 0.0,
    "CG": 0.107,
    "CT": 0.0,
    "GA": 0.069,
    "GC": 0.022,
    "GG": 1.0,
    "GT": 0.016,
    "TA": 0.0,
    "TC": 0.0,
    "TG": 0.039,
    "TT": 0.0,
}


def calculate_cfd_score(guide_seq: str, off_target_seq: str, pam: str) -> float:
    score = 1.0
    guide_seq = guide_seq.upper().replace("T", "U")
    off_target_seq = off_target_seq.upper()

    for i in range(min(len(guide_seq), len(off_target_seq))):
        if guide_seq[i] != off_target_seq[i]:
            key = f"r{guide_seq[i]}:d{off_target_seq[i]}"
            score *= CFD_SCORES.get(key, 1.0)

    pam_key = pam[1:] if len(pam) >= 2 else pam
    score *= PAM_SCORES.get(pam_key, 0.0)
    return round(score, 5)


def _segment_boundaries(length: int, max_mismatches: int) -> list[tuple[int, int]]:
    parts = max(1, max_mismatches + 1)
    seg_size = length // parts
    boundaries = []
    start = 0
    for i in range(parts):
        end = start + seg_size
        if i == parts - 1:
            end = length
        boundaries.append((start, end))
        start = end
    return boundaries


def _build_seed_indexes(bg_seq: str, guide_len: int, max_mismatches: int) -> list[dict[str, list[int]]]:
    boundaries = _segment_boundaries(guide_len, max_mismatches)
    indexes: list[dict[str, list[int]]] = []

    for seg_start, seg_end in boundaries:
        segment_index: dict[str, list[int]] = {}
        seg_len = seg_end - seg_start
        for pos in range(0, len(bg_seq) - guide_len + 1):
            segment = bg_seq[pos + seg_start : pos + seg_start + seg_len]
            if len(segment) != seg_len:
                continue
            segment_index.setdefault(segment, []).append(pos)
        indexes.append(segment_index)

    return indexes


def _candidate_positions_from_index(
    guide_seq: str,
    seed_indexes: list[dict[str, list[int]]],
    max_mismatches: int,
) -> Iterable[int]:
    boundaries = _segment_boundaries(len(guide_seq), max_mismatches)
    seen: set[int] = set()
    for (start, end), seg_index in zip(boundaries, seed_indexes):
        key = guide_seq[start:end]
        for pos in seg_index.get(key, []):
            if pos not in seen:
                seen.add(pos)
                yield pos


def _scan_off_targets_single_strand(
    guides: pd.DataFrame,
    bg_seq: str,
    max_mismatches: int,
    pam_len: int,
    strand_label: str,
) -> list[dict[str, int | str | float]]:
    results: list[dict[str, int | str | float]] = []
    if guides.empty:
        return results

    guide_len = len(guides["gRNA"].iloc[0])
    search_seq = bg_seq[: max(0, len(bg_seq) - pam_len)]
    seed_indexes = _build_seed_indexes(search_seq, guide_len, max_mismatches)

    for guide_seq in guides["gRNA"].unique():
        for pos in _candidate_positions_from_index(guide_seq, seed_indexes, max_mismatches):
            window = bg_seq[pos : pos + guide_len]
            mismatches = count_mismatches(guide_seq, window)
            if 0 < mismatches <= max_mismatches:
                pam_seq = bg_seq[pos + guide_len : pos + guide_len + pam_len]
                if len(pam_seq) != pam_len:
                    continue
                results.append(
                    {
                        "gRNA": guide_seq,
                        "OffTargetPos": pos,
                        "Mismatches": mismatches,
                        "TargetSeq": window,
                        "PAM": pam_seq,
                        "CFD_Score": calculate_cfd_score(guide_seq, window, pam_seq),
                        "Strand": strand_label,
                    }
                )
    return results


def find_off_targets_detailed(
    guides: pd.DataFrame,
    background_seq: str,
    max_mismatches: int = 2,
    pam: str = "NGG",
    include_reverse_strand: bool = True,
) -> pd.DataFrame:
    """Find potential off-targets using indexed candidate generation.

    Uses a segment-index strategy (pigeonhole principle) to avoid brute-force
    comparisons at every genomic position while preserving sensitivity for
    mismatch-bounded search.
    """
    if guides is None or guides.empty:
        return pd.DataFrame(columns=["gRNA", "OffTargetPos", "Mismatches", "TargetSeq", "PAM", "CFD_Score", "Strand"])

    bg_seq = _clean_seq(background_seq)[:1_000_000]
    pam_len = 3 if pam in ["NGG", "NAG"] else (2 if pam == "NG" else 4)

    results = _scan_off_targets_single_strand(guides, bg_seq, max_mismatches, pam_len, "+")

    if include_reverse_strand:
        rc_seq = str(Seq(bg_seq).reverse_complement())
        rc_results = _scan_off_targets_single_strand(guides, rc_seq, max_mismatches, pam_len, "-")
        for row in rc_results:
            # map reverse complement coordinate back to original background coordinates
            row["OffTargetPos"] = len(bg_seq) - int(row["OffTargetPos"]) - len(str(row["TargetSeq"])) - pam_len
        results.extend(rc_results)

    if not results:
        return pd.DataFrame(columns=["gRNA", "OffTargetPos", "Mismatches", "TargetSeq", "PAM", "CFD_Score", "Strand"])

    return (
        pd.DataFrame(results)
        .drop_duplicates(subset=["gRNA", "OffTargetPos", "TargetSeq", "PAM", "Strand"])
        .sort_values(["CFD_Score", "Mismatches"], ascending=[False, True])
        .reset_index(drop=True)
    )


def safe_translate(seq: str) -> str:
    extra = len(seq) % 3
    if extra:
        seq += "N" * (3 - extra)
    try:
        return str(Seq(seq).translate(to_stop=True))
    except Exception:
        return "Translation failed"


def simulate_protein_edit(
    seq: str,
    cut_index: int,
    edit_type: str = "del1",
    insert_base: str = "A",
    sub_from: str = "A",
    sub_to: str = "T",
) -> tuple[str, str, bool, bool]:
    seq = _clean_seq(seq)
    before = seq
    after = seq

    if edit_type == "del1":
        after = seq[:cut_index] + seq[cut_index + 1 :]
    elif edit_type == "insA":
        after = seq[:cut_index] + insert_base + seq[cut_index:]
    elif edit_type.startswith("del"):
        del_len = int(edit_type[3:])
        after = seq[:cut_index] + seq[cut_index + del_len :]
    elif edit_type.startswith("ins"):
        ins = edit_type[3:] or insert_base
        after = seq[:cut_index] + ins + seq[cut_index:]
    elif edit_type == "subAG" and cut_index + len(sub_from) <= len(seq):
        after = seq[:cut_index] + sub_to + seq[cut_index + len(sub_from) :]

    prot_before = safe_translate(before)
    prot_after = safe_translate(after)
    frameshift = (len(after) % 3) != (len(before) % 3)
    stop_lost = "*" in prot_after and "*" not in prot_before
    return prot_before, prot_after, frameshift, stop_lost


def indel_simulations(seq: str, cut_index: int) -> pd.DataFrame:
    seq = _clean_seq(seq)
    rows = []
    for n in (1, 2, 3):
        del_seq = seq[:cut_index] + seq[cut_index + n :]
        ins_seq = seq[:cut_index] + ("A" * n) + seq[cut_index:]
        rows.extend(
            [
                {"Edit": f"del{n}", "Protein": safe_translate(del_seq), "Frameshift": len(del_seq) % 3 != 0},
                {"Edit": f"ins{n}", "Protein": safe_translate(ins_seq), "Frameshift": len(ins_seq) % 3 != 0},
            ]
        )
    return pd.DataFrame(rows)
