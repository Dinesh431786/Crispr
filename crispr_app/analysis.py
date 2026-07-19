"""Core CRISPR analysis pipeline.

This module wires together the dedicated scientific engines:

* :mod:`scoring`    - calibrated on-target efficiency (Doench RS2/Azimuth-informed).
* :mod:`offtarget`  - per-site CFD + MIT/Hsu scores and aggregate specificity.
* :mod:`prime`      - PRIDICT2.0-informed pegRNA design.

It keeps the original public API (``find_gRNAs``, ``find_off_targets_detailed``,
``simulate_protein_edit``, ``indel_simulations``, ``design_prime_editing_pegRNAs``)
so the FastAPI layer is unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from Bio.Seq import Seq

try:  # Works both as top-level modules (uvicorn main:app) and as a package (tests).
    from offtarget import (
        CFD_SCORES,
        PAM_SCORES,
        aggregate_specificity,
        calculate_cfd_score,
        mit_hit_score,
    )
    from prime import design_prime_editing_pegRNAs
    from scoring import gc_fraction as _gc_fraction
    from scoring import on_target_score
    from models import predict_on_target, predict_interval
    from crisprscan import score_from_context as crisprscan_context
    from base_edit import assess as base_edit_assess
except ImportError:  # pragma: no cover - import-context fallback
    from .offtarget import (
        CFD_SCORES,
        PAM_SCORES,
        aggregate_specificity,
        calculate_cfd_score,
        mit_hit_score,
    )
    from .prime import design_prime_editing_pegRNAs
    from .scoring import gc_fraction as _gc_fraction
    from .scoring import on_target_score
    from .models import predict_on_target, predict_interval
    from .crisprscan import score_from_context as crisprscan_context
    from .base_edit import assess as base_edit_assess

__all__ = [
    "ScoringConfig",
    "hybrid_score",
    "ml_gRNA_score",
    "on_target_score",
    "check_pam",
    "find_gRNAs",
    "count_mismatches",
    "calculate_cfd_score",
    "find_off_targets_detailed",
    "summarize_specificity",
    "safe_translate",
    "simulate_protein_edit",
    "indel_simulations",
    "design_prime_editing_pegRNAs",
    "CFD_SCORES",
    "PAM_SCORES",
]


@dataclass(frozen=True)
class ScoringConfig:
    """Research-informed weights for the lightweight rule-based hybrid score."""

    gc_low: float = 0.40
    gc_high: float = 0.65
    seed_window: int = 8
    homopolymer_penalty: float = 0.15
    off_target_penalty: float = 0.06


def hybrid_score(guide: str, off_target_count: int = 0, cfg: ScoringConfig = ScoringConfig()) -> float:
    """Fast, interpretable rule-based score complementing the calibrated model."""
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


def ml_gRNA_score(guide: str, ngg_context: str | None = None) -> float:
    """Calibrated on-target efficiency (delegates to :mod:`scoring`)."""
    return on_target_score(guide, ngg_context)


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


def _pam_len(pam: str) -> int:
    return 2 if pam == "NG" else (4 if pam == "TTTV" else 3)


def _rank_by(df: pd.DataFrame, fitness) -> pd.DataFrame:
    """Sort by a goal-specific fitness vector without altering displayed scores."""
    return df.assign(_fit=np.asarray(fitness)).sort_values("_fit", ascending=False) \
             .drop(columns="_fit").reset_index(drop=True)


_STRATEGIES = ("balanced", "conservative", "robust", "optimistic")


def _strategy_value(df: pd.DataFrame, strategy: str) -> np.ndarray:
    """Ranking value = point score adjusted by the guide's conformal uncertainty.

    Turns the calibrated interval into a first-class *ranking* dimension:
    - balanced     : the ConsensusScore itself (default; uncertainty-agnostic)
    - conservative : score - half-width  (rank by the pessimistic bound)
    - robust       : score - 0.5*half-width  (penalise, don't fully subtract)
    - optimistic   : score + half-width  (rank by the best case)

    With normalized conformal the half-width varies per guide, so these genuinely
    reorder the list (they collapse to `balanced` only for a constant-width model).
    """
    base = df["ConsensusScore"].astype(float).to_numpy()
    if strategy == "balanced" or "CI_halfwidth" not in df.columns:
        return base
    hw = df["CI_halfwidth"].fillna(0.0).astype(float).to_numpy()
    if strategy == "conservative":
        return base - hw
    if strategy == "robust":
        return base - 0.5 * hw
    if strategy == "optimistic":
        return base + hw
    return base


def _base_edit_factor(guide: str, editor: str = "any") -> tuple[float, str]:
    """(base-edit fitness in [0,1], tag) for the chosen editor(s).

    Uses the full window-efficiency x editing-purity composite from
    :mod:`base_edit` (penalises bystander edits), not just window centrality.
    0 when the guide has no editable base for the chosen editor.
    """
    editors = ("CBE", "ABE") if editor == "any" else (editor,)
    best = max((base_edit_assess(guide, e) for e in editors),
               key=lambda a: a["be_score"], default=None)
    if best is None or not best["editable"]:
        return 0.0, ""
    tag = f"{best['editor'][0]}{best['target_pos']}"
    if best["bystanders"]:
        tag += f" +{len(best['bystanders'])}by"
    return float(best["be_score"]), tag


def find_gRNAs(
    dna_seq: str,
    pam: str = "NGG",
    guide_length: int = 20,
    min_gc: int = 40,
    max_gc: int = 70,
    add_5prime_g: bool = False,
    goal: str = "general",
    target_pos: int | None = None,
    editor: str = "any",
    ranking_strategy: str = "balanced",
) -> pd.DataFrame:
    sequence = dna_seq.upper().replace("\n", "").replace(" ", "")
    pam_len = _pam_len(pam)
    is_5prime_pam = pam == "TTTV"  # Cas12a: PAM is 5' of the protospacer
    guides: list[dict] = []

    use_crisprscan = pam == "NGG" and guide_length == 20 and not is_5prime_pam

    def _append_guide(strand: str, start: int, guide: str, pam_seq: str,
                      ngg_ctx: str | None, local_seq: str, local_i: int) -> None:
        gc = _gc_fraction(guide) * 100
        if min_gc <= gc <= max_gc and "TTTT" not in guide:
            g_out = guide if not add_5prime_g or guide.startswith("G") else f"G{guide[:-1]}"
            # Flanking context (6 nt 5', PAM+6 nt 3') — 3'-PAM enzymes only.
            if is_5prime_pam:
                up = down = ""
            else:
                up = local_seq[max(0, local_i - 6):local_i]
                down = local_seq[local_i + guide_length:local_i + guide_length + 9]
            on = predict_on_target(g_out, ngg_ctx, goal=goal, up=up, down=down)
            hy = hybrid_score(g_out)
            cs = crisprscan_context(local_seq, local_i, guide_length) if use_crisprscan else None
            # Consensus: blend the surrogate/learned score with the peer-reviewed
            # CRISPRscan model when its 35-mer context is available.
            if cs is not None:
                consensus = round(0.25 * hy + 0.35 * on + 0.40 * cs, 3)
            else:
                consensus = round(0.35 * hy + 0.65 * on, 3)
            row = {
                "Strand": strand,
                "Start": start,
                "gRNA": g_out,
                "PAM": pam_seq,
                "GC%": round(gc, 2),
                "HybridScore": hy,
                "MLScore": on,
                "OnTargetScore": on,
                "CRISPRScanScore": cs if cs is not None else "",
                "ConsensusScore": consensus,
            }
            # Calibrated interval from the goal-appropriate model, using the same
            # flanking context as the score (knockout -> OOF model, else cutting).
            ci = predict_interval(g_out, ngg_ctx, goal=("knockout" if goal == "knockout" else "general"),
                                  up=up, down=down)
            if ci is not None:
                row["CI_low"] = round(ci["low"] * 100)
                row["CI_high"] = round(ci["high"] * 100)
                row["CI_level"] = int(ci["coverage"] * 100)
                # Per-guide conformal half-width (0-1) -> the uncertainty signal
                # that drives uncertainty-aware ranking. Falls back to the clamped
                # displayed interval if the model predates normalized conformal.
                row["CI_halfwidth"] = round(float(ci.get("halfwidth", (ci["high"] - ci["low"]) / 2)), 3)
            guides.append(row)

    def _scan(seq: str, strand: str) -> None:
        n = len(seq)
        if is_5prime_pam:
            for i in range(n - guide_length - pam_len + 1):
                pam_seq = seq[i:i + pam_len]
                guide = seq[i + pam_len:i + pam_len + guide_length]
                if check_pam(pam_seq, pam):
                    start = i if strand == "+" else n - i - guide_length - pam_len
                    _append_guide(strand, start, guide, pam_seq, None, seq, i + pam_len)
        else:
            for i in range(n - guide_length - pam_len + 1):
                guide = seq[i:i + guide_length]
                pam_seq = seq[i + guide_length:i + guide_length + pam_len]
                if check_pam(pam_seq, pam):
                    ctx = seq[i + guide_length + pam_len:i + guide_length + pam_len + 1] or None
                    start = i if strand == "+" else n - i - guide_length - pam_len
                    _append_guide(strand, start, guide, pam_seq, ctx, seq, i)

    _scan(sequence, "+")
    _scan(str(Seq(sequence).reverse_complement()), "-")

    cols = ["Strand", "Start", "gRNA", "PAM", "GC%", "HybridScore", "MLScore", "OnTargetScore", "CRISPRScanScore", "ConsensusScore"]
    if not guides:
        return pd.DataFrame(columns=cols)
    df = pd.DataFrame(guides)

    strategy = ranking_strategy if ranking_strategy in _STRATEGIES else "balanced"
    # Uncertainty-aware ranking value (composes with the goal-specific factors
    # below; for 'general' it is the sole sort key). Displayed scores are unchanged.
    df["RankValue"] = _strategy_value(df, strategy)

    # Knock-in (HDR) mode: HDR efficiency falls off sharply with the distance
    # between the Cas9 cut and the intended edit site (~e-fold per ~10 bp), so
    # rank by cutting score weighted by proximity to target_pos. The displayed
    # ConsensusScore stays the (interpretable) cutting score; CutDist is exposed
    # separately. The cut is ~3 bp inside the protospacer from the PAM-proximal end.
    if goal == "knockin" and target_pos is not None:
        cut = np.where(df["Strand"] == "+", df["Start"] + guide_length - 3, df["Start"] + 3)
        df["CutDist"] = np.abs(cut - int(target_pos)).astype(int)
        fitness = df["RankValue"] * np.exp(-df["CutDist"] / 10.0)
        return _rank_by(df, fitness)

    # CRISPRi/a: dCas9 binds (does not cut) near the TSS; the effective window is
    # broad (~hundreds of bp), so rank by binding-site proximity with a wide decay.
    if goal == "crispri" and target_pos is not None:
        mid = df["Start"] + guide_length // 2
        df["TSSDist"] = np.abs(mid - int(target_pos)).astype(int)
        fitness = df["RankValue"] * np.exp(-df["TSSDist"] / 75.0)
        return _rank_by(df, fitness)

    # Base editing: only guides with a target base (C for CBE, A for ABE) inside
    # the activity window are usable; rank by on-target x window-centrality (peak
    # ~position 6). Guides with no editable base sink to the bottom.
    if goal == "baseedit":
        facs, tags = [], []
        for g in df["gRNA"]:
            fac, tag = _base_edit_factor(g, editor)
            facs.append(fac); tags.append(tag)
        df["BE_targets"] = tags
        fitness = df["RankValue"] * np.array(facs)
        return _rank_by(df, fitness)

    return df.sort_values("RankValue", ascending=False).reset_index(drop=True)


def _seq_identity(a: str, b: str) -> float:
    """Position-wise identity fraction (0-1) between two spacers, a proxy for
    library-recombination risk (near-duplicate guides recombine in pooled libraries)."""
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    return sum(1 for i in range(n) if a[i] == b[i]) / n


def design_multiplex(
    dna_seq: str,
    n_guides: int = 6,
    pam: str = "NGG",
    guide_length: int = 20,
    min_gc: int = 40,
    max_gc: int = 70,
    goal: str = "general",
    ranking_strategy: str = "balanced",
    diversity_weight: float = 0.5,
    min_spacing: int = 0,
    pool_size: int = 120,
) -> tuple[pd.DataFrame, dict]:
    """Greedy marginal-gain designer for a multiplex / pooled guide set.

    Picks ``n_guides`` that are individually strong yet collectively *diverse*:
    at each step it adds the candidate maximising

        marginal_gain = on_target_value - diversity_weight * max_identity_to_selected

    where on_target_value is the (strategy-adjusted) ranking value and
    max_identity_to_selected is the highest position-wise sequence identity to an
    already-chosen guide (a recombination-risk proxy). Optional ``min_spacing``
    forbids two guides whose cut sites sit within N bp, avoiding redundant cuts.
    Fully transparent: every pick reports its marginal gain and why it was chosen.
    """
    df = find_gRNAs(dna_seq, pam=pam, guide_length=guide_length, min_gc=min_gc,
                    max_gc=max_gc, goal=goal, ranking_strategy=ranking_strategy)
    if df.empty:
        return df, {"requested": n_guides, "selected": 0, "note": "no candidate guides"}

    pool = df.head(pool_size).reset_index(drop=True)
    value = pool["RankValue"].astype(float).to_numpy() if "RankValue" in pool.columns \
        else pool["ConsensusScore"].astype(float).to_numpy()
    starts = pool["Start"].astype(int).to_numpy()
    seqs = pool["gRNA"].tolist()

    chosen: list[int] = []
    remaining = set(range(len(pool)))
    records = []
    while len(chosen) < n_guides and remaining:
        best_i, best_gain, best_sim = None, -1e9, 0.0
        for i in list(remaining):
            if min_spacing > 0 and any(abs(int(starts[i]) - int(starts[j])) < min_spacing for j in chosen):
                continue
            sim = max((_seq_identity(seqs[i], seqs[j]) for j in chosen), default=0.0)
            gain = value[i] - diversity_weight * sim
            if gain > best_gain:
                best_i, best_gain, best_sim = i, gain, sim
        if best_i is None:  # everything left violates spacing
            break
        chosen.append(best_i)
        remaining.discard(best_i)
        g = pool.iloc[best_i]
        reason = (f"score {g['ConsensusScore']:.2f}" +
                  (f", {best_sim:.0%} max similarity to the set" if chosen[:-1] else ", seed (highest score)"))
        records.append({**g.to_dict(),
                        "MarginalGain": round(float(best_gain), 3),
                        "MaxSimilarity": round(float(best_sim), 3),
                        "Reason": reason})

    picked = pd.DataFrame(records)
    # Collective summary.
    if len(picked) >= 2:
        idents = [_seq_identity(records[a]["gRNA"], records[b]["gRNA"])
                  for a in range(len(records)) for b in range(a + 1, len(records))]
        max_ident = round(max(idents), 3)
        mean_div = round(1 - float(np.mean(idents)), 3)
    else:
        max_ident, mean_div = 0.0, 1.0
    sel_starts = [r["Start"] for r in records]
    summary = {
        "requested": n_guides,
        "selected": len(picked),
        "pool_size": int(len(pool)),
        "mean_on_target": round(float(pool.loc[[r_i for r_i in chosen], "ConsensusScore"].mean()), 3) if chosen else 0.0,
        "mean_diversity": mean_div,          # 1 - mean pairwise identity (higher = more diverse)
        "max_pairwise_identity": max_ident,  # recombination-risk proxy (lower = safer)
        "span_bp": int(max(sel_starts) - min(sel_starts)) if sel_starts else 0,
        "diversity_weight": diversity_weight,
        "min_spacing": min_spacing,
        "ranking_strategy": ranking_strategy,
    }
    return picked, summary


def count_mismatches(a: str, b: str) -> int:
    a = a.upper()
    b = b.upper()
    if len(a) != len(b):
        return 10 ** 9
    return sum(1 for x, y in zip(a, b) if x != y)


def _scan_strand(bg_seq: str, bg_arr: np.ndarray, guide_seq: str, guide_arr: np.ndarray,
                 guide_len: int, pam_len: int, pam: str, max_mismatches: int,
                 strand: str, bg_len_total: int) -> list[dict]:
    rows: list[dict] = []
    usable = len(bg_arr) - pam_len
    if usable <= guide_len:
        return rows
    windows = np.lib.stride_tricks.sliding_window_view(bg_arr[:usable], guide_len)
    mismatch_counts = np.sum(windows != guide_arr, axis=1)
    match_indices = np.where((mismatch_counts > 0) & (mismatch_counts <= max_mismatches))[0]

    for idx in match_indices:
        idx = int(idx)
        window_seq = bg_seq[idx:idx + guide_len]
        pam_seq = bg_seq[idx + guide_len:idx + guide_len + pam_len]
        if check_pam(pam_seq, pam):
            # Report position in original (+) coordinates.
            pos = idx if strand == "+" else bg_len_total - idx - guide_len - pam_len
            rows.append({
                "gRNA": guide_seq,
                "Strand": strand,
                "OffTargetPos": pos,
                "Mismatches": int(mismatch_counts[idx]),
                "TargetSeq": window_seq,
                "PAM": pam_seq,
                "CFD_Score": calculate_cfd_score(guide_seq, window_seq, pam_seq),
                "MIT_Score": mit_hit_score(guide_seq, window_seq),
            })
    return rows


def find_off_targets_detailed(guides: pd.DataFrame, background_seq: str, max_mismatches: int = 2, pam: str = "NGG") -> pd.DataFrame:
    """Vectorised, both-strand off-target search with CFD + MIT scoring.

    Uses NumPy sliding-window views for sub-second scanning of kilobase-scale
    backgrounds and scans BOTH the forward and reverse-complement strands so
    that off-targets on either strand are detected (a correctness improvement
    over forward-only scanning).
    """
    bg_seq = background_seq.upper().replace("\n", "").replace(" ", "")[:1_000_000]
    bg_len_total = len(bg_seq)
    bg_arr = np.frombuffer(bg_seq.encode(), dtype=np.int8)

    rc_seq = str(Seq(bg_seq).reverse_complement())
    rc_arr = np.frombuffer(rc_seq.encode(), dtype=np.int8)

    pam_len = 3 if pam in ("NGG", "NAG") else (2 if pam == "NG" else 4)

    results: list[dict] = []
    for guide_seq in guides["gRNA"].unique():
        guide_len = len(guide_seq)
        guide_arr = np.frombuffer(guide_seq.encode(), dtype=np.int8)
        results.extend(_scan_strand(bg_seq, bg_arr, guide_seq, guide_arr, guide_len,
                                    pam_len, pam, max_mismatches, "+", bg_len_total))
        results.extend(_scan_strand(rc_seq, rc_arr, guide_seq, guide_arr, guide_len,
                                    pam_len, pam, max_mismatches, "-", bg_len_total))

    return pd.DataFrame(results)


def summarize_specificity(off_targets: pd.DataFrame, guides: list[str]) -> pd.DataFrame:
    """Aggregate per-guide specificity (CRISPOR-style) from off-target hits."""
    rows = []
    for g in guides:
        sub = off_targets[off_targets["gRNA"] == g] if not off_targets.empty else off_targets
        cfd_list = list(sub["CFD_Score"]) if not sub.empty else []
        mit_list = list(sub["MIT_Score"]) if not sub.empty else []
        rows.append({
            "gRNA": g,
            "OffTargetCount": int(len(sub)),
            "CFD_Specificity": aggregate_specificity(cfd_list, scale=100.0),
            "MIT_Specificity": aggregate_specificity(mit_list, scale=1.0),
        })
    return pd.DataFrame(rows)


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
    seq = seq.upper().replace("\n", "").replace(" ", "")
    before = seq
    after = seq
    if edit_type == "del1":
        after = seq[:cut_index] + seq[cut_index + 1:]
    elif edit_type == "insA":
        after = seq[:cut_index] + insert_base + seq[cut_index:]
    elif edit_type.startswith("del"):
        del_len = int(edit_type[3:])
        after = seq[:cut_index] + seq[cut_index + del_len:]
    elif edit_type.startswith("ins"):
        ins = edit_type[3:] or insert_base
        after = seq[:cut_index] + ins + seq[cut_index:]
    elif edit_type == "subAG" and cut_index + len(sub_from) <= len(seq):
        after = seq[:cut_index] + sub_to + seq[cut_index + len(sub_from):]

    prot_before = safe_translate(before)
    prot_after = safe_translate(after)
    frameshift = (len(after) % 3) != (len(before) % 3)
    stop_lost = "*" in prot_after and "*" not in prot_before
    return prot_before, prot_after, frameshift, stop_lost


def indel_simulations(seq: str, cut_index: int) -> pd.DataFrame:
    rows = []
    for n in (1, 2, 3):
        del_seq = seq[:cut_index] + seq[cut_index + n:]
        ins_seq = seq[:cut_index] + ("A" * n) + seq[cut_index:]
        rows.extend([
            {"Edit": f"del{n}", "Protein": safe_translate(del_seq), "Frameshift": len(del_seq) % 3 != 0},
            {"Edit": f"ins{n}", "Protein": safe_translate(ins_seq), "Frameshift": len(ins_seq) % 3 != 0},
        ])
    return pd.DataFrame(rows)
