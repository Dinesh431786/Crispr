"""Prime-editing (pegRNA) design engine.

Designs SpCas9 pegRNAs (Spacer + Scaffold + RTT + PBS) for single-base
substitutions and scores them with determinants identified by high-throughput
prime-editing screens, in particular PRIDICT2.0 (Mathis et al., Nat. Biotechnol.
2024, doi:10.1038/s41587-024-02268-2) and the original prime-editing report
(Anzalone et al. 2019, Nature 576:149):

* Primer-binding site (PBS) melting temperature optimised toward ~37 deg C.
* PBS length ~13 nt and RTT length ~12 nt as favourable defaults.
* RTT should not begin with C at its first synthesised position (3' terminus of
  the pegRNA extension), which destabilises the edited flap.
* The edit should sit close to the nick but leave >=3 nt of 3' homology in the
  RTT downstream of the edit for efficient flap resolution.
* Moderate PBS GC content is preferred.

The returned ``Score`` is a calibrated [0, 1] heuristic for ranking candidates;
experimental validation remains essential.
"""

from __future__ import annotations

import math

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

PBS_TM_OPTIMUM = 37.0       # deg C (PRIDICT2.0)
PBS_TM_TOLERANCE = 7.0
PBS_LEN_OPTIMUM = 13
RTT_LEN_OPTIMUM = 12
MIN_3PRIME_HOMOLOGY = 3     # nt of RTT downstream of the edit


def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def _pbs_tm(seq: str) -> float:
    try:
        return float(mt.Tm_NN(seq))
    except Exception:
        s = seq.upper()
        return 2 * (s.count("A") + s.count("T")) + 4 * (s.count("G") + s.count("C"))


def _peg_score(pbs_tm: float, pbs_len: int, rtt_len: int,
               edit_pos_in_rtt: int, rtt_seq: str) -> float:
    """Calibrated pegRNA score in [0, 1] from PRIDICT2.0-informed determinants."""
    logit = 0.6

    # PBS Tm Gaussian reward around 37 deg C.
    logit += 1.2 * math.exp(-((pbs_tm - PBS_TM_OPTIMUM) / PBS_TM_TOLERANCE) ** 2)

    # Length preferences (mild quadratic penalties).
    logit -= 0.04 * (pbs_len - PBS_LEN_OPTIMUM) ** 2
    logit -= 0.03 * (rtt_len - RTT_LEN_OPTIMUM) ** 2

    # Edit should be near the nick but keep 3' homology in the RTT.
    homology_3p = rtt_len - edit_pos_in_rtt - 1
    if homology_3p < MIN_3PRIME_HOMOLOGY:
        logit -= 0.8
    logit -= 0.05 * max(0, edit_pos_in_rtt - 6)

    # Disfavour RTT beginning with C at its first synthesised position.
    if rtt_seq and rtt_seq[0] == "C":
        logit -= 0.5

    # Moderate PBS GC content (proxy via Tm already, light extra term on RTT GC).
    gc = (rtt_seq.count("G") + rtt_seq.count("C")) / max(len(rtt_seq), 1)
    logit -= 1.0 * (gc - 0.55) ** 2

    return round(_sigmoid(logit), 3)


def design_prime_editing_pegRNAs(
    dna_seq: str,
    target_pos: int,
    desired_edit: str,
    pbs_range: range = range(8, 18),
    rtt_range: range = range(10, 21),
) -> pd.DataFrame:
    """Design and rank candidate pegRNAs for a single-base substitution.

    ``target_pos`` is the 0-based index in ``dna_seq`` to edit to ``desired_edit``.
    Searches NGG PAMs within ~30 nt of the target on the (+) strand, places the
    Cas9 nick 3 bp 5' of the PAM, and enumerates PBS/RTT combinations.
    """
    dna_seq = dna_seq.upper()
    desired_edit = desired_edit.upper()
    results: list[dict] = []

    for i in range(max(0, target_pos - 30), min(len(dna_seq) - 23, target_pos + 10)):
        sub = dna_seq[i:i + 23]
        if len(sub) < 23:
            continue
        pam = sub[20:23]
        if not pam.endswith("GG"):
            continue

        spacer = sub[:20]
        nick_pos = i + 17  # 3 bp 5' of the NGG PAM
        dist_to_edit = target_pos - nick_pos
        if not (0 <= dist_to_edit <= 15):
            continue

        for pbs_len in pbs_range:
            if nick_pos - pbs_len < 0:
                continue
            pbs_seq = str(Seq(dna_seq[nick_pos - pbs_len:nick_pos]).reverse_complement())
            pbs_tm = _pbs_tm(pbs_seq)

            for rtt_len in rtt_range:
                if nick_pos + rtt_len >= len(dna_seq):
                    continue
                edit_idx_in_rtt = target_pos - nick_pos
                if not (0 <= edit_idx_in_rtt < rtt_len):
                    continue

                rtt_list = list(dna_seq[nick_pos:nick_pos + rtt_len])
                rtt_list[edit_idx_in_rtt] = desired_edit
                edited_rtt_seq = "".join(rtt_list)

                results.append({
                    "Spacer": spacer,
                    "PAM": pam,
                    "NickPos": nick_pos,
                    "PBS_Len": pbs_len,
                    "PBS_Seq": pbs_seq,
                    "PBS_Tm": round(pbs_tm, 1),
                    "RTT_Len": rtt_len,
                    "RTT_Seq": edited_rtt_seq,
                    "EditPosInRTT": edit_idx_in_rtt,
                    "Extension": edited_rtt_seq + pbs_seq,
                    "Score": _peg_score(pbs_tm, pbs_len, rtt_len,
                                        edit_idx_in_rtt, edited_rtt_seq),
                })

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values("Score", ascending=False).reset_index(drop=True)
    return df
