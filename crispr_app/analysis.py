from Bio.Seq import Seq
import pandas as pd
import streamlit as st

def hybrid_score(guide, off_target_count=0):
    gc = (guide.count('G') + guide.count('C')) / len(guide)
    seed = guide[-4:]
    score = 1.0
    if gc < 0.4 or gc > 0.7:
        score -= 0.2
    if "TTTT" in guide or "GGGG" in guide:
        score -= 0.1
    if seed.count("T") + seed.count("A") > 2:
        score -= 0.1
    if guide[-1] == "G":
        score += 0.05
    score -= 0.05 * off_target_count
    score = max(score, 0.0)
    score = min(score, 1.0)
    return round(score, 3)

def ml_gRNA_score(guide):
    gc = (guide.count('G') + guide.count('C')) / len(guide)
    score = 0.5
    if 0.40 < gc < 0.60:
        score += 0.2
    elif 0.35 < gc <= 0.40 or 0.60 <= gc < 0.65:
        score += 0.1
    seed = guide[-4:]
    if seed.count('T') + seed.count('A') <= 1:
        score += 0.1
    for b in "ATCG":
        if b*4 in guide:
            score -= 0.1
    if guide[-1] == "G":
        score += 0.05
    if guide[0] == "T":
        score -= 0.05
    score = max(0.0, min(score, 1.0))
    return round(score, 3)

def check_pam(pam_seq, pam):
    pam_seq = pam_seq.upper()
    if pam == "NGG":
        return len(pam_seq) == 3 and pam_seq[1:] == "GG"
    elif pam == "NAG":
        return len(pam_seq) == 3 and pam_seq[1:] == "AG"
    elif pam == "NG":
        return len(pam_seq) == 2 and pam_seq[1] == "G"
    elif pam == "TTTV":
        return len(pam_seq) == 4 and pam_seq[:3] == "TTT" and pam_seq[3] in "ACG"
    return False

@st.cache_data
def find_gRNAs(dna_seq, pam="NGG", guide_length=20, min_gc=40, max_gc=70, add_5prime_g=False):
    sequence = dna_seq.upper().replace("\n", "").replace(" ", "")
    pam_len = 2 if pam == "NG" else (4 if pam == "TTTV" else 3)
    guides = []
    # Forward strand
    for i in range(len(sequence) - guide_length - pam_len + 1):
        guide = sequence[i:i+guide_length]
        pam_seq = sequence[i+guide_length:i+guide_length+pam_len]
        if check_pam(pam_seq, pam):
            gc = (guide.count('G') + guide.count('C')) / guide_length * 100
            if min_gc <= gc <= max_gc and "TTTT" not in guide:
                g_out = guide
                if add_5prime_g and not guide.startswith("G"):
                    g_out = "G" + guide[:-1]
                guides.append({
                    "Strand": "+",
                    "Start": i,
                    "gRNA": g_out,
                    "PAM": pam_seq,
                    "GC%": round(gc,2),
                })
    # Reverse complement
    rc_sequence = str(Seq(sequence).reverse_complement())
    for i in range(len(rc_sequence) - guide_length - pam_len + 1):
        guide = rc_sequence[i:i+guide_length]
        pam_seq = rc_sequence[i+guide_length:i+guide_length+pam_len]
        if check_pam(pam_seq, pam):
            gc = (guide.count('G') + guide.count('C')) / guide_length * 100
            if min_gc <= gc <= max_gc and "TTTT" not in guide:
                g_out = guide
                if add_5prime_g and not guide.startswith("G"):
                    g_out = "G" + guide[:-1]
                guides.append({
                    "Strand": "-",
                    "Start": len(sequence)-i-guide_length-pam_len,
                    "gRNA": g_out,
                    "PAM": pam_seq,
                    "GC%": round(gc,2),
                })
    return pd.DataFrame(guides)

def count_mismatches(a, b):
    # Strict mismatch: must be same length and both uppercase
    a = a.upper()
    b = b.upper()
    if len(a) != len(b):
        return float('inf')
    return sum(1 for x, y in zip(a, b) if x != y)

# CFD mismatch scores from CRISPOR
CFD_SCORES = {
    'rA:dA': 1.0, 'rA:dC': 1.0, 'rA:dG': 0.857, 'rA:dT': 1.0,
    'rC:dA': 1.0, 'rC:dC': 0.913, 'rC:dG': 1.0, 'rC:dT': 1.0,
    'rG:dA': 1.0, 'rG:dC': 1.0, 'rG:dG': 0.714, 'rG:dT': 0.9,
    'rU:dA': 1.0, 'rU:dC': 0.957, 'rU:dG': 0.857, 'rU:dT': 1.0
}
PAM_SCORES = {
    'AA': 0.0, 'AC': 0.0, 'AG': 0.259, 'AT': 0.0,
    'CA': 0.0, 'CC': 0.0, 'CG': 0.107, 'CT': 0.0,
    'GA': 0.069, 'GC': 0.022, 'GG': 1.0, 'GT': 0.016,
    'TA': 0.0, 'TC': 0.0, 'TG': 0.039, 'TT': 0.0
}

def calculate_cfd_score(guide_seq, off_target_seq, pam):
    """
    Calculates the Cutting Frequency Determination (CFD) score for a given gRNA and off-target sequence.
    """
    score = 1.0
    guide_seq = guide_seq.upper().replace('T', 'U')
    off_target_seq = off_target_seq.upper()

    for i in range(len(guide_seq)):
        if guide_seq[i] != off_target_seq[i]:
            key = 'r' + guide_seq[i] + ':d' + off_target_seq[i]
            score *= CFD_SCORES.get(key, 1.0)

    pam_key = pam[1:]
    score *= PAM_SCORES.get(pam_key, 0.0)
    return score

@st.cache_data
def find_off_targets_detailed(guides, background_seq, max_mismatches=2, pam="NGG"):
    """
    Finds off-target sites for a list of gRNAs in a background sequence.
    Caches results to avoid re-computation for the same inputs.
    """
    results = []
    bg_seq = background_seq.upper().replace('\n', '').replace(' ', '')[:1_000_000]
    pam_len = 3 if pam in ["NGG", "NAG"] else (2 if pam == "NG" else 4)

    # If not in cache, compute
    for guide_seq in guides["gRNA"].unique():
        guide_len = len(guide_seq)
        for i in range(len(bg_seq) - guide_len - pam_len + 1):
            window = bg_seq[i:i+guide_len]
            mismatches = count_mismatches(guide_seq, window)
            if 0 < mismatches <= max_mismatches:
                pam_seq = bg_seq[i+guide_len:i+guide_len+pam_len]
                cfd_score = calculate_cfd_score(guide_seq, window, pam_seq)
                results.append({
                    "gRNA": guide_seq,
                    "OffTargetPos": i,
                    "Mismatches": mismatches,
                    "TargetSeq": window,
                    "PAM": pam_seq,
                    "CFD_Score": cfd_score
                })

    df_results = pd.DataFrame(results)
    return df_results

def safe_translate(seq):
    extra = len(seq) % 3
    if extra != 0:
        seq += "N" * (3 - extra)
    try:
        return str(Seq(seq).translate(to_stop=True))
    except Exception:
        return "Translation failed (invalid codon)"

def simulate_protein_edit(seq, cut_index, edit_type="del1", insert_base="A", sub_from="A", sub_to="T"):
    seq = seq.upper().replace('\n', '').replace(' ', '')
    before = seq
    after = seq
    if edit_type == "del1":
        after = seq[:cut_index] + seq[cut_index+1:]
    elif edit_type == "insA":
        after = seq[:cut_index] + insert_base + seq[cut_index:]
    elif edit_type.startswith("del"):
        del_len = int(edit_type[3:])
        after = seq[:cut_index] + seq[cut_index+del_len:]
    elif edit_type.startswith("ins"):
        insert_base = edit_type[3:]
        after = seq[:cut_index] + insert_base + seq[cut_index:]
    elif edit_type == "subAG":
        if cut_index + len(sub_from) <= len(seq):
            after = seq[:cut_index] + sub_to + seq[cut_index+len(sub_from):]
    try:
        prot_before = safe_translate(before)
        prot_after = safe_translate(after)
    except Exception:
        prot_before = prot_after = "Translation failed"
    stop_lost = "*" in prot_after and "*" not in prot_before
    frameshift = (len(after) % 3) != (len(before) % 3)
    return prot_before, prot_after, frameshift, stop_lost

def diff_proteins(before, after):
    result = []
    for i, (a, b) in enumerate(zip(before, after), 1):
        if a != b:
            result.append(f"Position {i}: {a} → {b}")
    if len(before) != len(after):
        result.append(f"Length changed: {len(before)} → {len(after)}")
    return "; ".join(result) if result else "No change"

def indel_simulations(seq, cut_index):
    results = []
    for n in [1,2,3]:
        del_seq = seq[:cut_index] + seq[cut_index+n:]
        del_prot = safe_translate(del_seq)
        del_fs = (len(del_seq)%3 != 0)
        results.append({
            "Edit": f"del{n}",
            "Protein": del_prot,
            "Frameshift": del_fs
        })
        ins_seq = seq[:cut_index] + ("A"*n) + seq[cut_index:]
        ins_prot = safe_translate(ins_seq)
        ins_fs = (len(ins_seq)%3 != 0)
        results.append({
            "Edit": f"ins{n}",
            "Protein": ins_prot,
            "Frameshift": ins_fs
        })
    return pd.DataFrame(results)

def predict_hdr_repair(seq, cut_pos):
    if cut_pos >= len(seq) - 1:
        return "uncertain"
    after = seq[:cut_pos] + seq[cut_pos+1:]
    before_prot = safe_translate(seq)
    after_prot = safe_translate(after)
    if len(after) % 3 != len(seq) % 3:
        if "*" in after_prot and "*" not in before_prot:
            return "frameshift/early stop"
        return "frameshift"
    else:
        if after_prot == before_prot:
            return "silent"
        else:
            return "in-frame"