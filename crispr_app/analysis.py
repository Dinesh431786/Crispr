from Bio.Seq import Seq
import pandas as pd

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

def find_off_targets_detailed(guides, background_seq, max_mismatches=2):
    # Improved: returns ALL windows up to max_mismatches, never just perfect matches!
    results = []
    bg_seq = background_seq.upper().replace('\n', '').replace(' ', '')[:1_000_000]
    for _, row in guides.iterrows():
        guide = row["gRNA"].upper()
        guide_len = len(guide)
        for i in range(len(bg_seq) - guide_len + 1):
            window = bg_seq[i:i+guide_len]
            if len(window) != guide_len:
                continue
            mismatches = count_mismatches(guide, window)
            if mismatches <= max_mismatches:
                results.append({
                    "gRNA": guide,
                    "OffTargetPos": i,
                    "Mismatches": mismatches,
                    "TargetSeq": window
                })
    return pd.DataFrame(results)

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