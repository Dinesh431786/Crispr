from Bio.Seq import Seq
from difflib import Differ
import pandas as pd

def find_gRNAs(dna_seq, pam="NGG", guide_length=20, min_gc=40, max_gc=70):
    sequence = dna_seq.upper().replace("\n", "").replace(" ", "")
    pam_len = len(pam)
    guides = []
    # Forward
    for i in range(len(sequence) - guide_length - pam_len + 1):
        guide = sequence[i:i+guide_length]
        pam_seq = sequence[i+guide_length:i+guide_length+pam_len]
        if check_pam(pam_seq, pam):
            gc = (guide.count('G') + guide.count('C')) / guide_length * 100
            if min_gc <= gc <= max_gc and "TTTT" not in guide:
                guides.append({"Strand": "+", "Start": i, "gRNA": guide, "PAM": pam_seq, "GC%": round(gc,2)})
    # Reverse
    rc_sequence = str(Seq(sequence).reverse_complement())
    for i in range(len(rc_sequence) - guide_length - pam_len + 1):
        guide = rc_sequence[i:i+guide_length]
        pam_seq = rc_sequence[i+guide_length:i+guide_length+pam_len]
        if check_pam(pam_seq, pam):
            gc = (guide.count('G') + guide.count('C')) / guide_length * 100
            if min_gc <= gc <= max_gc and "TTTT" not in guide:
                guides.append({"Strand": "-", "Start": len(sequence)-i-guide_length-pam_len, "gRNA": guide, "PAM": pam_seq, "GC%": round(gc,2)})
    return pd.DataFrame(guides)

def check_pam(pam_seq, pam):
    # Supports NGG, NAG, TTTV (V = A/C/G)
    if pam == "NGG":
        return pam_seq[1:] == "GG"
    elif pam == "NAG":
        return pam_seq[1:] == "AG"
    elif pam == "TTTV":
        return pam_seq[:3] == "TTT" and pam_seq[3] in "ACG"
    return False

def find_off_targets(guides, background_seq, max_mismatches=2):
    results = []
    bg_seq = background_seq.upper().replace('\n', '')[:1_000_000]
    for _, row in guides.iterrows():
        guide = row["gRNA"]
        matches = []
        for i in range(len(bg_seq) - len(guide) + 1):
            window = bg_seq[i:i+len(guide)]
            mismatches = sum(1 for a, b in zip(guide, window) if a != b)
            if mismatches <= max_mismatches:
                matches.append((i, mismatches))
        results.append({
            "gRNA": guide,
            "OffTargetHits": len(matches),
            "BestMismatch": min([m[1] for m in matches]) if matches else None
        })
    return pd.DataFrame(results)

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
        prot_before = str(Seq(before).translate(to_stop=True))
        prot_after = str(Seq(after).translate(to_stop=True))
    except Exception:
        prot_before = prot_after = "Translation failed"
    stop_lost = "*" in prot_after and "*" not in prot_before
    frameshift = (len(after) % 3) != (len(before) % 3)
    return prot_before, prot_after, frameshift, stop_lost

def diff_proteins(before, after):
    d = Differ()
    result = []
    for s in d.compare(before, after):
        if s.startswith('+'):
            result.append(f'**+{s[2:]}**')
        elif s.startswith('-'):
            result.append(f'~~{s[2:]}~~')
        elif s.startswith('?'):
            continue
        else:
            result.append(s)
    return ''.join(result)
