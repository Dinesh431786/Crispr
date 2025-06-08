from Bio.Seq import Seq
from difflib import Differ
import pandas as pd

### --- ADVANCED gRNA SCORING (Consensus/ML + rules) --- ###
def consensus_guide_score(guide):
    # Real: consensus from multiple scoring methods, can be expanded
    score = score_guide_advanced(guide)
    # Optionally: integrate external ML models, CRSeek, Azimuth, etc.
    return score

def score_guide_advanced(guide):
    gc = (guide.count('G') + guide.count('C')) / len(guide)
    seed = guide[-4:]  # PAM-proximal, most important
    bad_seed = seed.count('T') + seed.count('A')
    score = 1.0
    if gc < 0.4 or gc > 0.7:
        score -= 0.3
    if bad_seed >= 2:
        score -= 0.2
    if "TTTT" in guide or "GGGG" in guide or "AAAA" in guide or "CCCC" in guide:
        score -= 0.2
    if guide[-1] == "G":
        score += 0.1
    if guide[0] == "T":
        score -= 0.1
    # Can expand: secondary structure, repeat, self-complementarity
    return round(max(score, 0.0), 2)

def check_pam(pam_seq, pam):
    if pam == "NGG":
        return pam_seq[1:] == "GG"
    elif pam == "NAG":
        return pam_seq[1:] == "AG"
    elif pam == "TTTV":
        return pam_seq[:3] == "TTT" and pam_seq[3] in "ACG"
    return False

def find_gRNAs(dna_seq, pam="NGG", guide_length=20, min_gc=40, max_gc=70):
    sequence = dna_seq.upper().replace("\n", "").replace(" ", "")
    pam_len = len(pam)
    guides = []
    for i in range(len(sequence) - guide_length - pam_len + 1):
        guide = sequence[i:i+guide_length]
        pam_seq = sequence[i+guide_length:i+guide_length+pam_len]
        if check_pam(pam_seq, pam):
            gc = (guide.count('G') + guide.count('C')) / guide_length * 100
            if min_gc <= gc <= max_gc and "TTTT" not in guide:
                guides.append({
                    "Strand": "+",
                    "Start": i,
                    "gRNA": guide,
                    "PAM": pam_seq,
                    "GC%": round(gc,2),
                    "ActivityScore": score_guide_advanced(guide)
                })
    rc_sequence = str(Seq(sequence).reverse_complement())
    for i in range(len(rc_sequence) - guide_length - pam_len + 1):
        guide = rc_sequence[i:i+guide_length]
        pam_seq = rc_sequence[i+guide_length:i+guide_length+pam_len]
        if check_pam(pam_seq, pam):
            gc = (guide.count('G') + guide.count('C')) / guide_length * 100
            if min_gc <= gc <= max_gc and "TTTT" not in guide:
                guides.append({
                    "Strand": "-",
                    "Start": len(sequence)-i-guide_length-pam_len,
                    "gRNA": guide,
                    "PAM": pam_seq,
                    "GC%": round(gc,2),
                    "ActivityScore": score_guide_advanced(guide)
                })
    return pd.DataFrame(guides)

### --- OFF-TARGETS --- ###
def find_off_targets_detailed(guides, background_seq, max_mismatches=2):
    results = []
    bg_seq = background_seq.upper().replace('\n', '')[:1_000_000]
    for _, row in guides.iterrows():
        guide = row["gRNA"]
        details = []
        for i in range(len(bg_seq) - len(guide) + 1):
            window = bg_seq[i:i+len(guide)]
            mismatches = sum(1 for a, b in zip(guide, window) if a != b)
            if mismatches <= max_mismatches:
                details.append({
                    "Position": i,
                    "Mismatches": mismatches,
                    "Sequence": window
                })
        results.append({
            "gRNA": guide,
            "OffTargetCount": len(details),
            "Details": details
        })
    # Flatten for Streamlit
    flat = []
    for result in results:
        for detail in result["Details"]:
            flat.append({
                "gRNA": result["gRNA"],
                "OffTargetPos": detail["Position"],
                "Mismatches": detail["Mismatches"],
                "TargetSeq": detail["Sequence"]
            })
    return pd.DataFrame(flat)

### --- PROTEIN EDIT SIMULATION --- ###
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
    elif edit_type == "insa":
        after = seq[:cut_index] + insert_base + seq[cut_index:]
    elif edit_type.startswith("del"):
        del_len = int(edit_type[3:])
        after = seq[:cut_index] + seq[cut_index+del_len:]
    elif edit_type.startswith("ins"):
        insert_base = edit_type[3:]
        after = seq[:cut_index] + insert_base + seq[cut_index:]
    elif edit_type == "subag":
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

def indel_simulations(seq, cut_index):
    # Simulate −3, −2, −1, +1, +2, +3 edits
    results = []
    for n in [1,2,3]:
        # Deletions
        del_seq = seq[:cut_index] + seq[cut_index+n:]
        del_prot = safe_translate(del_seq)
        del_fs = (len(del_seq)%3 != 0)
        results.append({
            "Edit": f"del{n}",
            "Protein": del_prot,
            "Frameshift": del_fs
        })
        # Insertions (A)
        ins_seq = seq[:cut_index] + ("A"*n) + seq[cut_index:]
        ins_prot = safe_translate(ins_seq)
        ins_fs = (len(ins_seq)%3 != 0)
        results.append({
            "Edit": f"ins{n}",
            "Protein": ins_prot,
            "Frameshift": ins_fs
        })
    return pd.DataFrame(results)

### --- PROTEIN DOMAIN ANNOTATION (can use Uniprot/PFAM APIs in future) --- ###
def annotate_protein_domains(seq):
    # Real: Placeholder for Uniprot/PFAM API annotation
    # Here, simply returns the length and a mock domain. You can extend with true annotation via API.
    protein = safe_translate(seq)
    if len(protein) > 0 and protein != "Translation failed (invalid codon)":
        return pd.DataFrame([{
            "Domain": "ExampleDomain",
            "StartAA": 10,
            "EndAA": 50,
            "Function": "Put real API call or annotation here",
        }])
    return None
