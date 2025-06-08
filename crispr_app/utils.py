import re
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord

def validate_sequence(dna_seq, allow_n=False):
    seq = dna_seq.upper().replace('\n', '').replace(' ', '')
    if len(seq) < 23:
        return False, "Sequence too short. Minimum 23bp required."
    pattern = r'[ATCGN]+' if allow_n else r'[ATCG]+'
    if not re.fullmatch(pattern, seq):
        return False, "Invalid characters in DNA sequence. Only A/T/C/G{} allowed.".format("/N" if allow_n else "")
    return True, seq

def load_fasta(uploaded_file):
    try:
        records = list(SeqIO.parse(uploaded_file, "fasta"))
        if not records:
            return None, "No FASTA records found."
        if len(records) > 1:
            return str(records[0].seq), "Warning: Multiple FASTA records found. Only first record loaded."
        return str(records[0].seq), ""
    except Exception as e:
        # Try reading as plain text if not FASTA
        try:
            uploaded_file.seek(0)
            seq = uploaded_file.read().decode('utf-8').strip()
            if seq and re.fullmatch(r'[ATCGN\n\r ]+', seq.upper()):
                return seq.replace('\n', '').replace('\r', '').replace(' ', ''), ""
        except Exception:
            pass
        return None, f"Error parsing FASTA or plain sequence: {e}"

def visualize_guide_location(dna_seq, guide, start_index, pam_seq=None, strand='+'):
    features = [
        GraphicFeature(start=start_index, end=start_index + len(guide), label=f"gRNA {guide}", color="#ffd700", strand=+1 if strand == '+' else -1)
    ]
    if pam_seq:
        features.append(GraphicFeature(
            start=start_index + len(guide), end=start_index + len(guide) + len(pam_seq),
            label=f"PAM {pam_seq}", color="#9acd32", strand=+1 if strand == '+' else -1
        ))
    record = GraphicRecord(sequence=dna_seq, features=features)
    ax, _ = record.plot(figure_width=10)
    return ax
