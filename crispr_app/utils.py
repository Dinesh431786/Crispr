import re
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord

def validate_sequence(dna_seq):
    seq = dna_seq.upper().replace('\n', '').replace(' ', '')
    if len(seq) < 23:
        return False, "Sequence too short. Minimum 23bp required."
    if not re.fullmatch(r'[ATCG]+', seq):
        return False, "Invalid characters in DNA sequence. Only A/T/C/G allowed."
    return True, ""

def load_fasta(uploaded_file):
    try:
        records = list(SeqIO.parse(uploaded_file, "fasta"))
        if not records:
            return None, "No FASTA records found."
        return str(records[0].seq), ""
    except Exception as e:
        return None, f"Error parsing FASTA: {e}"

def visualize_guide_location(dna_seq, guide, start_index):
    features = [GraphicFeature(start=start_index, end=start_index + len(guide), label=guide, color="#ffd700")]
    record = GraphicRecord(sequence=dna_seq, features=features)
    ax, _ = record.plot(figure_width=10)
    return ax
