import re
import io
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord
import plotly.graph_objects as go

def validate_sequence(dna_seq, allow_n=False):
    """
    Validates a DNA sequence string.
    Returns (True, cleaned_sequence) or (False, error_message).
    """
    if dna_seq is None:
        return False, "No DNA sequence provided."
    seq = dna_seq.upper().replace('\n', '').replace(' ', '')
    if len(seq) < 23:
        return False, "Sequence too short. Minimum 23bp required."
    pattern = r'[ATCGN]+' if allow_n else r'[ATCG]+'
    if not re.fullmatch(pattern, seq):
        return False, f"Invalid characters in DNA sequence. Only A/T/C/G{'/N' if allow_n else ''} allowed."
    return True, seq

def load_fasta(uploaded_file):
    """
    Loads a FASTA file or plain DNA sequence from an uploaded file object.
    Supports Streamlit's file_uploader (binary or text).
    Returns (sequence, message/warning).
    """
    try:
        if uploaded_file is None:
            return "", "No file uploaded."
        uploaded_file.seek(0)
        file_bytes = uploaded_file.read()

        # Try decoding if bytes
        if isinstance(file_bytes, bytes):
            text = file_bytes.decode("utf-8", errors="replace")
        else:
            text = file_bytes  # already str

        # Try FASTA parsing
        fasta_io = io.StringIO(text)
        records = list(SeqIO.parse(fasta_io, "fasta"))
        if records and len(records) > 0:
            if len(records) > 1:
                return str(records[0].seq), "Warning: Multiple FASTA records found. Only first record loaded."
            return str(records[0].seq), ""
        
        # If not FASTA, check for plain sequence
        seq = text.strip()
        if seq and re.fullmatch(r'[ATCGN\n\r ]+', seq.upper()):
            return seq.replace('\n', '').replace('\r', '').replace(' ', ''), ""
        return "", "No FASTA records found and not a plain sequence."

    except Exception as e:
        return "", f"Error parsing FASTA or plain sequence: {e}"

def visualize_guide_location(dna_seq, guide, start_index, pam_seq=None, strand='+'):
    """
    Draws a gRNA location (plus optional PAM) on the DNA using dna_features_viewer.
    """
    try:
        features = [
            GraphicFeature(
                start=start_index,
                end=start_index + len(guide),
                label=f"gRNA {guide}",
                color="#ffd700",
                strand=+1 if strand == '+' else -1
            )
        ]
        if pam_seq:
            features.append(
                GraphicFeature(
                    start=start_index + len(guide),
                    end=start_index + len(guide) + len(pam_seq),
                    label=f"PAM {pam_seq}",
                    color="#9acd32",
                    strand=+1 if strand == '+' else -1
                )
            )
        record = GraphicRecord(sequence=dna_seq, features=features)
        ax, _ = record.plot(figure_width=10)
        return ax
    except Exception as e:
        print(f"Visualization error: {e}")
        return None

def plot_protein_domains(protein, domains, cut_site=None):
    """
    Visualizes protein domains (as rectangles) on a 1D plot, with optional cut-site indicator.
    """
    fig = go.Figure()
    fig.add_shape(type="rect", x0=0, x1=len(protein), y0=0.4, y1=0.6, fillcolor="lightgray", line_width=0)
    if domains is not None and not domains.empty:
        for _, d in domains.iterrows():
            fig.add_shape(type="rect",
                          x0=d["StartAA"], x1=d["EndAA"], y0=0.3, y1=0.7,
                          fillcolor="orange", line_color="darkorange", opacity=0.7)
            fig.add_annotation(x=(d["StartAA"]+d["EndAA"])/2, y=0.8, text=d["Domain"], showarrow=False, font_size=10)
    if cut_site is not None:
        fig.add_shape(type="line", x0=cut_site, x1=cut_site, y0=0.2, y1=0.8, line_color="red", line_width=3)
    fig.update_yaxes(visible=False)
    fig.update_layout(height=160, width=650, margin=dict(l=10, r=10, t=10, b=10))
    return fig