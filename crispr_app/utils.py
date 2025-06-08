import re
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord
import plotly.graph_objects as go

### --- SEQUENCE VALIDATION --- ###
def validate_sequence(dna_seq, allow_n=False):
    seq = dna_seq.upper().replace('\n', '').replace(' ', '')
    if len(seq) < 23:
        return False, "Sequence too short. Minimum 23bp required."
    pattern = r'[ATCGN]+' if allow_n else r'[ATCG]+'
    if not re.fullmatch(pattern, seq):
        return False, f"Invalid characters in DNA sequence. Only A/T/C/G{'/N' if allow_n else ''} allowed."
    return True, seq

### --- FASTA LOADING --- ###
def load_fasta(uploaded_file):
    try:
        records = list(SeqIO.parse(uploaded_file, "fasta"))
        if not records:
            return None, "No FASTA records found."
        if len(records) > 1:
            return str(records[0].seq), "Warning: Multiple FASTA records found. Only first record loaded."
        return str(records[0].seq), ""
    except Exception as e:
        try:
            uploaded_file.seek(0)
            seq = uploaded_file.read().decode('utf-8').strip()
            if seq and re.fullmatch(r'[ATCGN\n\r ]+', seq.upper()):
                return seq.replace('\n', '').replace('\r', '').replace(' ', ''), ""
        except Exception:
            pass
        return None, f"Error parsing FASTA or plain sequence: {e}"

### --- DNA GUIDE VISUALIZATION --- ###
def visualize_guide_location(dna_seq, guide, start_index, pam_seq=None, strand='+'):
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

### --- ADVANCED PROTEIN DOMAIN VISUALIZATION (Plotly) --- ###
def plot_protein_domains(protein, domains, cut_site=None):
    fig = go.Figure()
    # Protein backbone
    fig.add_shape(type="rect", x0=0, x1=len(protein), y0=0.4, y1=0.6, fillcolor="lightgray", line_width=0)
    # Domains
    for _, d in domains.iterrows():
        fig.add_shape(type="rect",
                      x0=d["StartAA"], x1=d["EndAA"], y0=0.3, y1=0.7,
                      fillcolor="orange", line_color="darkorange", opacity=0.7)
        fig.add_annotation(x=(d["StartAA"]+d["EndAA"])/2, y=0.8, text=d["Domain"], showarrow=False, font_size=10)
    # Cut site
    if cut_site is not None:
        fig.add_shape(type="line", x0=cut_site, x1=cut_site, y0=0.2, y1=0.8, line_color="red", line_width=3)
    fig.update_yaxes(visible=False)
    fig.update_layout(height=160, width=650, margin=dict(l=10, r=10, t=10, b=10))
    return fig
