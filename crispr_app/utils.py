import io
import re
from Bio import SeqIO

def load_fasta(uploaded_file):
    """
    Accepts Streamlit file_uploader object (binary or text) or StringIO.
    Returns (sequence, message/warning).
    """
    try:
        # Always rewind to the start
        uploaded_file.seek(0)
        
        # Try to read the file as bytes
        file_bytes = uploaded_file.read()
        
        # Try to decode as UTF-8 (works for Streamlit uploads)
        if isinstance(file_bytes, bytes):
            text = file_bytes.decode("utf-8")
        else:
            text = file_bytes  # already str
        
        # Try to parse as FASTA
        fasta_io = io.StringIO(text)
        records = list(SeqIO.parse(fasta_io, "fasta"))
        if records and len(records) > 0:
            if len(records) > 1:
                return str(records[0].seq), "Warning: Multiple FASTA records found. Only first record loaded."
            return str(records[0].seq), ""
        
        # If not FASTA, check if it's a plain sequence
        seq = text.strip()
        if seq and re.fullmatch(r'[ATCGN\n\r ]+', seq.upper()):
            return seq.replace('\n', '').replace('\r', '').replace(' ', ''), ""
        return None, "No FASTA records found and not a plain sequence."
    
    except Exception as e:
        return None, f"Error parsing FASTA or plain sequence: {e}"
