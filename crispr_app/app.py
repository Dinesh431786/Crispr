import streamlit as st
import pandas as pd
from utils import validate_sequence, load_fasta, visualize_guide_location
from analysis import find_gRNAs, find_off_targets_detailed, simulate_protein_edit, diff_proteins, indel_simulations

GUIDE_TYPES = {
    "Cas9 NGG": "NGG",
    "Cas9 NAG": "NAG",
    "Cas12a TTTV": "TTTV"
}
EDIT_TYPES = {
    "Delete 1 bp": "del1",
    "Insert A": "insA",
    "Delete 3 bp": "del3",
    "Insert G": "insG",
    "Substitute A→T": "subAG"
}

st.set_page_config(page_title="🧬 CRISPR Lab", layout="wide")
st.title("🧬 CRISPR Lab: gRNA Designer & Impact Analyzer")

with st.sidebar:
    st.header("🧬 Sequence Input")
    uploaded = st.file_uploader("Upload .fasta file", type=["fasta", "fa", "txt"])
    dna_seq = st.text_area("Or paste DNA sequence:", height=150)
    if uploaded:
        loaded_seq, err = load_fasta(uploaded)
        if err:
            st.error(err)
        else:
            dna_seq = loaded_seq
    pam_label = st.selectbox("Select PAM Site", options=list(GUIDE_TYPES.keys()))
    pam = GUIDE_TYPES[pam_label]
    guide_length = st.slider("Guide Length", 18, 25, 20)
    min_gc = st.slider("Minimum GC%", 30, 60, 40)
    max_gc = st.slider("Maximum GC%", 60, 80, 70)
    bg_seq = st.text_area("Background DNA (for off-target)", height=100)
    max_mismatches = st.slider("Max Mismatches (Off-target)", 0, 4, 2)

if "ai_response" not in st.session_state:
    st.session_state.ai_response = ""

if st.button("🔍 Find gRNAs"):
    valid, msg = validate_sequence(dna_seq)
    if not valid:
        st.error(msg)
    else:
        with st.spinner("Finding gRNAs..."):
            df = find_gRNAs(dna_seq, pam, guide_length, min_gc, max_gc)
        if df.empty:
            st.warning("No valid gRNAs found. Try changing PAM or sequence.")
        else:
            st.success(f"✅ Found {len(df)} gRNAs")
            st.dataframe(df)
            csv = df.to_csv(index=False)
            st.download_button("⬇️ Download gRNAs", data=csv, file_name="guides.csv", mime="text/csv")

            tab1, tab2, tab3, tab4 = st.tabs([
                "🔍 Off-targets (detailed)",
                "🧬 Simulation & Indel",
                "🤖 AI Explain",
                "🖼️ Visualization"
            ])

            with tab1:
                if bg_seq.strip():
                    ot_df = find_off_targets_detailed(df, bg_seq, max_mismatches)
                    if ot_df.empty:
                        st.info("No off-targets found in background DNA.")
                    else:
                        st.dataframe(ot_df)
                        st.download_button("⬇️ Download Off-targets", data=ot_df.to_csv(index=False), file_name="offtargets.csv", mime="text/csv")
                else:
                    st.info("Paste background DNA above to enable off-target search.")

            with tab2:
                gRNA_list = df["gRNA"].tolist()
                gRNA_choice = st.selectbox("Choose gRNA", gRNA_list)
                cut_index = dna_seq.upper().find(gRNA_choice)
                if cut_index != -1:
                    edit_label = st.selectbox("Edit Type", options=list(EDIT_TYPES.keys()))
                    edit_type = EDIT_TYPES[edit_label]
                    if edit_type == "subAG":
                        sub_from = st.text_input("Substitute from", value="A")
                        sub_to = st.text_input("To", value="T")
                    else:
                        sub_from = sub_to = ""
                    prot_before, prot_after, fs, stop = simulate_protein_edit(
                        dna_seq, cut_index + guide_length, edit_type, sub_from=sub_from, sub_to=sub_to
                    )
                    st.markdown(f"**Before:** `{prot_before}`")
                    st.markdown(f"**After:** `{prot_after}`")
                    st.markdown(f"**Diff:** {diff_proteins(prot_before, prot_after)}")
                    st.markdown(f"**Frameshift:** {'Yes' if fs else 'No'}")
                    st.markdown(f"**Premature Stop:** {'Yes' if stop else 'No'}")
                    report = f"""
Protein Before Edit: {prot_before}
Protein After Edit: {prot_after}
Frameshift: {'Yes' if fs else 'No'}
Premature Stop Codon: {'Yes' if stop else 'No'}
"""
                    st.download_button("⬇️ Download Protein Report", report, file_name="protein_report.txt")

                    # Show indel simulation table
                    st.subheader("Indel Simulations (del/ins 1–3bp):")
                    indel_df = indel_simulations(dna_seq, cut_index + guide_length)
                    st.dataframe(indel_df)
                    st.download_button("⬇️ Download Indel Results", data=indel_df.to_csv(index=False), file_name="indel_simulation.csv", mime="text/csv")
                else:
                    st.warning("Selected gRNA not found in sequence.")

            with tab3:
                gRNA_list = df["gRNA"].tolist()
                gRNA_choice = gRNA_list[0] if gRNA_list else ""
                ai_backend = st.selectbox("AI Backend", ["Gemini", "OpenAI"], key="ai_backend")
                api_key = st.text_input("AI API Key", type="password", key="ai_api_key")
                gene_info = st.text_area("Describe the edit context", value=f"Editing at {gRNA_choice}", key="gene_info")
                
                if st.button("Ask AI", key="ask_ai_button"):
                    # Replace this with a real API call if needed
                    st.session_state.ai_response = f"Edit context for {gRNA_choice}: {gene_info} (AI explanation would appear here.)"
                
                # Always show the response if it exists
                if st.session_state.ai_response:
                    st.info(st.session_state.ai_response)
                else:
                    st.caption("Enter key and click 'Ask AI' for real-time explanation.")

            with tab4:
                gRNA_list = df["gRNA"].tolist()
                gRNA_choice = gRNA_list[0] if gRNA_list else ""
                cut_index = dna_seq.upper().find(gRNA_choice)
                if cut_index != -1:
                    ax = visualize_guide_location(dna_seq, gRNA_choice, cut_index)
                    st.pyplot(ax.figure)
                else:
                    st.info("Guide position could not be visualized.")
