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

# Sidebar: Sequence input
with st.sidebar:
    st.header("🧬 Sequence Input")
    uploaded = st.file_uploader("Upload .fasta file", type=["fasta", "fa", "txt"])
    dna_seq = st.text_area("Or paste DNA sequence:", height=150, key="dna_seq_area")
    if uploaded:
        loaded_seq, err = load_fasta(uploaded)
        if err:
            st.error(err)
        else:
            dna_seq = loaded_seq
    pam_label = st.selectbox("Select PAM Site", options=list(GUIDE_TYPES.keys()), key="pam_site")
    pam = GUIDE_TYPES[pam_label]
    guide_length = st.slider("Guide Length", 18, 25, 20, key="guide_length")
    min_gc = st.slider("Minimum GC%", 30, 60, 40, key="min_gc")
    max_gc = st.slider("Maximum GC%", 60, 80, 70, key="max_gc")
    bg_seq = st.text_area("Background DNA (for off-target)", height=100, key="bg_seq")
    max_mismatches = st.slider("Max Mismatches (Off-target)", 0, 4, 2, key="max_mismatches")

# Session state setup
for key in [
    "df_guides", "ai_response", "offtarget_df",
    "selected_gRNA", "selected_edit", "sim_result", "sim_indel"
]:
    if key not in st.session_state:
        st.session_state[key] = None

if st.button("🔍 Find gRNAs", key="find_grnas"):
    valid, msg = validate_sequence(dna_seq)
    if not valid:
        st.error(msg)
        st.session_state.df_guides = None
    else:
        with st.spinner("Finding gRNAs..."):
            st.session_state.df_guides = find_gRNAs(
                dna_seq, pam, guide_length, min_gc, max_gc
            )
        st.session_state.offtarget_df = None
        st.session_state.ai_response = ""
        st.session_state.sim_result = None
        st.session_state.sim_indel = None

# Main logic starts here
df = st.session_state.df_guides
if df is not None and not df.empty:
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

    # --- Off-targets tab ---
    with tab1:
        if bg_seq.strip():
            if st.button("Run Off-target Search", key="run_offtarget"):
                st.session_state.offtarget_df = find_off_targets_detailed(
                    df, bg_seq, max_mismatches
                )
            ot_df = st.session_state.offtarget_df
            if ot_df is not None:
                if ot_df.empty:
                    st.info("No off-targets found in background DNA.")
                else:
                    st.dataframe(ot_df)
                    st.download_button("⬇️ Download Off-targets", data=ot_df.to_csv(index=False), file_name="offtargets.csv", mime="text/csv")
            else:
                st.info("Click 'Run Off-target Search' to analyze off-targets.")
        else:
            st.info("Paste background DNA above to enable off-target search.")

    # --- Simulation & Indel tab ---
    with tab2:
        gRNA_list = df["gRNA"].tolist()
        if gRNA_list:
            if st.session_state.selected_gRNA not in gRNA_list:
                st.session_state.selected_gRNA = gRNA_list[0]
            st.session_state.selected_gRNA = st.selectbox(
                "Choose gRNA", gRNA_list,
                index=gRNA_list.index(st.session_state.selected_gRNA),
                key="choose_gRNA"
            )
            if st.session_state.selected_edit not in EDIT_TYPES:
                st.session_state.selected_edit = list(EDIT_TYPES.keys())[0]
            st.session_state.selected_edit = st.selectbox(
                "Edit Type",
                list(EDIT_TYPES.keys()),
                index=list(EDIT_TYPES.keys()).index(st.session_state.selected_edit),
                key="choose_edit"
            )
            if EDIT_TYPES[st.session_state.selected_edit] == "subAG":
                sub_from = st.text_input("Substitute from", value="A", key="sub_from")
                sub_to = st.text_input("To", value="T", key="sub_to")
            else:
                sub_from = ""
                sub_to = ""

            if st.button("Simulate Edit", key="simulate_edit"):
                cut_index = dna_seq.upper().find(st.session_state.selected_gRNA)
                if cut_index != -1:
                    prot_before, prot_after, fs, stop = simulate_protein_edit(
                        dna_seq, cut_index + guide_length, EDIT_TYPES[st.session_state.selected_edit],
                        sub_from=sub_from, sub_to=sub_to
                    )
                    st.session_state.sim_result = (prot_before, prot_after, fs, stop)
                    st.session_state.sim_indel = indel_simulations(dna_seq, cut_index + guide_length)
                else:
                    st.warning("Selected gRNA not found in sequence.")

            if st.session_state.sim_result:
                prot_before, prot_after, fs, stop = st.session_state.sim_result
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

            if st.session_state.sim_indel is not None:
                st.subheader("Indel Simulations (del/ins 1–3bp):")
                st.dataframe(st.session_state.sim_indel)
                st.download_button("⬇️ Download Indel Results", data=st.session_state.sim_indel.to_csv(index=False), file_name="indel_simulation.csv", mime="text/csv")
        else:
            st.warning("No gRNAs available for simulation.")

    # --- AI Explain tab ---
    with tab3:
        gRNA_list = df["gRNA"].tolist()
        gRNA_choice = gRNA_list[0] if gRNA_list else ""
        ai_backend = st.selectbox("AI Backend", ["Gemini", "OpenAI"], key="ai_backend")
        api_key = st.text_input("AI API Key", type="password", key="ai_api_key")
        gene_info = st.text_area("Describe the edit context", value=f"Editing at {gRNA_choice}", key="gene_info")
        if st.button("Ask AI", key="ask_ai_button"):
            # Replace this with real API call if needed!
            st.session_state.ai_response = f"Edit context for {gRNA_choice}: {gene_info} (AI explanation would appear here.)"
        if st.session_state.ai_response:
            st.info(st.session_state.ai_response)
        else:
            st.caption("Enter key and click 'Ask AI' for real-time explanation.")

    # --- Visualization tab ---
    with tab4:
        gRNA_list = df["gRNA"].tolist()
        gRNA_choice = gRNA_list[0] if gRNA_list else ""
        cut_index = dna_seq.upper().find(gRNA_choice)
        if cut_index != -1:
            ax = visualize_guide_location(dna_seq, gRNA_choice, cut_index)
            st.pyplot(ax.figure)
        else:
            st.info("Guide position could not be visualized.")

else:
    st.info("Paste a DNA sequence and click 'Find gRNAs' to begin.")

