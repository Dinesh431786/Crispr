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

st.set_page_config(page_title="🧬 CRISPR Lab NextGen", layout="wide")
st.title("🧬 CRISPR Lab NextGen: gRNA Designer & Impact Analyzer")

# Sidebar: Sequence input & settings
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
    edit_offset = st.slider("Edit Offset from PAM (for sim)", 0, guide_length, guide_length, help="Where the cut happens relative to gRNA start (e.g., Cas9 is 3bp upstream of PAM).", key="edit_offset")
    st.header("🔑 AI Settings")
    st.session_state.ai_backend = st.selectbox("AI Backend", ["Gemini", "OpenAI"], key="ai_backend_sidebar")
    st.session_state.api_key = st.text_input("Gemini/OpenAI API Key", type="password", key="api_key_sidebar")

# Session state setup for ALL
for key in [
    "df_guides", "ai_response", "offtarget_df",
    "selected_gRNA", "selected_edit", "sim_result", "sim_indel", "guide_scores"
]:
    if key not in st.session_state:
        st.session_state[key] = None

# Find gRNAs button and results
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
        st.session_state.guide_scores = None

# Main logic starts here
df = st.session_state.df_guides
if df is not None and not df.empty:
    st.success(f"✅ Found {len(df)} gRNAs")
    st.dataframe(df)
    csv = df.to_csv(index=False)
    st.download_button("⬇️ Download gRNAs", data=csv, file_name="guides.csv", mime="text/csv")

    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "🔍 Off-targets (detailed)",
        "🧬 Simulation & Indel",
        "🤖 AI Explain",
        "🖼️ Visualization",
        "⭐ gRNA Ranking"
    ])

    # --- Off-targets tab ---
    with tab1:
        if bg_seq.strip():
            if st.button("Run Off-target Search", key="run_offtarget"):
                st.session_state.offtarget_df = find_off_targets_detailed(
                    df, bg_seq, max_mismatches
                )
                # --- Advanced: gRNA specificity scoring
                score_dict = {}
                ot_df = st.session_state.offtarget_df
                for guide in df["gRNA"]:
                    hits = ot_df[ot_df["gRNA"] == guide]
                    score = 1.0 if hits.empty else 1.0 / (1 + hits["Mismatches"].sum())
                    score_dict[guide] = round(score, 3)
                st.session_state.guide_scores = score_dict
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
            selected_gRNA = st.selectbox(
                "Choose gRNA",
                gRNA_list,
                index=gRNA_list.index(st.session_state.selected_gRNA),
                key="choose_gRNA"
            )
            st.session_state.selected_gRNA = selected_gRNA

            edit_options = list(EDIT_TYPES.keys())
            if st.session_state.selected_edit not in edit_options:
                st.session_state.selected_edit = edit_options[0]
            selected_edit = st.selectbox(
                "Edit Type",
                edit_options,
                index=edit_options.index(st.session_state.selected_edit),
                key="choose_edit"
            )
            st.session_state.selected_edit = selected_edit

            # Extra for substitution
            if EDIT_TYPES[selected_edit] == "subAG":
                sub_from = st.text_input("Substitute from", value="A", key="sub_from")
                sub_to = st.text_input("To", value="T", key="sub_to")
            else:
                sub_from = ""
                sub_to = ""

            if st.button("Simulate Edit", key="simulate_edit"):
                cut_index = dna_seq.upper().find(st.session_state.selected_gRNA)
                if cut_index != -1:
                    prot_before, prot_after, fs, stop = simulate_protein_edit(
                        dna_seq,
                        cut_index + edit_offset,
                        EDIT_TYPES[selected_edit],
                        sub_from=sub_from, sub_to=sub_to
                    )
                    st.session_state.sim_result = (prot_before, prot_after, fs, stop)
                    st.session_state.sim_indel = indel_simulations(dna_seq, cut_index + edit_offset)
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
        gene_info = st.text_area("Describe the edit context", value=f"Editing at {gRNA_choice}", key="gene_info")
        if st.button("Ask AI", key="ask_ai_button"):
            # Example stub: integrate your Gemini/OpenAI API here using st.session_state.api_key
            st.session_state.ai_response = (
                f"Edit context for {gRNA_choice}: {gene_info} "
                f"(AI explanation would appear here. Backend: {st.session_state.ai_backend})"
            )
        if st.session_state.ai_response:
            st.info(st.session_state.ai_response)
        else:
            st.caption("Enter description and click 'Ask AI'.")

    # --- Visualization tab ---
    with tab4:
        gRNA_list = df["gRNA"].tolist()
        gRNA_choice = st.session_state.selected_gRNA if st.session_state.selected_gRNA in gRNA_list else (gRNA_list[0] if gRNA_list else "")
        cut_index = dna_seq.upper().find(gRNA_choice)
        if cut_index != -1:
            ax = visualize_guide_location(dna_seq, gRNA_choice, cut_index)
            st.pyplot(ax.figure)
        else:
            st.info("Guide position could not be visualized.")

    # --- gRNA Ranking tab (advanced) ---
    with tab5:
        if st.session_state.guide_scores:
            st.subheader("gRNA Specificity Ranking (higher=better, fewer off-targets)")
            spec_df = pd.DataFrame([
                {"gRNA": g, "SpecificityScore": s} for g, s in st.session_state.guide_scores.items()
            ]).sort_values("SpecificityScore", ascending=False)
            st.dataframe(spec_df)
        else:
            st.info("Run Off-target Search to see gRNA specificity scores.")

else:
    st.info("Paste a DNA sequence and click 'Find gRNAs' to begin.")

