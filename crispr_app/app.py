import streamlit as st
import pandas as pd
from utils import validate_sequence, load_fasta
from analysis import (
    find_gRNAs, find_off_targets_detailed,
    simulate_protein_edit, diff_proteins, indel_simulations,
    hybrid_score, ml_gRNA_score,
)

# --- Streamlit page setup ---
st.set_page_config(page_title="🧬 CRISPR gRNA Designer", layout="wide")

# --- Wizard navigation ---
steps = [
    "1️⃣ DNA Input",
    "2️⃣ gRNA Design",
    "3️⃣ Edit Simulation",
    "4️⃣ Off-targets (optional)",
    "5️⃣ AI Explain (optional)"
]
step = st.sidebar.radio("Workflow Steps", steps)

# --- Session state for cross-step data ---
for k in (
    "dna_seq", "df_guides", "df_display", "offtargets",
    "guide_scores", "selected_gRNA", "sim_result", "sim_indel"
):
    st.session_state.setdefault(k, None)

# --- Step 1: DNA Input ---
if step == steps[0]:
    st.title("🧬 Step 1: Enter Your DNA Sequence")
    st.markdown("Upload or paste your DNA. You can use any species or construct!")

    uploaded = st.file_uploader("Upload FASTA or text file", type=["fasta", "fa", "txt"])
    dna_seq = st.text_area("Or paste DNA sequence below:", height=120, key="input_dna_seq")
    if uploaded:
        seq, err = load_fasta(uploaded)
        if err:
            st.error(err)
        else:
            dna_seq = seq
            st.success("Loaded sequence from file!")

    st.session_state.dna_seq = dna_seq.strip()
    if st.button("Proceed to gRNA Design ➡️", type="primary"):
        if not dna_seq or len(dna_seq) < 23:
            st.error("Enter at least 23 bp of DNA to proceed.")
        else:
            st.session_state.dna_seq = dna_seq.strip()
            st.success("DNA accepted! Move to Step 2 (gRNA Design) in sidebar.")

# --- Step 2: gRNA Design ---
elif step == steps[1]:
    st.title("🧬 Step 2: Find gRNAs")
    st.markdown("Tune your design parameters, then click **Find gRNAs**.")

    # --- Sidebar params ---
    st.sidebar.header("gRNA Design Settings")
    pam_label = st.sidebar.selectbox("PAM", ["Cas9 NGG", "Cas9 NAG", "Cas9 NG", "Cas12a TTTV"], key="pam")
    GUIDE_TYPES = {"Cas9 NGG": "NGG", "Cas9 NAG": "NAG", "Cas9 NG": "NG", "Cas12a TTTV": "TTTV"}
    pam = GUIDE_TYPES[pam_label]

    u6_toggle = st.sidebar.toggle("U6 Promoter (add G at 5’ if needed)", value=False, key="u6_toggle")
    guide_len = st.sidebar.slider("Guide length", 18, 25, 20, key="guide_len")
    min_gc = st.sidebar.slider("Min GC %", 30, 60, 40, key="min_gc")
    max_gc = st.sidebar.slider("Max GC %", 60, 80, 70, key="max_gc")

    st.info("**Note:** Hybrid Score = empirical lab rules (GC, seed, homopolymers, PAM, etc). ML Score is 'ML-inspired', not a trained model.")

    # --- Main action ---
    if st.button("🔍 Find gRNAs"):
        ok, msg = validate_sequence(st.session_state.dna_seq)
        if not ok:
            st.error(msg)
            st.session_state.df_guides = None
        else:
            with st.spinner("Searching for gRNAs..."):
                df = find_gRNAs(
                    st.session_state.dna_seq, pam, guide_len, min_gc, max_gc
                )
                if df.empty:
                    st.error("No gRNAs found with these settings.")
                else:
                    # Scoring (if missing)
                    if "HybridScore" not in df.columns:
                        df["HybridScore"] = [hybrid_score(g) for g in df.gRNA]
                    if "MLScore" not in df.columns:
                        df["MLScore"] = [ml_gRNA_score(g) for g in df.gRNA]
                    df["ConsensusScore"] = ((df["HybridScore"] + df["MLScore"]) / 2).clip(upper=1.0)
                    # U6 toggle
                    if u6_toggle:
                        df["gRNA"] = df["gRNA"].apply(lambda g: g if g.startswith("G") else "G" + g[:-1])
                    st.session_state.df_guides = df
                    st.success(f"{len(df)} gRNAs found!")
    df = st.session_state.df_guides
    if df is not None and not df.empty:
        show_cols = ["Strand", "Start", "gRNA", "PAM", "GC%", "HybridScore", "MLScore", "ConsensusScore"]
        st.dataframe(df[show_cols], use_container_width=True, height=400)
        st.download_button("⬇️ Download gRNAs CSV", df[show_cols].to_csv(index=False), "guides.csv")

# --- Step 3: Edit Simulation ---
elif step == steps[2]:
    st.title("🧬 Step 3: Simulate Edit/Indel Effects")
    df = st.session_state.df_guides
    if df is None or df.empty:
        st.warning("Design gRNAs first! Go back to Step 2.")
        st.stop()
    u6_toggle = st.session_state.get("u6_toggle", False)
    # Select gRNA
    g_list = df.gRNA.tolist()
    st.session_state.selected_gRNA = st.selectbox("Choose a gRNA for simulation", g_list, key="sel_gRNA_sim")
    # Edit type
    edit_options = {
        "Delete 1 bp": "del1", "Insert A": "insA",
        "Delete 3 bp": "del3", "Insert G": "insG", "Substitute A→T": "subAG"
    }
    st.session_state.selected_edit = st.selectbox("Edit type", list(edit_options), key="edit_type")
    sub_from = sub_to = ""
    if edit_options[st.session_state.selected_edit] == "subAG":
        sub_from = st.text_input("Sub FROM", "A")
        sub_to = st.text_input("Sub TO", "T")
    # Simulate button
    if st.button("Simulate Edit"):
        idx = st.session_state.dna_seq.upper().find(st.session_state.selected_gRNA)
        if idx == -1:
            st.error("gRNA not found in sequence!")
        else:
            st.session_state.sim_result = simulate_protein_edit(
                st.session_state.dna_seq, idx + st.session_state.get("guide_len", 20),  # edit_offset
                edit_options[st.session_state.selected_edit], sub_from=sub_from, sub_to=sub_to
            )
            st.session_state.sim_indel = indel_simulations(
                st.session_state.dna_seq, idx + st.session_state.get("guide_len", 20)
            )
    # Show results
    if st.session_state.sim_result:
        before, after, fs, stop = st.session_state.sim_result
        st.markdown(f"**Before protein:** <span style='color:green'>{before}</span>", unsafe_allow_html=True)
        st.markdown(f"**After protein:** <span style='color:blue'>{after}</span>", unsafe_allow_html=True)
        st.markdown(f"**Diff:** {diff_proteins(before, after)}")
        st.write("Frameshift:", fs, "| Premature stop:", stop)
    if st.session_state.sim_indel is not None:
        st.subheader("±1–3 bp indel simulation")
        st.dataframe(st.session_state.sim_indel, use_container_width=True, height=250)

# --- Step 4: Off-targets (optional) ---
elif step == steps[3]:
    st.title("🧬 Step 4: Off-target Scan (Optional)")
    df = st.session_state.df_guides
    if df is None or df.empty:
        st.warning("Design gRNAs first! Go back to Step 2.")
        st.stop()
    st.sidebar.header("Off-target Scan Settings")
    bg_seq = st.sidebar.text_area("Background DNA (off-target)", key="bg_seq")
    max_mm = st.sidebar.slider("Max mismatches", 0, 4, 2, key="max_mm")
    if st.button("Scan off-targets"):
        result_from_find = find_off_targets_detailed(df, bg_seq, max_mm)
        st.session_state.offtargets = result_from_find
    ot_df = st.session_state.offtargets
    if ot_df is not None:
        st.dataframe(ot_df, use_container_width=True, height=350)
        st.download_button("⬇️ Download off-targets", ot_df.to_csv(index=False), "offtargets.csv")

# --- Step 5: AI Explain (Optional) ---
elif step == steps[4]:
    st.title("🧬 Step 5: AI Explain (Optional)")
    st.info("Let Gemini/OpenAI summarize or explain your results! Paste your API key to activate.")
    # Advanced users only: AI backend/settings panel
    st.sidebar.header("AI Explain Settings")
    ai_backend = st.sidebar.selectbox("AI Backend", ["Gemini", "OpenAI"], key="ai_backend")
    gemini_model = st.sidebar.selectbox(
        "Gemini Model", [
            "gemini-1.5-flash-latest", "gemini-1.5-pro-latest",
            "gemini-pro", "gemini-1.0-pro-latest"
        ], key="gemini_model"
    ) if ai_backend == "Gemini" else None
    api_key = st.sidebar.text_input("API Key", type="password", key="api_key")
    st.info("AI summary/report feature coming soon! 🚧")
    # You can connect the AI feature as you had it before here.

# ---- END ----

st.markdown("""
<style>
/* Cleaner cards & spacing */
.stButton>button {font-size: 1.1rem; padding: 0.6em 1.2em;}
.stDataFrame, .stTable {background: #fffaf4;}
</style>
""", unsafe_allow_html=True)
