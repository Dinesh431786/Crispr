import streamlit as st
import pandas as pd
from utils import validate_sequence, load_fasta, visualize_guide_location
from analysis import (
    find_gRNAs,
    find_off_targets_detailed,
    simulate_protein_edit,
    diff_proteins,
    indel_simulations,
    predict_hdr_repair,
    hybrid_score,
    ml_gRNA_score,
)
import datetime

# --- Score Explanations ---
SCORE_SUMMARY = """
#### Understanding the Scores

| Score Name       | What It Means                                                | Range      | How to Use                            |
|------------------|-------------------------------------------------------------|------------|---------------------------------------|
| **Hybrid Score** | Lab/practical score: GC%, homopolymers, seed region, off-targets | 0.0–1.0    | Higher is better. Aim for >0.85       |
| **ML Score**     | Data-driven/ML: Patterns from large CRISPR screens          | 0.0–1.0    | Higher is better. Aim for >0.7        |
| **Consensus Score** | Average of Hybrid + ML (simple and robust for most users) | 0.0–1.0    | Higher is better. Aim for >0.75       |

**Best gRNAs score high across all metrics! For simplicity, sort by Consensus Score.**
"""

# --- App Setup ---
st.set_page_config(page_title="🧬 CRISPR Lab NextGen", layout="wide")
st.title("🧬 CRISPR Lab NextGen – gRNA Designer & Impact Analyzer")
st.markdown(SCORE_SUMMARY)

# --- Sidebar Inputs ---
with st.sidebar:
    st.header("🧬 Sequence Input")
    uploaded = st.file_uploader("Upload .fasta", type=["fasta", "fa", "txt"])
    dna_seq = st.text_area("Or paste DNA sequence:", height=150, key="dna_seq")
    if uploaded:
        seq, err = load_fasta(uploaded)
        if err:
            st.error(err)
        else:
            dna_seq = seq

    pam_label = st.selectbox("PAM", ["Cas9 NGG", "Cas9 NAG", "Cas12a TTTV"], key="pam")
    GUIDE_TYPES = {"Cas9 NGG":"NGG","Cas9 NAG":"NAG","Cas12a TTTV":"TTTV"}
    pam = GUIDE_TYPES[pam_label]

    guide_len = st.slider("Guide length", 18, 25, 20, key="guide_len")
    min_gc = st.slider("Min GC %", 30, 60, 40, key="min_gc")
    max_gc = st.slider("Max GC %", 60, 80, 70, key="max_gc")

    bg_seq = st.text_area("Background DNA (off-target)", height=100, key="bg_seq")
    max_mm = st.slider("Max mismatches", 0, 4, 2, key="max_mm")

    edit_offset = st.slider(
        "Edit offset from PAM", 0, guide_len, guide_len, key="edit_offset",
        help="Cas9 cut ≈ 3 bp upstream of PAM; set as needed."
    )

    st.header("🤖 AI Settings")
    ai_backend = st.selectbox("AI Backend", ["Gemini", "OpenAI"], key="ai_backend_sidebar")
    if ai_backend == "Gemini":
        gemini_model = st.selectbox(
            "Gemini Model",
            ["gemini-1.5-flash-latest", "gemini-1.5-pro-latest"],
            key="gemini_model_sidebar",
        )
    api_key = st.text_input("API Key", type="password", key="api_key_sidebar")
    if api_key and len(api_key.strip()) > 10:
        st.success(f"{ai_backend} API initialized!", icon="✅")

# --- Session State Initialization ---
for k in (
    "df_guides", "offtargets", "guide_scores",
    "selected_gRNA", "selected_edit",
    "sim_result", "sim_indel", "ai_response"
):
    st.session_state.setdefault(k, None)

# --- gRNA Search Button ---
if st.button("🔍 Find gRNAs"):
    ok, msg = validate_sequence(dna_seq)
    if not ok:
        st.error(msg)
        st.session_state.df_guides = None
    else:
        with st.spinner("Searching gRNAs…"):
            df = find_gRNAs(dna_seq, pam, guide_len, min_gc, max_gc)
            df["HybridScore"] = [hybrid_score(g) for g in df.gRNA]
            df["MLScore"] = [ml_gRNA_score(g) for g in df.gRNA]
            df["ConsensusScore"] = (df["HybridScore"] + df["MLScore"]) / 2
            # Sort by Consensus Score by default
            df = df.sort_values("ConsensusScore", ascending=False).reset_index(drop=True)
            st.session_state.df_guides = df
        st.session_state.update(
            offtargets=None, guide_scores=None,
            sim_result=None, sim_indel=None, ai_response=None
        )

df = st.session_state.df_guides
if df is None or df.empty:
    st.info("Paste DNA & click **Find gRNAs** to begin.")
    st.stop()

st.success(f"✅ {len(df)} gRNAs found")
st.dataframe(df, use_container_width=True)
st.download_button("⬇️ Download gRNAs CSV", df.to_csv(index=False), "guides.csv")

# --- Tabs for Workflow ---
tab_ot, tab_sim, tab_ai, tab_vis, tab_rank = st.tabs([
    "Off-targets", "Simulation & Indel", "AI Explain", "Visualization", "Ranking"
])

with tab_ot:
    if st.button("Scan off-targets"):
        st.session_state.offtargets = find_off_targets_detailed(df, bg_seq, max_mm)
        st.session_state.guide_scores = {
            g: round(
                1.0 / (1 + st.session_state.offtargets[st.session_state.offtargets.gRNA == g].Mismatches.sum()),
                3
            )
            for g in df.gRNA
        }
    if st.session_state.offtargets is not None:
        if st.session_state.offtargets.empty:
            st.info("No off-targets within given mismatches.")
        else:
            st.dataframe(st.session_state.offtargets, use_container_width=True)
            st.download_button("⬇️ Download off-targets", st.session_state.offtargets.to_csv(index=False), "offtargets.csv")

with tab_sim:
    sel = st.selectbox("gRNA", df.gRNA.tolist(), key="sel_gRNA")
    edit_type = st.selectbox("Edit type", list(analysis.EDIT_TYPES.keys()), key="sel_edit")
    st.session_state.selected_gRNA = sel
    st.session_state.selected_edit = edit_type
    sub_from = sub_to = ""
    if analysis.EDIT_TYPES[edit_type] == "subAG":
        sub_from = st.text_input("Sub FROM", "A")
        sub_to = st.text_input("Sub TO", "T")

    if st.button("Simulate"):
        idx = dna_seq.upper().find(sel)
        if idx == -1:
            st.error("gRNA not found in sequence!")
        else:
            st.session_state.sim_result = simulate_protein_edit(
                dna_seq, idx + edit_offset,
                analysis.EDIT_TYPES[edit_type],
                sub_from=sub_from, sub_to=sub_to
            )
            st.session_state.sim_indel = indel_simulations(dna_seq, idx + edit_offset)

    if st.session_state.sim_result:
        before, after, fs, stop = st.session_state.sim_result
        st.markdown(f"**Before:** `{before}` → **After:** `{after}`")
        st.markdown(f"**Diff:** {diff_proteins(before, after)}")
        st.write(f"Frameshift: {fs} | Premature stop: {stop}")
    if st.session_state.sim_indel is not None:
        st.subheader("±1–3 bp Indel Simulation")
        st.dataframe(st.session_state.sim_indel, use_container_width=True)

with tab_ai:
    st.header("🔍 AI Explanation")
    st.write("AI analysis will include all three scoring metrics for comprehensive gRNA evaluation.")
    prompt = build_gemini_prompt()  # reuse helper from before
    if st.button("Ask AI"):
        if not api_key or len(api_key.strip()) < 10:
            st.error("Provide a valid API key.")
        else:
            st.session_state.ai_response = run_ai_prompt(api_key, ai_backend, prompt)
    if st.session_state.ai_response:
        st.info(st.session_state.ai_response)

with tab_vis:
    idx = dna_seq.upper().find(st.session_state.selected_gRNA or "")
    if idx >= 0:
        ax = visualize_guide_location(dna_seq, st.session_state.selected_gRNA, idx)
        st.pyplot(ax.figure)
    else:
        st.info("gRNA not found for visualization.")

with tab_rank:
    if st.session_state.guide_scores:
        rank_df = pd.DataFrame([
            {"gRNA": g, "Specificity": s}
            for g, s in st.session_state.guide_scores.items()
        ]).sort_values("Specificity", ascending=False)
        st.dataframe(rank_df, use_container_width=True)
    else:
        st.info("Run off-targets scan to rank guides.")
