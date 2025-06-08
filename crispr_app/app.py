import streamlit as st
import pandas as pd
from utils import (
    validate_sequence,
    load_fasta,
    visualize_guide_location,
)
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

SCORE_SUMMARY = """
#### Understanding the Scores

| Score Name      | What It Means                                                | Range      | How to Use                            |
|-----------------|-------------------------------------------------------------|------------|---------------------------------------|
| Hybrid Score    | Lab-rule score: GC%, homopolymers, seed region, off-targets | 0.0–1.0    | Higher is better. Aim for >0.85       |
| ML Score        | Data-driven (AI/ML): Patterns from large CRISPR screens     | 0.0–1.0    | Higher is better. Aim for >0.7        |
| Consensus Score | The average of Hybrid & ML for balanced ranking             | 0.0–1.0    | Highest values are most reliable      |

**Best gRNAs score high in all three!**
"""

SCORE_EXPLAIN = """
**Hybrid Score:**  
Practical laboratory rule-based score, range: **0.0 (poor) to 1.0 (ideal)**  
Considers GC%, homopolymer runs, seed region, off-target count, and terminal base.  
**Interpretation:** Higher = more reliable guide for most experiments.

**ML Score:**  
Machine-learning inspired meta-score (**0.0–1.0**).  
Based on published patterns in large gRNA screens (GC%, homopolymers, seed, position, etc).  
**Interpretation:** Higher = more likely to work according to CRISPR data trends.

**Consensus Score:**  
Combines the Hybrid Score and ML Score into a single metric by averaging them:Consensus Score = (Hybrid Score + ML Score)

This provides a balanced view, giving equal importance to laboratory rules and machine-learning predictions.  
**Interpretation:** Use guides with the highest Consensus Scores for best results!
"""

st.set_page_config(page_title="🧬 CRISPR Lab NextGen", layout="wide")
st.title("🧬 CRISPR Lab NextGen – gRNA Designer & Impact Analyzer")

st.markdown(SCORE_SUMMARY)
st.markdown(SCORE_EXPLAIN)

# ---- Sidebar ----
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
    GUIDE_TYPES = {
        "Cas9 NGG": "NGG",
        "Cas9 NAG": "NAG",
        "Cas12a TTTV": "TTTV",
    }
    pam = GUIDE_TYPES[pam_label]
    guide_len = st.slider("Guide length", 18, 25, 20, key="guide_len")
    min_gc = st.slider("Min GC %", 30, 60, 40, key="min_gc")
    max_gc = st.slider("Max GC %", 60, 80, 70, key="max_gc")
    bg_seq = st.text_area("Background DNA (off-target)", height=100, key="bg_seq")
    max_mm = st.slider("Max mismatches", 0, 4, 2, key="max_mm")
    edit_offset = st.slider(
        "Edit offset from PAM",
        0,
        guide_len,
        guide_len,
        key="edit_offset",
        help="Cas9 cut ≈ 3 bp upstream of PAM; set as needed.",
    )

    st.header("🤖 AI Explain Settings")
    ai_backend = st.selectbox("AI Backend", ["Gemini", "OpenAI"], key="ai_backend_sidebar")
    gemini_model = "gemini-1.5-flash-latest"
    if ai_backend == "Gemini":
        gemini_model = st.selectbox(
            "Gemini Model",
            ["gemini-1.5-flash-latest", "gemini-1.5-pro-latest"],
            key="gemini_model_sidebar",
        )
    api_key = st.text_input("API Key", type="password", key="api_key_sidebar")
    if api_key and len(api_key.strip()) > 10:
        st.success(f"{ai_backend} API initialized!", icon="✅")

# ---- Session state holders ----
for k in (
    "df_guides",
    "offtargets",
    "guide_scores",
    "selected_gRNA",
    "selected_edit",
    "sim_result",
    "sim_indel",
    "ai_response",
    "gemini_report",
):
    st.session_state.setdefault(k, None)

# ---- gRNA search ----
if st.button("🔍 Find gRNAs"):
    ok, msg = validate_sequence(dna_seq)
    if not ok:
        st.error(msg)
        st.session_state.df_guides = None
    else:
        with st.spinner("Searching gRNAs…"):
            st.session_state.df_guides = find_gRNAs(
                dna_seq, pam, guide_len, min_gc, max_gc
            )
        st.session_state.update(
            offtargets=None,
            guide_scores=None,
            sim_result=None,
            sim_indel=None,
            ai_response="",
            gemini_report=None,
        )

df = st.session_state.df_guides
if df is None or df.empty:
    st.info("Paste DNA & click **Find gRNAs** to begin.")
    st.stop()

# ---- Add scores columns ----
if "HybridScore" not in df.columns or "MLScore" not in df.columns or "ConsensusScore" not in df.columns:
    df["HybridScore"] = [hybrid_score(g) for g in df.gRNA]
    df["MLScore"] = [ml_gRNA_score(g) for g in df.gRNA]
    df["ConsensusScore"] = (df["HybridScore"] + df["MLScore"]) / 2

st.success(f"✅ {len(df)} gRNAs found")
st.dataframe(df, use_container_width=True)
st.download_button("⬇️ Download gRNAs CSV", df.to_csv(index=False), "guides.csv")

# ---- ONE CLICK GEMINI REPORT ----
st.markdown("---")
st.header("📄 One Click Gemini Report")
st.markdown(SCORE_SUMMARY)

def build_gemini_prompt():
    context_parts = [
        SCORE_SUMMARY,
        "### Hybrid and ML Score Logic\n",
        SCORE_EXPLAIN,
        "\n\n### gRNA Candidates Table (top 10 shown)\n",
        df[["gRNA", "HybridScore", "MLScore", "ConsensusScore"]].head(10).to_csv(sep="|", index=False),
    ]
    ot_df = st.session_state.offtargets
    if ot_df is not None and not ot_df.empty:
        off_target_summary = ot_df.groupby("gRNA")["Mismatches"].count().reset_index()
        context_parts.append("\n\n### Off-target Summary\n")
        context_parts.append(off_target_summary.to_csv(sep="|", index=False))
    sim_res = st.session_state.sim_result
    if sim_res:
        before, after, fs, stop = sim_res
        context_parts.append("\n\n### Simulation Result\n")
        context_parts.append(f"Before protein: {before}\n")
        context_parts.append(f"After protein: {after}\n")
        context_parts.append(f"Frameshift: {fs} | Premature stop: {stop}\n")
    context_str = "\n".join(context_parts)
    prompt = (
        context_str
        + "\n\nSummarize the above results for a CRISPR scientist, highlighting: "
          "1. Which guides have the highest reliability and why. "
          "2. Any off-target risks. "
          "3. Editing simulation impact. "
          "4. Additional tips for experiment design."
    )
    return prompt

if st.button("📄 Generate Gemini Report"):
    ai_backend = st.session_state.get("ai_backend_sidebar", "Gemini")
    api_key = st.session_state.get("api_key_sidebar", "")
    gemini_model = st.session_state.get("gemini_model_sidebar", "gemini-1.5-flash-latest")
    if not api_key or len(api_key.strip()) < 10:
        st.error("Enter a valid API key in the sidebar.")
    else:
        prompt = build_gemini_prompt()
        try:
            if ai_backend == "Gemini":
                import google.generativeai as genai
                genai.configure(api_key=api_key)
                model = genai.GenerativeModel(gemini_model)
                result = model.generate_content(prompt)
                st.session_state.gemini_report = result.text if hasattr(result, "text") else str(result)
            else:
                import openai
                openai.api_key = api_key
                resp = openai.ChatCompletion.create(
                    model="gpt-3.5-turbo",
                    messages=[
                        {"role": "system", "content": "You are a CRISPR genome editing expert."},
                        {"role": "user", "content": prompt},
                    ],
                )
                st.session_state.gemini_report = resp.choices[0].message.content
        except Exception as e:
            import traceback
            st.session_state.gemini_report = "API error:\n" + traceback.format_exc(limit=2)

if st.session_state.gemini_report:
    st.markdown(SCORE_SUMMARY)
    st.subheader("Gemini AI Report")
    st.info(st.session_state.gemini_report)

# ---- TABS (Off-targets, Simulation, AI, Visualization, Ranking) ----
tab_ot, tab_sim, tab_ai, tab_vis, tab_rank = st.tabs(
    ["Off-targets", "Simulation & Indel", "AI Explain", "Visualization", "Ranking"]
)

with tab_ot:
    if not bg_seq.strip():
        st.info("Provide background DNA in sidebar for off-target scanning.")
    else:
        if st.button("Scan off-targets"):
            st.session_state.offtargets = find_off_targets_detailed(
                df, bg_seq, max_mm
            )
            scores = {
                g: round(
                    1.0
                    if st.session_state.offtargets[
                        st.session_state.offtargets.gRNA == g
                    ].empty
                    else 1.0
                    / (
                        1
                        + st.session_state.offtargets[
                            st.session_state.offtargets.gRNA == g
                        ]
                        .Mismatches.sum()
                    ),
                    3,
                )
                for g in df.gRNA
            }
            st.session_state.guide_scores = scores
        if st.session_state.offtargets is not None:
            if st.session_state.offtargets.empty:
                st.info("No off-targets within given mismatches.")
            else:
                st.dataframe(st.session_state.offtargets, use_container_width=True)
                st.download_button(
                    "⬇️ Download off-targets",
                    st.session_state.offtargets.to_csv(index=False),
                    "offtargets.csv",
                )

with tab_sim:
    g_list = df.gRNA.tolist()
    st.session_state.selected_gRNA = st.selectbox(
        "gRNA", g_list, key="sel_gRNA"
    )
    EDIT_TYPES = {
        "Delete 1 bp": "del1",
        "Insert A": "insA",
        "Delete 3 bp": "del3",
        "Insert G": "insG",
        "Substitute A→T": "subAG",
    }
    st.session_state.selected_edit = st.selectbox(
        "Edit type", list(EDIT_TYPES), key="sel_edit"
    )
    sub_from = sub_to = ""
    if EDIT_TYPES[st.session_state.selected_edit] == "subAG":
        sub_from = st.text_input("Sub FROM", "A")
        sub_to = st.text_input("Sub TO", "T")

    if st.button("Simulate"):
        idx = dna_seq.upper().find(st.session_state.selected_gRNA)
        if idx == -1:
            st.error("gRNA not found in sequence!")
        else:
            st.session_state.sim_result = simulate_protein_edit(
                dna_seq,
                idx + edit_offset,
                EDIT_TYPES[st.session_state.selected_edit],
                sub_from=sub_from,
                sub_to=sub_to,
            )
            st.session_state.sim_indel = indel_simulations(
                dna_seq, idx + edit_offset
            )

    if st.session_state.sim_result:
        before, after, fs, stop = st.session_state.sim_result
        st.markdown(f"**Before protein:** `{before}`")
        st.markdown(f"**After protein:** `{after}`")
        st.markdown(f"**Diff:** {diff_proteins(before, after)}")
        st.write("Frameshift:", fs, "| Premature stop:", stop)
    if st.session_state.sim_indel is not None:
        st.subheader("±1–3 bp indel simulation")
        st.dataframe(st.session_state.sim_indel, use_container_width=True)

with tab_ai:
    st.header("AI Explain (Gemini / OpenAI)")
    context_parts = [
        SCORE_SUMMARY,
        "### Hybrid and ML Score Logic\n",
        SCORE_EXPLAIN,
        "\n\n### gRNA Candidates Table (top 10 shown)\n",
        df[["gRNA", "HybridScore", "MLScore", "ConsensusScore"]].head(10).to_csv(sep="|", index=False),
    ]
    ot_df = st.session_state.offtargets
    if ot_df is not None and not ot_df.empty:
        off_target_summary = ot_df.groupby("gRNA")["Mismatches"].count().reset_index()
        context_parts.append("\n\n### Off-target Summary\n")
        context_parts.append(off_target_summary.to_csv(sep="|", index=False))
    sim_res = st.session_state.sim_result
    if sim_res:
        before, after, fs, stop = sim_res
        context_parts.append("\n\n### Simulation Result\n")
        context_parts.append(f"Before protein: {before}\n")
        context_parts.append(f"After protein: {after}\n")
        context_parts.append(f"Frameshift: {fs} | Premature stop: {stop}\n")
    context_str = "\n".join(context_parts)
    with st.expander("🔎 See full context sent to AI (for debugging)", expanded=False):
        st.code(context_str)
    user_notes = st.text_area(
        "Add any specific questions or notes for AI (optional):", 
        "", 
        key="ai_notes"
    )
    prompt = (
        context_str
        + "\n\n"
        + (user_notes.strip() if user_notes else "")
        + "\n\nSummarize the above results for a CRISPR scientist, highlighting: "
          "1. Which guides have the highest reliability and why. "
          "2. Any off-target risks. "
          "3. Editing simulation impact. "
          "4. Additional tips for experiment design."
    )
    if st.button("Ask AI"):
        ai_backend = st.session_state.get("ai_backend_sidebar", "Gemini")
        api_key = st.session_state.get("api_key_sidebar", "")
        gemini_model = st.session_state.get("gemini_model_sidebar", "gemini-1.5-flash-latest")
        if not api_key or len(api_key.strip()) < 10:
            st.error("Enter a valid API key in the sidebar.")
        else:
            try:
                if ai_backend == "Gemini":
                    import google.generativeai as genai
                    genai.configure(api_key=api_key)
                    model = genai.GenerativeModel(gemini_model)
                    result = model.generate_content(prompt)
                    st.session_state.ai_response = result.text if hasattr(result, "text") else str(result)
                else:
                    import openai
                    openai.api_key = api_key
                    resp = openai.ChatCompletion.create(
                        model="gpt-3.5-turbo",
                        messages=[
                            {"role": "system", "content": "You are a CRISPR genome editing expert."},
                            {"role": "user", "content": prompt},
                        ],
                    )
                    st.session_state.ai_response = resp.choices[0].message.content
            except Exception as e:
                import traceback
                st.session_state.ai_response = "API error:\n" + traceback.format_exc(limit=2)
    if st.session_state.ai_response:
        st.markdown(SCORE_SUMMARY)
        st.info(st.session_state.ai_response)

with tab_vis:
    idx = dna_seq.upper().find(st.session_state.selected_gRNA)
    if idx != -1:
        ax = visualize_guide_location(dna_seq, st.session_state.selected_gRNA, idx)
        st.pyplot(ax.figure)
    else:
        st.info("gRNA not found for visualization.")

with tab_rank:
    if st.session_state.guide_scores:
        rank_df = (
            pd.DataFrame(
                [
                    {"gRNA": g, "Specificity": s}
                    for g, s in st.session_state.guide_scores.items()
                ]
            )
            .sort_values("Specificity", ascending=False)
            .reset_index(drop=True)
        )
        st.dataframe(rank_df, use_container_width=True)
    else:
        st.info("Run off-target scan to get specificity ranking.")

