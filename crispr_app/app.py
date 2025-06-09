import streamlit as st
import pandas as pd
from utils import (
    validate_sequence,
    load_fasta,
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

| Score Name      | What It Means                                                                                                   | Range      | How to Use                                              |
|-----------------|----------------------------------------------------------------------------------------------------------------|------------|---------------------------------------------------------|
| Hybrid Score    | Lab-rule score: GC%, homopolymers, seed region, off-targets. **Rule-based, not ML.**                           | 0.0–1.0    | **Excellent:** >0.85, **Recommended:** >0.8             |
| ML Score        | Machine learning-inspired (not trained ML): Based on patterns from large CRISPR screens, using rules for GC%, seed region, homopolymers, etc. | 0.0–1.0    | **Excellent:** >0.7, **Recommended:** >0.65             |
| Consensus Score | Average of Hybrid & ML for balanced ranking                                                                    | 0.0–1.0    | **Excellent:** >0.85, **Recommended:** >0.8             |

**Aim for Consensus Score >0.85 (“Excellent”). Guides >0.8 are also “Recommended”. Lower scores may work but are not ideal.**

**Note:**  
Hybrid Score is calculated from accepted laboratory rules (GC content, homopolymers, seed region, off-targets, terminal base).  
It is **not** from machine learning or experimental data.  
ML Score is “machine learning inspired,” based on features found in published ML models, but **not** the output of a trained ML algorithm.
"""

SCORE_EXPLAIN = """
**Hybrid Score:**  
Calculated using laboratory rule-based factors (GC%, homopolymers, seed region, off-target count, terminal base).  
Range: 0.0 (poor) to 1.0 (excellent). Higher = more reliable guide.
This is a purely rule-based score, not derived from ML or experimental training.

**ML Score:**  
Based on large published CRISPR screen data using AI/ML patterns (GC%, homopolymers, seed, position, etc).  
Range: 0.0 (poor) to 1.0 (excellent). Higher = more likely to work in practice.
This score is “ML-inspired,” not from a trained ML model; it is a deterministic rule-based score.

**Consensus Score:**  
Consensus Score = (Hybrid Score + ML Score) / 2  
Averages both lab rules and ML predictions for a balanced rank.  
Higher = best chance of experimental success.
"""

st.set_page_config(page_title="🧬 CRISPR Guide RNA Designer – For Plants, Humans, Microbes & More", layout="wide")
st.title("🧬 CRISPR Guide RNA Designer")
st.markdown("#### <span style='color:#22a35d;'>For Plants, Humans, Microbes – For All DNA</span>", unsafe_allow_html=True)

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

    # U6 toggle up for better UX
    u6_g_toggle = st.toggle(
        "U6 Promoter (add G at 5’ if needed)", value=False, key="u6_g_toggle",
        help="If ON, adds a leading 'G' to each gRNA if not already present (for U6/T7 promoters)."
    )

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
            [
                "gemini-1.5-flash-latest",
                "gemini-1.5-pro-latest",
                "gemini-pro",
                "gemini-1.0-pro-latest",
            ],
            key="gemini_model_sidebar",
        )
    api_key = st.text_input("API Key", type="password", key="api_key_sidebar")
    if api_key and len(api_key.strip()) > 10:
        st.success(f"{ai_backend} API initialized!", icon="✅")

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

# Helper: U6 promoter G-adding wrapper
def u6_g_mod(seq):
    if seq.startswith("G"):
        return seq
    return "G" + seq

def apply_u6_toggle_to_df(df, u6_toggle):
    df_mod = df.copy()
    if u6_toggle:
        # Replace the gRNA column in-place (NO duplicate column bug)
        df_mod["gRNA"] = df_mod["gRNA"].apply(u6_g_mod)
    return df_mod

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

if "HybridScore" not in df.columns or "MLScore" not in df.columns or "ConsensusScore" not in df.columns:
    df["HybridScore"] = [hybrid_score(g) for g in df.gRNA]
    df["MLScore"] = [ml_gRNA_score(g) for g in df.gRNA]
    df["ConsensusScore"] = ((df["HybridScore"] + df["MLScore"]) / 2).clip(upper=1.0)

u6_toggle = st.session_state.get("u6_g_toggle", False)
df_display = apply_u6_toggle_to_df(df, u6_toggle)

st.success(f"✅ {len(df_display)} gRNAs found")
st.dataframe(df_display, use_container_width=True)
st.download_button("⬇️ Download gRNAs CSV", df_display.to_csv(index=False), "guides.csv")

st.markdown("---")
st.header("📄 One Click Gemini Report")

def build_gemini_prompt():
    context_parts = [
        SCORE_SUMMARY,
        "### Score Logic Explanation (for AI only)\n",
        SCORE_EXPLAIN,
        "\n\n### gRNA Candidates Table (top 10 shown)\n",
        df_display[["gRNA", "HybridScore", "MLScore", "ConsensusScore"]].head(10).to_csv(sep="|", index=False),
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
            error_str = str(e)
            # FRIENDLY Gemini error reporting
            if "API key not valid" in error_str or "API_KEY_INVALID" in error_str:
                st.error("❌ Your Gemini API key is invalid or this model is not enabled for your account/project. Please double-check your key and model selection.")
            elif "model not found" in error_str or "not supported" in error_str:
                st.error("❌ The selected Gemini model is not available. Try selecting 'gemini-pro' in the sidebar.")
            else:
                st.error(f"Gemini API error: {error_str}")
            st.session_state.gemini_report = ""

if st.session_state.gemini_report:
    st.subheader("Gemini AI Report")
    st.info(st.session_state.gemini_report)

tab_ot, tab_sim, tab_ai, tab_rank = st.tabs(
    ["Off-targets", "Simulation & Indel", "AI Explain", "Ranking"]
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
    # The gRNA display for selection (should show with U6 G if toggle is on)
    if u6_toggle:
        g_list_display = [u6_g_mod(g) for g in df.gRNA.tolist()]
    else:
        g_list_display = df.gRNA.tolist()
    gRNA_display_to_seq = {u6_g_mod(g) if u6_toggle else g: g for g in df.gRNA.tolist()}

    st.session_state.selected_gRNA = st.selectbox(
        "gRNA", g_list_display, key="sel_gRNA"
    )
    # Map back to original sequence for all simulations/analysis!
    gRNA_for_analysis = gRNA_display_to_seq[st.session_state.selected_gRNA]

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
        idx = dna_seq.upper().find(gRNA_for_analysis)
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
        SCORE_EXPLAIN,
        "\n\n### gRNA Candidates Table (top 10 shown)\n",
        df_display[["gRNA", "HybridScore", "MLScore", "ConsensusScore"]].head(10).to_csv(sep="|", index=False),
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
            error_str = str(e)
            if "API key not valid" in error_str or "API_KEY_INVALID" in error_str:
                st.error("❌ Your Gemini API key is invalid or this model is not enabled for your account/project. Please double-check your key and model selection.")
            elif "model not found" in error_str or "not supported" in error_str:
                st.error("❌ The selected Gemini model is not available. Try selecting 'gemini-pro' in the sidebar.")
            else:
                st.error(f"Gemini API error: {error_str}")
            st.session_state.ai_response = ""
    if st.session_state.ai_response:
        st.info(st.session_state.ai_response)

with tab_rank:
    if st.session_state.guide_scores:
        rank_df = (
            pd.DataFrame(
                [
                    {"gRNA": u6_g_mod(g) if u6_toggle else g, "Specificity": s}
                    for g, s in st.session_state.guide_scores.items()
                ]
            )
            .sort_values("Specificity", ascending=False)
            .reset_index(drop=True)
        )
        st.dataframe(rank_df, use_container_width=True)
    else:
        st.info("Run off-target scan to get specificity ranking.")
