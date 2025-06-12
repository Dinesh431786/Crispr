import streamlit as st
import pandas as pd
from utils import validate_sequence, load_fasta
from analysis import (
    find_gRNAs, find_off_targets_detailed, simulate_protein_edit, diff_proteins,
    indel_simulations, hybrid_score, ml_gRNA_score
)

# ---- App Config ----
st.set_page_config(page_title="🧬 CRISPR Guide RNA Designer", layout="wide")

# ---- HEADER ----
st.title("🧬 CRISPR Guide RNA Designer")
st.markdown("#### Fast, simple tool for **plants, humans, microbes, and any DNA!**")

# ---- SIDEBAR: Sequence Upload ----
with st.sidebar:
    st.header("1️⃣ Upload Sequence")
    uploaded = st.file_uploader("Upload .fasta", type=["fasta", "fa", "txt"])
    dna_seq = st.text_area("Or paste DNA sequence:", height=150, key="dna_seq")
    if uploaded:
        seq, err = load_fasta(uploaded)
        if err:
            st.error(err)
        else:
            dna_seq = seq

    pam_label = st.selectbox("PAM", ["Cas9 NGG", "Cas9 NAG", "Cas9 NG", "Cas12a TTTV"], key="pam")
    GUIDE_TYPES = {"Cas9 NGG": "NGG", "Cas9 NAG": "NAG", "Cas9 NG": "NG", "Cas12a TTTV": "TTTV"}
    pam = GUIDE_TYPES[pam_label]
    u6_g_toggle = st.toggle("U6 Promoter (add G at 5’ if needed)", value=False)
    guide_len = st.slider("Guide length", 18, 25, 20)
    min_gc = st.slider("Min GC %", 30, 60, 40)
    max_gc = st.slider("Max GC %", 60, 80, 70)

# ---- MAIN ----

# Session state for data
for k in ("df_guides", "offtargets", "guide_scores", "sim_result", "sim_indel", "ai_response", "gemini_report"):
    st.session_state.setdefault(k, None)

# ---- STEP 1: Find gRNAs ----
st.header("Step 1: Guide RNA Discovery")
st.info("**Hybrid Score** = lab rules; **ML Score** = pattern-based. Aim for Consensus >0.8.")

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
        st.session_state.update(offtargets=None, guide_scores=None, sim_result=None, sim_indel=None, ai_response="", gemini_report=None)

df = st.session_state.df_guides
if df is None or df.empty:
    st.info("Paste DNA & click **Find gRNAs** to begin.")
    st.stop()

if "HybridScore" not in df.columns or "MLScore" not in df.columns or "ConsensusScore" not in df.columns:
    df["HybridScore"] = [hybrid_score(g) for g in df.gRNA]
    df["MLScore"] = [ml_gRNA_score(g) for g in df.gRNA]
    df["ConsensusScore"] = ((df["HybridScore"] + df["MLScore"]) / 2).clip(upper=1.0)

def u6_g_mod(seq):
    if seq.startswith("G"):
        return seq
    return "G" + seq

def apply_u6_toggle_to_df(df, u6_toggle):
    df_mod = df.copy()
    if u6_toggle and "gRNA" in df_mod.columns:
        df_mod["gRNA"] = df_mod["gRNA"].apply(u6_g_mod)
    return df_mod

u6_toggle = u6_g_toggle
df_display = apply_u6_toggle_to_df(df, u6_toggle)

st.success(f"✅ {len(df_display)} gRNAs found")
st.dataframe(df_display, use_container_width=True)
st.download_button("⬇️ Download gRNAs CSV", df_display.to_csv(index=False), "guides.csv")

# ---- STEP 2: Off-Target Analysis & Specificity Ranking ----
st.header("Step 2: Off-Target Analysis & Specificity Ranking")
bg_seq = st.text_area("Background DNA (off-target)", height=100, key="bg_seq")
max_mm = st.slider("Max mismatches", 0, 4, 2, key="max_mm_slider")

if st.button("Scan off-targets"):
    result_from_find = find_off_targets_detailed(df, bg_seq, max_mm)
    # handle series case
    if isinstance(result_from_find, pd.Series):
        ot_df = result_from_find.to_frame().T
    else:
        ot_df = result_from_find
    st.session_state.offtargets = ot_df

    # Specificity score
    scores = {}
    if ot_df is not None and not ot_df.empty and "gRNA" in ot_df.columns and "Mismatches" in ot_df.columns:
        for g in df.gRNA:
            subset = ot_df[ot_df["gRNA"] == g]
            if subset.empty:
                scores[g] = 1.0
            else:
                scores[g] = round(1.0 / (1 + subset["Mismatches"].sum()), 3)
    else:
        scores = {g: 1.0 for g in df.gRNA}
        if ot_df is not None and not ot_df.empty:
            st.error("Off-target results missing required columns ('gRNA', 'Mismatches').")
    st.session_state.guide_scores = scores

ot_df = st.session_state.offtargets
if ot_df is not None:
    if ot_df.empty:
        st.info("No off-targets within given mismatches.")
    else:
        st.dataframe(ot_df, use_container_width=True)
        st.download_button("⬇️ Download off-targets", ot_df.to_csv(index=False), "offtargets.csv")

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
    st.subheader("gRNA Specificity Ranking")
    st.dataframe(rank_df, use_container_width=True)
else:
    st.info("Run off-target scan to get specificity ranking.")

# ---- STEP 3: Edit/Indel Simulation ----
st.header("Step 3: Simulate Edit/Indel Effects")

if u6_toggle:
    g_list_display = [u6_g_mod(g) for g in df.gRNA.tolist()]
else:
    g_list_display = df.gRNA.tolist()
gRNA_display_to_seq = {u6_g_mod(g) if u6_toggle else g: g for g in df.gRNA.tolist()}

selected_gRNA = st.selectbox("Choose a gRNA for simulation", g_list_display, key="sel_gRNA")
gRNA_for_analysis = gRNA_display_to_seq[selected_gRNA]

EDIT_TYPES = {
    "Delete 1 bp": "del1",
    "Insert A": "insA",
    "Delete 3 bp": "del3",
    "Insert G": "insG",
    "Substitute A→T": "subAG",
}
selected_edit = st.selectbox("Edit type", list(EDIT_TYPES), key="sel_edit")
edit_offset = st.slider("Edit offset from PAM", 0, guide_len, guide_len, key="edit_offset")
sub_from = sub_to = ""
if EDIT_TYPES[selected_edit] == "subAG":
    sub_from = st.text_input("Sub FROM", "A")
    sub_to = st.text_input("Sub TO", "T")

if st.button("Simulate Edit"):
    idx = dna_seq.upper().find(gRNA_for_analysis)
    if idx == -1:
        st.error("gRNA not found in sequence!")
    else:
        st.session_state.sim_result = simulate_protein_edit(
            dna_seq, idx + edit_offset, EDIT_TYPES[selected_edit], sub_from=sub_from, sub_to=sub_to,
        )
        st.session_state.sim_indel = indel_simulations(dna_seq, idx + edit_offset)

if st.session_state.sim_result:
    before, after, fs, stop = st.session_state.sim_result
    st.markdown(f"**Before protein:** <span style='color:green'>{before}</span>", unsafe_allow_html=True)
    st.markdown(f"**After protein:** <span style='color:blue'>{after}</span>", unsafe_allow_html=True)
    st.markdown(f"<b>Diff:</b> {diff_proteins(before, after)}", unsafe_allow_html=True)
    st.write("Frameshift:", fs, "| Premature stop:", stop)
if st.session_state.sim_indel is not None:
    st.markdown("**±1–3 bp indel simulation**")
    st.dataframe(st.session_state.sim_indel, use_container_width=True)

# ---- STEP 4: AI Report (Gemini/OpenAI) ----
st.header("Step 4: Scientific AI Summary (Gemini/OpenAI)")

with st.expander("AI Backend Settings", expanded=False):
    ai_backend = st.selectbox("AI Backend", ["Gemini", "OpenAI"])
    gemini_model = st.selectbox("Gemini Model", [
        "gemini-1.5-flash-latest",
        "gemini-1.5-pro-latest",
        "gemini-pro",
        "gemini-1.0-pro-latest",
    ]) if ai_backend == "Gemini" else ""
    api_key = st.text_input("API Key", type="password")
    ai_prompt = st.text_area("AI prompt (optional):", "Summarize gRNA results and potential experimental outcomes.")

# -- Compose prompt with REAL data for AI --
context_parts = []

if df is not None and not df.empty:
    context_parts.append("## Top gRNAs and scores:\n")
    context_parts.append(df[["gRNA", "ConsensusScore", "HybridScore", "MLScore"]].head(10).to_csv(index=False))
if ot_df is not None and not ot_df.empty:
    context_parts.append("\n## Off-target summary (first 10):\n")
    context_parts.append(ot_df.head(10).to_csv(index=False))
if st.session_state.sim_result:
    before, after, fs, stop = st.session_state.sim_result
    context_parts.append("\n## Edit simulation result:\n")
    context_parts.append(f"Before protein: {before}\nAfter protein: {after}\nFrameshift: {fs}, Premature stop: {stop}")
if ai_prompt.strip():
    context_parts.append(f"\nUser notes: {ai_prompt.strip()}")

full_prompt = "\n".join(context_parts) + "\n\nSummarize the CRISPR results above for a scientist. Highlight best guides, off-target risks, and predicted edit effects."

with st.expander("See prompt/context sent to AI", expanded=False):
    st.code(full_prompt)

if st.button("Generate AI/Gemini Report"):
    if not api_key or len(api_key.strip()) < 10:
        st.error("Enter a valid API key in the sidebar.")
    else:
        try:
            with st.spinner("Contacting Gemini/OpenAI..."):
                if ai_backend == "Gemini":
                    import google.generativeai as genai
                    genai.configure(api_key=api_key)
                    model = genai.GenerativeModel(gemini_model)
                    result = model.generate_content(full_prompt)
                    st.session_state.gemini_report = result.text if hasattr(result, "text") else str(result)
                else:
                    import openai
                    openai.api_key = api_key
                    resp = openai.ChatCompletion.create(
                        model="gpt-3.5-turbo",
                        messages=[
                            {"role": "system", "content": "You are a CRISPR genome editing expert."},
                            {"role": "user", "content": full_prompt},
                        ],
                    )
                    st.session_state.gemini_report = resp.choices[0].message.content
        except Exception as e:
            error_str = str(e)
            if "API key not valid" in error_str or "API_KEY_INVALID" in error_str:
                st.error("❌ Your Gemini API key is invalid or this model is not enabled for your account/project. Please double-check your key and model selection.")
            elif "model not found" in error_str or "not supported" in error_str:
                st.error("❌ The selected Gemini model is not available. Try selecting 'gemini-pro' in the sidebar.")
            else:
                st.error(f"Gemini/OpenAI API error: {error_str}")
            st.session_state.gemini_report = ""

if st.session_state.gemini_report:
    st.subheader("Gemini/OpenAI Report")
    st.info(st.session_state.gemini_report)

st.markdown("---")
st.caption("Built by [YourName], 2024. For demo/research use only.")
