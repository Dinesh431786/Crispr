
import streamlit as st
import pandas as pd
from utils import validate_sequence, load_fasta
from analysis import (
    find_gRNAs, find_off_targets_detailed,
    simulate_protein_edit, diff_proteins, indel_simulations,
    hybrid_score, ml_gRNA_score,
)

st.set_page_config(page_title="🧬 CRISPR gRNA Designer", layout="wide")

steps = [
    "1️⃣ DNA Input",
    "2️⃣ gRNA Design",
    "3️⃣ Edit Simulation",
    "4️⃣ Off-targets & Ranking",
    "5️⃣ AI Explain / Gemini Report"
]
step = st.sidebar.radio("🚦 Workflow Steps", steps)

# --- Session State Defaults ---
for k in (
    "dna_seq", "df_guides", "offtargets", "guide_scores",
    "selected_gRNA", "sim_result", "sim_indel", "gemini_report", "ai_response"
):
    st.session_state.setdefault(k, None)

# ---------- STEP 1: DNA INPUT ----------
if step == steps[0]:
    st.title("🧬 Step 1: Upload or Paste Your DNA Sequence")
    st.info("DNA can be plant, animal, microbe, or synthetic—anything with A/T/C/G.")

    with st.expander("❓ What formats work?", expanded=False):
        st.write("Paste plain DNA or upload a FASTA/text file (single sequence per file).")

    uploaded = st.file_uploader("Upload .fasta or .txt", type=["fasta", "fa", "txt"])
    dna_seq = st.text_area("Or paste DNA sequence here:", height=120, key="input_dna_seq")

    if uploaded:
        seq, err = load_fasta(uploaded)
        if err:
            st.error(err)
        else:
            dna_seq = seq
            st.success("Loaded DNA from file!")
    st.session_state.dna_seq = (dna_seq or "").strip()

    if st.button("✅ Proceed to gRNA Design"):
        ok, msg = validate_sequence(st.session_state.dna_seq)
        if not ok:
            st.error(msg)
        else:
            st.success("DNA accepted! Select Step 2 (gRNA Design) in sidebar.")

# ---------- STEP 2: gRNA DESIGN ----------
elif step == steps[1]:
    st.title("🧬 Step 2: Design & Score gRNAs")
    if not st.session_state.dna_seq:
        st.warning("Please upload or paste a DNA sequence in Step 1.")
        st.stop()

    st.markdown(f"""
        <div style="background-color:#f4fbfa; border-radius:8px; padding:15px;">
        <b>Design Summary:</b><br>
        <span style="color:#3475b5;">DNA length:</span> {len(st.session_state.dna_seq)} bp<br>
        <span style="color:#5ac980;">Recommended guide length:</span> 20 bp
        </div>
        """, unsafe_allow_html=True)

    with st.expander("⚙️ Show Advanced Settings", expanded=False):
        pam_label = st.selectbox("PAM", ["Cas9 NGG", "Cas9 NAG", "Cas9 NG", "Cas12a TTTV"], key="pam")
        GUIDE_TYPES = {"Cas9 NGG": "NGG", "Cas9 NAG": "NAG", "Cas9 NG": "NG", "Cas12a TTTV": "TTTV"}
        pam = GUIDE_TYPES[pam_label]
        u6_toggle = st.toggle("U6 Promoter (add G at 5’ if needed)", value=False, key="u6_toggle")
        guide_len = st.slider("Guide length", 18, 25, 20, key="guide_len")
        min_gc = st.slider("Min GC %", 30, 60, 40, key="min_gc")
        max_gc = st.slider("Max GC %", 60, 80, 70, key="max_gc")

    st.caption("**Hybrid Score:** Lab rules; **ML Score:** ML-inspired rules; **Consensus:** Average of both.")

    if st.button("🔍 Find gRNAs", type="primary"):
        ok, msg = validate_sequence(st.session_state.dna_seq)
        if not ok:
            st.error(msg)
            st.session_state.df_guides = None
        else:
            with st.spinner("Finding gRNAs..."):
                df = find_gRNAs(
                    st.session_state.dna_seq, pam, guide_len, min_gc, max_gc
                )
                if df.empty:
                    st.error("No gRNAs found for these settings.")
                else:
                    df["HybridScore"] = [hybrid_score(g) for g in df.gRNA]
                    df["MLScore"] = [ml_gRNA_score(g) for g in df.gRNA]
                    df["ConsensusScore"] = ((df["HybridScore"] + df["MLScore"]) / 2).clip(upper=1.0)
                    if u6_toggle:
                        df["gRNA"] = df["gRNA"].apply(lambda g: g if g.startswith("G") else "G" + g[:-1])
                    st.session_state.df_guides = df
                    st.success(f"Found {len(df)} gRNAs! See results below.")

    df = st.session_state.df_guides
    if df is not None and not df.empty:
        best = df[df.ConsensusScore >= 0.8]
        st.subheader("🌟 Top gRNA Candidates")
        st.dataframe(best[["Strand", "Start", "gRNA", "PAM", "GC%", "ConsensusScore"]], use_container_width=True)
        st.download_button("⬇️ Download all gRNAs as CSV", df.to_csv(index=False), "guides.csv")
        with st.expander("🔬 Show All Guides", expanded=False):
            st.dataframe(df, use_container_width=True)

# ---------- STEP 3: EDIT SIMULATION ----------
elif step == steps[2]:
    st.title("🧬 Step 3: Simulate Edit / Indel Effects")
    df = st.session_state.df_guides
    if df is None or df.empty:
        st.warning("Please find gRNAs first (go back to Step 2).")
        st.stop()
    g_list = df.gRNA.tolist()
    if not g_list:
        st.warning("No gRNAs available.")
        st.stop()
    st.session_state.selected_gRNA = st.selectbox("Choose a gRNA for simulation", g_list, key="sel_gRNA_sim")
    edit_options = {
        "Delete 1 bp": "del1", "Insert A": "insA",
        "Delete 3 bp": "del3", "Insert G": "insG", "Substitute A→T": "subAG"
    }
    st.session_state.selected_edit = st.selectbox("Edit type", list(edit_options), key="edit_type")
    sub_from = sub_to = ""
    if edit_options[st.session_state.selected_edit] == "subAG":
        sub_from = st.text_input("Sub FROM", "A")
        sub_to = st.text_input("Sub TO", "T")
    guide_len = st.session_state.df_guides["gRNA"].str.len().iloc[0]
    edit_offset = st.slider("Edit offset from PAM", 0, guide_len, guide_len, key="edit_offset")

    if st.button("Simulate Edit"):
        idx = (st.session_state.dna_seq or "").upper().find(st.session_state.selected_gRNA)
        if idx == -1:
            st.error("gRNA not found in sequence! Try another guide or check your DNA.")
        else:
            st.session_state.sim_result = simulate_protein_edit(
                st.session_state.dna_seq, idx + edit_offset,
                edit_options[st.session_state.selected_edit], sub_from=sub_from, sub_to=sub_to
            )
            st.session_state.sim_indel = indel_simulations(
                st.session_state.dna_seq, idx + edit_offset
            )
    # Show result visually
    if st.session_state.sim_result:
        before, after, fs, stop = st.session_state.sim_result
        st.markdown(f"<b>Before protein:</b> <span style='color:#238b21'>{before}</span>", unsafe_allow_html=True)
        st.markdown(f"<b>After protein:</b> <span style='color:#2239d7'>{after}</span>", unsafe_allow_html=True)
        st.markdown(f"<b>Diff:</b> {diff_proteins(before, after)}")
        st.write("Frameshift:", fs, "| Premature stop:", stop)
    if st.session_state.sim_indel is not None:
        st.subheader("±1–3 bp indel simulation")
        st.dataframe(st.session_state.sim_indel, use_container_width=True)

# ---------- STEP 4: OFF-TARGETS & RANKING ----------
elif step == steps[3]:
    st.title("🧬 Step 4: Off-target Analysis & Specificity Ranking")
    df = st.session_state.df_guides
    if df is None or df.empty:
        st.warning("Design gRNAs first! Go back to Step 2.")
        st.stop()
    with st.expander("⚙️ Advanced Off-target Settings", expanded=False):
        bg_seq = st.text_area("Background DNA (for off-target search)", key="bg_seq")
        max_mm = st.slider("Max mismatches", 0, 4, 2, key="max_mm")
    if st.button("🔬 Scan off-targets"):
        if not bg_seq or len(bg_seq.strip()) < 10:
            st.error("Please paste a realistic background DNA (at least 10 nt)!")
            st.session_state.offtargets = pd.DataFrame()
        else:
            result_from_find = find_off_targets_detailed(df, bg_seq, max_mm)
            st.session_state.offtargets = result_from_find
            # Score
            scores = {}
            if not result_from_find.empty and "gRNA" in result_from_find.columns and "Mismatches" in result_from_find.columns:
                for g in df.gRNA:
                    subset = result_from_find[result_from_find["gRNA"] == g]
                    if subset.empty:
                        scores[g] = 1.0
                    else:
                        scores[g] = round(1.0 / (1 + subset["Mismatches"].sum()), 3)
            else:
                scores = {g: 1.0 for g in df.gRNA}
            st.session_state.guide_scores = scores

    ot_df = st.session_state.offtargets
    if ot_df is not None:
        st.subheader("Off-target Results")
        if ot_df.empty:
            st.info("No off-targets found! Paste a gRNA (or variant) as background to test matching.")
        else:
            st.dataframe(ot_df, use_container_width=True)
            st.download_button("⬇️ Download off-targets", ot_df.to_csv(index=False), "offtargets.csv")
    # Specificity ranking
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
        st.subheader("gRNA Specificity Ranking")
        st.dataframe(rank_df, use_container_width=True)
    else:
        st.info("Run off-target scan to get specificity ranking.")

# ---------- STEP 5: AI EXPLAIN / GEMINI REPORT ----------
elif step == steps[4]:
    st.title("🧬 Step 5: AI Explain & Gemini Report")
    st.info("Ask Gemini/OpenAI for a scientific summary or experiment advice! Enter your API key below.")
    with st.expander("⚙️ AI Backend Settings", expanded=False):
        ai_backend = st.selectbox("AI Backend", ["Gemini", "OpenAI"], key="ai_backend")
        gemini_model = "gemini-1.5-flash-latest"
        if ai_backend == "Gemini":
            gemini_model = st.selectbox(
                "Gemini Model", [
                    "gemini-1.5-flash-latest", "gemini-1.5-pro-latest",
                    "gemini-pro", "gemini-1.0-pro-latest"
                ], key="gemini_model"
            )
        api_key = st.text_input("API Key", type="password", key="api_key")
    # Prompt for AI
    prompt = st.text_area("AI prompt (optional):", "Summarize gRNA results and potential experimental outcomes.")
    if st.button("📄 Generate AI/Gemini Report"):
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
            if "API key not valid" in error_str or "API_KEY_INVALID" in error_str:
                st.error("❌ Your API key is invalid or not enabled for this model.")
            elif "model not found" in error_str or "not supported" in error_str:
                st.error("❌ The selected model is not available. Try another one in settings.")
            else:
                st.error(f"API error: {error_str}")
            st.session_state.gemini_report = ""
    if st.session_state.gemini_report:
        st.subheader("Gemini/OpenAI Report")
        st.info(st.session_state.gemini_report)

# ---- Minor styling: improve card/section visuals ----
st.markdown("""
<style>
.stButton>button {font-size: 1.1rem; padding: 0.6em 1.2em;}
.stDataFrame, .stTable {background: #fffaf4;}
.stRadio>div>label {font-weight: bold;}
.stExpander {background: #f7fafd !important;}
[data-testid="stSidebar"] {background-color: #f4fbfa;}
</style>
""", unsafe_allow_html=True)