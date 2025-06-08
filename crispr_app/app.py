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
)

# ──────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────
GUIDE_TYPES = {
    "Cas9 NGG": "NGG",
    "Cas9 NAG": "NAG",
    "Cas12a TTTV": "TTTV",
}
EDIT_TYPES = {
    "Delete 1 bp": "del1",
    "Insert A": "insA",
    "Delete 3 bp": "del3",
    "Insert G": "insG",
    "Substitute A→T": "subAG",
}

st.set_page_config(page_title="🧬 CRISPR Lab NextGen", layout="wide")
st.title("🧬 CRISPR Lab NextGen – gRNA Designer & Impact Analyzer")

# ──────────────────────────────────────────────────────────
# Sidebar – sequence & AI settings
# ──────────────────────────────────────────────────────────
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

    pam_label = st.selectbox("PAM", list(GUIDE_TYPES.keys()), key="pam")
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

    # AI settings
    st.header("🔑 AI Settings")
    st.session_state.ai_backend = st.selectbox(
        "Backend", ["Gemini", "OpenAI"], key="ai_backend"
    )
    st.session_state.api_key = st.text_input(
        "API Key", type="password", key="api_key"
    )
    # Instant confirmation
    if st.session_state.api_key and len(st.session_state.api_key.strip()) > 10:
        st.success(f"{st.session_state.ai_backend} API initialized!")

# ──────────────────────────────────────────────────────────
# Initialise session-state holders
# ──────────────────────────────────────────────────────────
for k in (
    "df_guides",
    "offtargets",
    "guide_scores",
    "selected_gRNA",
    "selected_edit",
    "sim_result",
    "sim_indel",
    "ai_response",
):
    st.session_state.setdefault(k, None)

# ──────────────────────────────────────────────────────────
# gRNA search
# ──────────────────────────────────────────────────────────
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
        # reset downstream state
        st.session_state.update(
            offtargets=None,
            guide_scores=None,
            sim_result=None,
            sim_indel=None,
            ai_response="",
        )

df = st.session_state.df_guides
if df is None or df.empty:
    st.info("Paste DNA & click **Find gRNAs** to begin.")
    st.stop()

# ──────────────────────────────────────────────────────────
# Display gRNAs
# ──────────────────────────────────────────────────────────
st.success(f"✅ {len(df)} gRNAs found")
st.dataframe(df, use_container_width=True)
st.download_button("⬇️ Download gRNAs CSV", df.to_csv(index=False), "guides.csv")

tab_ot, tab_sim, tab_ai, tab_vis, tab_rank = st.tabs(
    ["Off-targets", "Simulation & Indel", "AI Explain", "Visualization", "Ranking"]
)

# ──────────────────────────────────────────────────────────
# Off-target tab
# ──────────────────────────────────────────────────────────
with tab_ot:
    if not bg_seq.strip():
        st.info("Provide background DNA in sidebar for off-target scanning.")
    else:
        if st.button("Scan off-targets"):
            st.session_state.offtargets = find_off_targets_detailed(
                df, bg_seq, max_mm
            )
            # simple specificity score
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

# ──────────────────────────────────────────────────────────
# Simulation & Indel tab
# ──────────────────────────────────────────────────────────
with tab_sim:
    g_list = df.gRNA.tolist()
    st.session_state.selected_gRNA = st.selectbox(
        "gRNA", g_list, key="sel_gRNA"
    )
    st.session_state.selected_edit = st.selectbox(
        "Edit type", list(EDIT_TYPES), key="sel_edit"
    )

    # extra fields for substitution
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

# ──────────────────────────────────────────────────────────
# AI Explain tab
# ──────────────────────────────────────────────────────────
with tab_ai:
    ctx = st.text_area(
        "Describe the edit context",
        f"Editing at {st.session_state.selected_gRNA}",
        key="ai_ctx",
    )
    if st.button("Ask AI"):
        key = st.session_state.api_key
        if not key or len(key.strip()) < 10:
            st.error("Enter a valid API key in sidebar.")
        else:
            try:
                if st.session_state.ai_backend == "Gemini":
                    import google.generativeai as genai

                    genai.configure(api_key=key)
                    model = genai.GenerativeModel("gemini-pro")
                    st.session_state.ai_response = model.generate_content(ctx).text
                else:  # OpenAI
                    import openai

                    openai.api_key = key
                    resp = openai.ChatCompletion.create(
                        model="gpt-3.5-turbo",
                        messages=[
                            {"role": "system", "content": "You are a bioinformatics expert."},
                            {"role": "user", "content": ctx},
                        ],
                    )
                    st.session_state.ai_response = resp.choices[0].message.content
            except Exception as e:
                import traceback

                st.session_state.ai_response = (
                    "API error:\n" + traceback.format_exc(limit=2)
                )
    if st.session_state.ai_response:
        st.info(st.session_state.ai_response)

# ──────────────────────────────────────────────────────────
# Visualization tab
# ──────────────────────────────────────────────────────────
with tab_vis:
    idx = dna_seq.upper().find(st.session_state.selected_gRNA)
    if idx != -1:
        ax = visualize_guide_location(dna_seq, st.session_state.selected_gRNA, idx)
        st.pyplot(ax.figure)
    else:
        st.info("gRNA not found for visualization.")

# ──────────────────────────────────────────────────────────
# Ranking tab
# ──────────────────────────────────────────────────────────
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
