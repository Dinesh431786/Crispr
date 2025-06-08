import streamlit as st
import pandas as pd
from utils import (
    validate_sequence,
    load_fasta,
    visualize_guide_location,
    plot_protein_domains
)
from analysis import (
    find_gRNAs,
    find_off_targets_detailed,
    simulate_protein_edit,
    diff_proteins,
    indel_simulations,
    hybrid_score,
    ml_gRNA_score,        # NEW: real ML model (even simple RF for demo)
    predict_hdr_repair,   # NEW: repair logic
    annotate_protein_domains,
)
import datetime

# --- Gemini imports
import google.generativeai as genai

st.set_page_config(page_title="🧬 CRISPR Lab NextGen", layout="wide")
st.title("🧬 CRISPR Lab NextGen – Advanced AI-Powered gRNA Designer")

# Sidebar – Input
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
    pam = {"Cas9 NGG": "NGG", "Cas9 NAG": "NAG", "Cas12a TTTV": "TTTV"}[pam_label]
    guide_len = st.slider("Guide length", 18, 25, 20, key="guide_len")
    min_gc = st.slider("Min GC %", 30, 60, 40, key="min_gc")
    max_gc = st.slider("Max GC %", 60, 80, 70, key="max_gc")
    bg_seq = st.text_area("Background DNA (off-target, optional)", height=100, key="bg_seq")
    max_mm = st.slider("Max mismatches", 0, 4, 2, key="max_mm")
    edit_offset = st.slider(
        "Edit offset from PAM", 0, guide_len, guide_len, key="edit_offset",
        help="Cas9 cut ≈ 3 bp upstream of PAM; set as needed."
    )
    st.header("🤖 Gemini AI")
    gemini_key = st.text_input("Gemini API Key", type="password", key="gemini_api_key")
    if gemini_key and len(gemini_key.strip()) > 10:
        st.success("Gemini API initialized!")

for k in (
    "df_guides", "offtargets", "guide_scores", "selected_gRNA", "selected_edit",
    "sim_result", "sim_indel", "protein_domains", "offtarget_counts", "domain_penalties",
    "gemini_summary"
):
    st.session_state.setdefault(k, None)

# gRNA search
if st.button("🔍 Find gRNAs"):
    ok, seq_or_msg = validate_sequence(dna_seq)
    if not ok:
        st.error(seq_or_msg)
        st.session_state.df_guides = None
    else:
        with st.spinner("Searching gRNAs…"):
            guides = find_gRNAs(seq_or_msg, pam, guide_len, min_gc, max_gc)
            guides["OffTargetCount"] = 0
            guides["DomainPenalty"] = 0.0
            guides["HybridScore"] = 0.0
            guides["MLScore"] = 0.0
            guides["HDRType"] = ""
            st.session_state.offtarget_counts = {}
            st.session_state.domain_penalties = {}

            # --- OFF-TARGET
            if bg_seq and len(bg_seq.strip()) >= guide_len + 3:
                off_targets_df = find_off_targets_detailed(guides, bg_seq, max_mm)
                st.session_state.offtargets = off_targets_df
                for g in guides.gRNA:
                    count = len(off_targets_df[off_targets_df.gRNA == g])
                    guides.loc[guides.gRNA == g, "OffTargetCount"] = count
                    st.session_state.offtarget_counts[g] = count
            else:
                st.session_state.offtargets = None

            # --- PROTEIN DOMAIN PENALTY (HMMER, local)
            domain_df = None
            if len(guides) > 0:
                domain_df = annotate_protein_domains(seq_or_msg)
                st.session_state.protein_domains = domain_df
                for idx, row in guides.iterrows():
                    cut_pos = row["Start"] + edit_offset
                    domain_penalty = 0.0
                    if domain_df is not None and not domain_df.empty and cut_pos:
                        for _, d in domain_df.iterrows():
                            if d["StartAA"] <= cut_pos//3 <= d["EndAA"]:
                                domain_penalty = 1.0
                    guides.at[idx, "DomainPenalty"] = domain_penalty
                    st.session_state.domain_penalties[row["gRNA"]] = domain_penalty

            # --- TRUE GUIDE SCORING: Hybrid + ML + Penalties
            for idx, row in guides.iterrows():
                guides.at[idx, "MLScore"] = ml_gRNA_score(row["gRNA"])
                guides.at[idx, "HybridScore"] = hybrid_score(
                    row["gRNA"],
                    off_target_count=row["OffTargetCount"],
                    domain_penalty=row["DomainPenalty"]
                )
                # --- REPAIR LOGIC
                cut_pos = row["Start"] + edit_offset
                guides.at[idx, "HDRType"] = predict_hdr_repair(dna_seq, cut_pos)
            st.session_state.df_guides = guides
        st.session_state.update(
            guide_scores=None, sim_result=None, sim_indel=None, gemini_summary=None,
        )

df = st.session_state.df_guides
if df is None or df.empty:
    st.info("Paste DNA & click **Find gRNAs** to begin.")
    st.stop()

st.success(f"✅ {len(df)} gRNAs found. Showing top-scoring guides.")
st.dataframe(
    df.sort_values("HybridScore", ascending=False),
    use_container_width=True,
    hide_index=True
)
st.download_button("⬇️ Download gRNAs CSV", df.to_csv(index=False), "guides.csv")

tab_ot, tab_sim, tab_vis, tab_rank, tab_gemini, tab_report = st.tabs(
    ["Off-targets", "Simulation & Indel", "Visualization", "Ranking", "Gemini AI Summary", "Project Report"]
)

with tab_ot:
    if st.session_state.offtargets is None or st.session_state.offtargets.empty:
        st.info("No off-targets or background DNA provided.")
    else:
        st.dataframe(st.session_state.offtargets, use_container_width=True)
        st.download_button(
            "⬇️ Download off-targets",
            st.session_state.offtargets.to_csv(index=False),
            "offtargets.csv"
        )

EDIT_TYPES = {
    "Delete 1 bp": "del1",
    "Insert A": "insA",
    "Delete 3 bp": "del3",
    "Insert G": "insG",
    "Substitute A→T": "subAG",
}
with tab_sim:
    g_list = df.gRNA.tolist()
    st.session_state.selected_gRNA = st.selectbox("gRNA", g_list, key="sel_gRNA")
    edit_label = st.selectbox(
        "Edit type", list(EDIT_TYPES.keys()), key="edit_label"
    )
    st.session_state.selected_edit = EDIT_TYPES[edit_label]
    sub_from = sub_to = ""
    if st.session_state.selected_edit == "subAG":
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
                st.session_state.selected_edit,
                sub_from=sub_from,
                sub_to=sub_to,
            )
            st.session_state.sim_indel = indel_simulations(
                dna_seq, idx + edit_offset
            )
            st.session_state.protein_domains = annotate_protein_domains(dna_seq)

    if st.session_state.sim_result:
        before, after, fs, stop = st.session_state.sim_result
        st.markdown(f"**Before protein:** `{before}`")
        st.markdown(f"**After protein:** `{after}`")
        st.markdown(f"**Diff:** {diff_proteins(before, after)}")
        st.write("Frameshift:", fs, "| Premature stop:", stop)
        st.markdown(f"**Repair Type:** `{predict_hdr_repair(dna_seq, idx + edit_offset)}`")
    if st.session_state.sim_indel is not None:
        st.subheader("±1–3 bp indel simulation")
        st.dataframe(st.session_state.sim_indel, use_container_width=True)
    if st.session_state.protein_domains is not None:
        st.subheader("Protein Domain Annotation (local, real)")
        st.dataframe(st.session_state.protein_domains)
        # Show domain/cut viz
        if not st.session_state.protein_domains.empty:
            cut_aa = idx // 3 if idx else 0
            fig = plot_protein_domains(
                before, st.session_state.protein_domains, cut_aa
            )
            st.plotly_chart(fig, use_container_width=True)

with tab_vis:
    idx = dna_seq.upper().find(st.session_state.selected_gRNA)
    if idx != -1:
        ax = visualize_guide_location(
            dna_seq,
            st.session_state.selected_gRNA,
            idx
        )
        st.pyplot(ax.figure)
    else:
        st.info("gRNA not found for visualization.")

with tab_rank:
    if "HybridScore" in df.columns:
        rank_df = (
            df[["gRNA", "HybridScore", "MLScore", "OffTargetCount", "DomainPenalty", "HDRType"]]
            .sort_values("HybridScore", ascending=False)
            .reset_index(drop=True)
        )
        st.dataframe(rank_df, use_container_width=True)
    else:
        st.info("Run gRNA search to get ranking.")

with tab_gemini:
    st.header("Gemini AI Summary – Smart Project Insights")
    if st.button("Get Gemini AI Summary"):
        # Build rich context
        gRNA = st.session_state.selected_gRNA
        guide_row = df[df["gRNA"] == gRNA].iloc[0]
        sim_result = st.session_state.sim_result
        indel_df = st.session_state.sim_indel
        offtargets_df = st.session_state.offtargets
        domain_df = st.session_state.protein_domains

        prompt = f"""
        # CRISPR Edit Summary
        gRNA: {gRNA}
        HybridScore: {guide_row['HybridScore']}
        MLScore: {guide_row['MLScore']}
        OffTargetCount: {guide_row['OffTargetCount']}
        DomainPenalty: {guide_row['DomainPenalty']}
        HDRType: {guide_row['HDRType']}
        PAM: {guide_row['PAM']}
        Strand: {guide_row['Strand']}
        GC%: {guide_row['GC%']}
        Simulation: {sim_result}
        Protein Domains: {domain_df.to_dict(orient='records') if domain_df is not None else None}
        Indel Effects: {indel_df.to_dict(orient='records') if indel_df is not None else None}
        Off-targets: {offtargets_df[offtargets_df.gRNA == gRNA].to_dict(orient='records') if offtargets_df is not None else None}
        Generate a detailed, accurate summary including biological risks, likely effect, and experiment advice.
        """

        try:
            genai.configure(api_key=gemini_key)
            model = genai.GenerativeModel("models/gemini-1.5-flash-latest")
            response = model.generate_content(prompt)
            st.session_state.gemini_summary = response.text if hasattr(response, "text") else str(response)
        except Exception as e:
            st.session_state.gemini_summary = f"Gemini API error: {e}"

    if st.session_state.gemini_summary:
        st.info(st.session_state.gemini_summary)

with tab_report:
    st.header("📝 Project Report & Export")
    st.markdown(f"**Date:** {datetime.date.today()}")
    st.markdown("**Guide Table**")
    st.dataframe(df, use_container_width=True)
    st.markdown("**Protein Domains**")
    if st.session_state.protein_domains is not None:
        st.dataframe(st.session_state.protein_domains)
    st.download_button(
        "Download gRNA Table (CSV)",
        df.to_csv(index=False),
        "guides.csv"
    )
    if st.session_state.offtargets is not None:
        st.download_button(
            "Download Off-Target Table (CSV)",
            st.session_state.offtargets.to_csv(index=False),
            "offtargets.csv"
        )
