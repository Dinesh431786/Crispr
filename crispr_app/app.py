import streamlit as st
import subprocess
import tempfile
import pandas as pd
import shutil
from utils import validate_sequence, load_fasta, visualize_guide_location
from analysis import simulate_protein_edit, diff_proteins, find_off_targets

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

def ai_explain(edit_text, api_key, ai_backend):
    try:
        if ai_backend == "Gemini":
            import google.generativeai as genai
            genai.configure(api_key=api_key)
            model = genai.GenerativeModel('gemini-pro')
            response = model.generate_content(edit_text)
            return response.text
        elif ai_backend == "OpenAI":
            import openai
            openai.api_key = api_key
            completion = openai.ChatCompletion.create(
                model="gpt-4",
                messages=[{"role": "user", "content": edit_text}]
            )
            return completion['choices'][0]['message']['content']
    except Exception as e:
        return f"AI call failed: {e}"

st.set_page_config(page_title="🧬 CRISPR Lab", layout="wide")
st.title("🧬 CRISPR Lab: gRNA Designer & Impact Analyzer")

with st.sidebar:
    st.header("🧬 Sequence Input")
    uploaded = st.file_uploader("Upload .fasta file", type=["fasta", "fa", "txt"])
    dna_seq = st.text_area("Or paste DNA sequence:", height=150)

    if uploaded:
        loaded_seq, err = load_fasta(uploaded)
        if err:
            st.error(err)
        else:
            dna_seq = loaded_seq

    pam_label = st.selectbox("Select PAM Site", options=list(GUIDE_TYPES.keys()))
    pam = GUIDE_TYPES[pam_label]
    guide_length = st.slider("Guide Length", 18, 25, 20)
    bg_seq = st.text_area("Background DNA (for off-target)", height=100)

if st.button("🔍 Find gRNAs"):
    valid, msg = validate_sequence(dna_seq)
    if not valid:
        st.error(msg)
    elif shutil.which("guidemaker") is None:
        st.error("GuideMaker not found. Install it first.")
    else:
        with tempfile.NamedTemporaryFile("w+", suffix=".fasta", delete=False) as fasta_file:
            fasta_file.write(">input_seq\n" + dna_seq.upper().replace("\n", ""))
            fasta_file.flush()

            with tempfile.NamedTemporaryFile("w+", suffix=".csv", delete=False) as output_file:
                command = [
                    "guidemaker",
                    "--input_fasta", fasta_file.name,
                    "--pamseq", pam,
                    "--guidelength", str(guide_length),
                    "--pam_orientation", "3prime",
                    "--output", output_file.name
                ]

                try:
                    with st.spinner("Running GuideMaker..."):
                        subprocess.run(command, check=True)
                    df = pd.read_csv(output_file.name)
                    if df.empty:
                        st.warning("No valid gRNAs found. Try changing PAM or sequence.")
                    else:
                        st.success(f"✅ Found {len(df)} gRNAs")
                        st.dataframe(df)
                        csv = df.to_csv(index=False)
                        st.download_button("⬇️ Download gRNAs", data=csv, file_name="guides.csv", mime="text/csv")

                        tab1, tab2, tab3, tab4 = st.tabs(["🔍 Off-target", "🧬 Simulation", "🤖 AI Explain", "🖼️ Visualization"])

                        with tab1:
                            if bg_seq.strip():
                                ot_df = find_off_targets(df, bg_seq)
                                st.dataframe(ot_df)
                                st.download_button("⬇️ Download Off-target Results", data=ot_df.to_csv(index=False), file_name="offtarget.csv", mime="text/csv")
                            else:
                                st.info("Paste background DNA above to enable off-target search.")

                        with tab2:
                            gRNA_choice = st.selectbox("Choose gRNA", df["gRNA"].tolist())
                            cut_index = dna_seq.upper().find(gRNA_choice)
                            if cut_index != -1:
                                edit_label = st.selectbox("Edit Type", options=list(EDIT_TYPES.keys()))
                                edit_type = EDIT_TYPES[edit_label]
                                if edit_type == "subAG":
                                    sub_from = st.text_input("Substitute from", value="A")
                                    sub_to = st.text_input("To", value="T")
                                else:
                                    sub_from = sub_to = ""
                                prot_before, prot_after, fs, stop = simulate_protein_edit(
                                    dna_seq, cut_index + guide_length, edit_type, sub_from=sub_from, sub_to=sub_to
                                )
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
                                st.download_button("⬇️ Download Report", report, file_name="protein_report.txt")
                            else:
                                st.warning("Selected gRNA not found in sequence.")

                        with tab3:
                            ai_backend = st.selectbox("AI Backend", ["Gemini", "OpenAI"])
                            api_key = st.text_input("AI API Key", type="password")
                            gene_info = st.text_area("Describe the edit context", value=f"Editing at {gRNA_choice}")
                            if st.button("Ask AI"):
                                if api_key:
                                    explanation = ai_explain(gene_info, api_key, ai_backend)
                                    st.info(explanation)
                                else:
                                    st.warning("Enter API key first.")

                        with tab4:
                            ax = visualize_guide_location(dna_seq, gRNA_choice, cut_index)
                            st.pyplot(ax.figure)

                except subprocess.CalledProcessError as e:
                    st.error("GuideMaker failed to run.")
                    st.code(str(e))
