# 🧬 CRISPR Guide RNA Designer

A **fast, user-friendly web app** for designing, scoring, and analyzing CRISPR guide RNAs (gRNAs) — **no coding required**.

---

## 🚀 Try It Now

- [Launch the app on Streamlit](https://crispr-voxelta.streamlit.app/)  
- **Free. Open source. No login required.**

# 🧬 CRISPR Guide RNA Designer

**A fast, user-friendly web app for designing, scoring, and analyzing CRISPR guide RNAs (gRNAs).**



## 🎯 Unique Features

- **Zero setup:** Paste a DNA sequence or upload a FASTA file
- **Multiple PAMs:** Cas9 (NGG, NAG) & Cas12a (TTTV) supported
- **Hybrid & ML-inspired scoring:** Rule-based plus data-inspired consensus
- **Off-target scan:** Use any DNA as a custom background to spot risk sites
- **U6 promoter toggle:** One-click “add G at 5’” for U6/T7 promoters
- **Indel/protein simulation:** Visualize the effect of edits
- **AI reporting:** One-click Gemini/OpenAI-powered gRNA summaries
- **Download results:** CSV export for guides and off-targets
- **Modern UX:** Streamlit-based, works on desktop & mobile
- **MIT Licensed:** Use, share, fork, or modify

---

## 👥 Who is this for?

- Molecular, plant, or biomedical researchers
- Academic labs and classroom use
- Biotech & R&D teams (fast pilot CRISPR projects)
- Students, DIY bio, and open-science community

---

## 🛠️ How to Use

1. **Open the app** (see link above)
2. **Paste DNA** or **upload FASTA**
3. **Select PAM/parameters** in the sidebar
4. *(Optional)* Toggle **U6 Promoter** for gRNA with 5' G
5. Click **Find gRNAs**
6. Review and **download results**
7. *(Optional)* Paste background DNA to scan off-targets
8. *(Optional)* Run indel/protein simulation or AI-powered report

---

## 📊 Scoring Methodology

**Hybrid Score**: Based on established lab rules (GC content, homopolymers, seed region, off-target penalty, terminal base).  
**ML-inspired Score**: Derived from features found in ML studies of gRNA efficacy (but not a trained ML model).  
**Consensus Score**: Balanced average of the two for ranking.

> **Note:** Scores help prioritize guides, but do not replace experimental validation.

---

## 📝 Installation (For Local Use)

```bash
git clone https://github.com/Dinesh431786/Crispr.git
cd crispr/crispr_app
pip install -r requirements.txt
streamlit run app.py

🔑 AI API Keys
Gemini (Google): Get API key

OpenAI: Get API key

Paste your key in the app sidebar for AI-powered explanations.

🧪 Example FASTA
fasta
Copy
Edit
>TestGene
ATGAGTCTGCTCTTCGCGTTGGAGTGAAATCTGAGATGATGGGTTGAAATCGCAGTTCGACCTGAACTTTTATCTGCTCTTCGCGTTGAGCGGACCGTGGGAAGTTTCGCGTTGATCAGTTCTTCTGCTCTTCGCGTTTAAGCCTTGCGTTGTTTATCTGCTCTTCGCGTTTATCAGCCTGGGCGTTGATCTTTTATCTGCTCTTCGCGTTAACGGAAGCCGG
🙋 FAQ
Is it free?
Yes, open source and free for all.

Does it use real ML?
No, scoring is based on published ML findings but is rule-based.

Species?
Works for any DNA (human, plant, animal, microbe, synthetic).

AI required?
All core features work without AI. AI summary is optional with API key.

🤝 Contributing
Pull requests, feedback, and bug reports welcome!
See CONTRIBUTING.md for details.

👨‍🔬 Author
Dinesh K — design, code, and documentation
GitHub

⚖️ License
MIT License — free for all use.
