
🧬 CRISPR Guide RNA Designer

A fast, user-friendly web app for designing, scoring, and analyzing CRISPR guide RNAs (gRNAs).

🚀 Try It Now 👉 Launch the app on Streamlit 👈


🎯 Unique Features & USP

Simple: Paste a DNA sequence or upload a FASTA file—no coding or setup required.

Smart Scoring: Each guide is scored with both lab rules and data-inspired heuristics (Hybrid, ML-inspired, and Consensus scores).

Off-Target Awareness: Scan your custom background DNA for potential off-targets.

AI Summary: One-click Gemini/OpenAI-powered guide report.

Free & Open Source: Use, share, or modify—no restrictions.


👥 Who is this for? (ICP)

Molecular biology, plant, and biomedical researchers

Academic labs and teaching settings

Biotech startups and R&D teams

Students, beginners, and DIY bio enthusiasts



🛠️ How to Use

1. Open the Streamlit app.

2. Paste your DNA sequence or upload a FASTA file.

3. Set CRISPR/PAM parameters as needed.

4. Click Find gRNAs to view candidate guides.

5. Review scores and download your results.

6. (Optional) Provide background DNA to scan for off-targets.

7. (Optional) Simulate edits or request an AI-powered summary.


📊 Scoring Methodology

Our gRNA scoring combines rule-based and machine learning-inspired approaches to prioritize reliable guides:

Hybrid Score (Rule-based):
Calculated from biochemical rules including GC content (optimal 40–70%), homopolymer avoidance, seed region composition, off-target penalties, and terminal base preferences. These are grounded in established CRISPR research and provide interpretable, experimentally validated guidance.

ML-inspired Score:
Based on features identified in large CRISPR screening datasets using machine learning. Although no trained ML model is directly used here, this score mimics key predictive factors like nucleotide composition and positional effects derived from published studies.

Consensus Score:
The average of Hybrid and ML-inspired scores, offering a balanced metric that combines experimental rules with data-driven insights for better guide prioritization.


> Note: These scores assist in guide selection but do not replace experimental validation, which is essential for confirming guide efficiency and specificity.


📝 Installation (for local use)

git clone https://github.com/Dinesh431786/Crispr.git
cd crispr/crispr_app
pip install -r requirements.txt
streamlit run app.py

🔑 License

MIT License – Free for all use.

🙋 FAQ

Is this tool free to use?
Yes, it’s open source and free for all users.

Does this use actual machine learning?
No, scores are rule-based and inspired by ML findings, not derived from trained models.

Is this suitable for crops/animals/bacteria?
Yes, it works for any DNA input, but always validate experimentally.

🤝 Contributing

Pull requests, suggestions, and issue reports are welcome! See CONTRIBUTING.md for guidelines.

👨‍🔬 Author

Dinesh K – design, code, and documentation

