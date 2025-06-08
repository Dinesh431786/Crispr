🧬 CRISPR Guide RNA Designer
A fast, user-friendly web app for designing, scoring, and analyzing CRISPR guide RNAs (gRNAs).

🚀 Try It Now
👉 Launch the app on Streamlit 👈

🎯 Unique Features & USP
Simple: Paste a DNA sequence or upload a FASTA file—no coding or setup required.

Smart Scoring: Each guide is scored with both lab rules and data-inspired heuristics (Hybrid, ML, and Consensus scores).

Off-Target Awareness: Scan your custom background DNA for potential off-targets.

AI Summary: One-click Gemini/OpenAI-powered guide report.

Free & Open Source: Use, share, or modify—no restrictions.

👥 Who is this for? (ICP)
Molecular biology, plant, and biomedical researchers.

Academic labs and teaching settings.

Biotech startups and R&D teams.

Students, beginners, and DIY bio enthusiasts.

🛠️ How to Use
Open the Streamlit app.

Paste your DNA sequence or upload a FASTA file.

Set CRISPR/PAM parameters as needed.

Click Find gRNAs to view candidate guides.

Review scores and download your results.

(Optional) Provide background DNA to scan for off-targets.

(Optional) Simulate edits or request an AI-powered summary.

📊 How are scores calculated?
Score Name	What It Means	Range	How to Use
Hybrid Score	Rule-based: GC%, homopolymers, seed region, off-targets. No ML or training data.	0.0–1.0	>0.85 = Excellent, >0.8 = Good
ML Score	Machine learning-inspired: Based on published gRNA screen patterns. No ML model actually used.	0.0–1.0	>0.7 = Excellent, >0.65 = Good
Consensus Score	Average of Hybrid & ML for balanced ranking. Higher means higher confidence.	0.0–1.0	>0.85 = Excellent, >0.8 = Good

See the app for detailed explanations of each score.

📝 Installation (for local use)
bash
Copy
Edit
git clone https://github.com/yourusername/crispr-guide-designer.git
cd crispr-guide-designer
pip install -r requirements.txt
streamlit run app.py
🔑 License
MIT License – Free for all use.

🙋 FAQ
Is this tool free to use?
Yes, it’s open source and free for all users.

Does this use actual machine learning?
No, scores are rule-based and inspired by ML findings, not trained models.

Is this suitable for crops/animals/bacteria?
Yes, it works for any DNA input, but always validate experimentally.

🤝 Contributing
Pull requests, suggestions, and issue reports are welcome!
See CONTRIBUTING.md for guidelines.

👨‍🔬 Author
Dinesh K – design, code, and documentation

Find your best CRISPR guides—no hassle, no complexity!

“For researchers, students, and anyone wanting reliable CRISPR guide design in seconds.”

