# CRISPR Lab: Advanced gRNA Designer & Impact Analyzer

This app finds CRISPR gRNAs, predicts activity, analyzes off-targets, simulates protein impact—including indels—and connects to AI for edit explanations.

## Features
- gRNA discovery (Cas9/Cas12a, custom PAM/GC)
- Activity scoring
- Detailed off-target mapping
- Protein simulation (robust, codon-correct)
- Bulk indel simulation
- AI edit explanation tab (Gemini/OpenAI ready)
- Downloadable tables and reports

## To run:
```bash
pip install -r requirements.txt
streamlit run app.py
