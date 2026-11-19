# CodonOptimizer Python

A browser-based codon optimization demonstrator built with **Python** and **Streamlit**.

## Live app

Live App:
https://codonoptimizer-python-dj7wxecjjfwupdjqgre8xn.streamlit.app/

## Features

- Paste a coding DNA sequence or a single FASTA record
- Select *E. coli*, *S. cerevisiae* (yeast), or human as the target host
- Select preferred synonymous codons using host-specific demonstration weights
- Verify that the translated amino-acid sequence is unchanged
- Calculate GC percentage and a simplified CAI-like score
- Review codon substitutions
- Download the optimized sequence as FASTA

## Run locally

Python 3.9 or newer is recommended.

```bash
python -m venv .venv
# Windows:
.venv\\Scripts\\activate
# macOS/Linux:
source .venv/bin/activate

pip install -r requirements.txt
streamlit run app.py
```

The app opens in your browser, normally at http://localhost:8501.

## Deploy on Streamlit Community Cloud

1. Open https://share.streamlit.io/
2. Sign in with GitHub.
3. Choose repository `kanishk2721-art/CodonOptimizer-python`.
4. Select branch `main` and main file path `app.py`.
5. Deploy.

Streamlit Community Cloud runs the Python application. GitHub Pages is for static sites and is not the correct hosting platform for this Streamlit app.

## Important scientific limitations

The codon preference weights included in `codon_optimizer.py` are illustrative demonstration values, not a curated or validated organism-specific codon-usage dataset. The CAI-like value is simplified and must not be represented as a validated CAI measurement. A high score does not guarantee higher protein expression.

This prototype does not model mRNA secondary structure, codon-pair effects, translational kinetics, regulatory motifs, protein folding, toxicity, or host-specific expression context. For research use, replace the demonstration weights with a documented, versioned reference dataset and validate designs experimentally.

## Input assumptions

- DNA alphabet only: A, T, G, C
- One coding sequence at a time
- Standard start codon ATG
- Complete codons and a terminal stop codon (TAA, TAG, or TGA)
- No internal stop codons

## Project structure

```
CodonOptimizer-python/
├── app.py
├── codon_optimizer.py
├── requirements.txt
└── README.md
```
