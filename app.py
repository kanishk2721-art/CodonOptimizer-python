import streamlit as st
import pandas as pd
from codon_optimizer import (
    optimize_sequence, translate_dna, gc_content, calculate_cai,
    validate_sequence, parse_fasta, CodonOptimizationError
)

st.set_page_config(page_title="CodonOptimizer Python", page_icon="🧬", layout="wide")

st.markdown("""
<style>
.block-container {padding-top: 2.2rem; max-width: 1180px;}
.hero {padding: 1.7rem 2rem; border-radius: 18px; background: linear-gradient(120deg,#14243b,#203c50); color:white; margin-bottom:1.3rem;}
.hero h1 {font-size:2.25rem; margin:0 0 .35rem 0;}
.hero p {color:#d7e5ef; margin:0;}
.small-note {font-size:.88rem;color:#697586;}
</style>
<div class="hero">
<h1>🧬 CodonOptimizer</h1>
<p>Python-powered synonymous codon optimization and sequence analysis</p>
</div>
""", unsafe_allow_html=True)

st.write("Enter a coding DNA sequence or FASTA record, choose a target host, and inspect a synonymous-codon redesign.")
host = st.selectbox("Target host organism", ["Escherichia coli", "Yeast (S. cerevisiae)", "Human"])
input_mode = st.radio("Input format", ["Raw DNA sequence", "FASTA"], horizontal=True)
default_sequence = "ATGGCTGAATTTGGTAACTAA"
raw = st.text_area("DNA sequence", value=default_sequence, height=150,
                   placeholder="Paste A, T, G and C sequence here...")
optimize_button = st.button("Optimize sequence", type="primary", use_container_width=True)

if optimize_button:
    try:
        sequence = parse_fasta(raw) if input_mode == "FASTA" else raw
        sequence = validate_sequence(sequence)
        result = optimize_sequence(sequence, host)
        translated_original = translate_dna(sequence)
        translated_optimized = translate_dna(result["optimized_sequence"])

        st.success("Optimization completed. The amino-acid sequence was preserved.")
        a,b,c,d = st.columns(4)
        a.metric("Original GC", f"{gc_content(sequence):.2f}%")
        b.metric("Optimized GC", f"{gc_content(result['optimized_sequence']):.2f}%")
        c.metric("Codons changed", f"{result['codons_changed']} / {result['codons_total']}")
        d.metric("Simplified CAI", f"{calculate_cai(result['optimized_sequence'], host):.3f}")

        left,right = st.columns(2)
        with left:
            st.subheader("Original DNA")
            st.code(sequence, language=None)
        with right:
            st.subheader("Optimized DNA")
            st.code(result["optimized_sequence"], language=None)

        st.subheader("Translation verification")
        st.write("Original protein")
        st.code(translated_original, language=None)
        st.write("Optimized protein")
        st.code(translated_optimized, language=None)
        if translated_original == translated_optimized:
            st.success("Verified: both sequences translate to the same amino-acid sequence.")
        else:
            st.error("Translation differs. Do not use this result.")

        rows = []
        for i, (old, new) in enumerate(zip(result["original_codons"], result["optimized_codons"]), start=1):
            if old != new:
                rows.append({"Position": i, "Original codon": old, "Optimized codon": new,
                             "Amino acid": result["amino_acids"][i-1]})
        st.subheader("Codon substitutions")
        if rows:
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
        else:
            st.info("No codons required substitution for this host preference table.")

        fasta = f">CodonOptimizer_python|host={host.replace(' ', '_')}\n"
        fasta += "\n".join(result["optimized_sequence"][i:i+60]
                           for i in range(0, len(result["optimized_sequence"]), 60))
        st.download_button("Download optimized FASTA", fasta, file_name="optimized_sequence.fasta",
                           mime="text/plain", use_container_width=True)

        st.caption("The host preference tables in this educational version are illustrative. The displayed CAI is a simplified index, not a validated species-specific prediction.")
    except CodonOptimizationError as exc:
        st.error(str(exc))
    except Exception as exc:
        st.error(f"Unable to process this input: {exc}")

with st.expander("About the method and limitations"):
    st.markdown("""
    - Replaces codons only with synonymous codons, so the encoded protein is retained.
    - Accepts coding DNA containing A, T, G and C; whitespace is removed.
    - Requires a complete sequence of triplets and a terminal stop codon.
    - Host preferences are transparent demonstration weights, not validated codon-usage reference datasets.
    - CAI is simplified; it is not a measured expression value or a guarantee of protein yield.
    - This version does not model mRNA structure, codon-pair effects, ribosome kinetics, regulatory motifs or protein folding.
    - Verify sequences and expression experimentally before research or production use.
    """)

st.divider()
st.caption("CodonOptimizer Python • Streamlit application • Browser-based interface, Python processing")
