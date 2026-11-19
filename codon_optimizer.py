"""Educational synonymous codon optimizer.

The host preference values below are illustrative weights to demonstrate the
workflow. They are NOT validated organism-specific codon usage datasets.
"""
from __future__ import annotations

from collections import defaultdict
import math
import re
from typing import Dict, List

class CodonOptimizationError(ValueError):
    """Raised for invalid coding DNA or unsupported input."""

CODON_TABLE: Dict[str, str] = {
    # Phenylalanine / Leucine
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    # Isoleucine / Methionine / Valine
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    # Serine / Proline / Threonine / Alanine
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "AGT":"S","AGC":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    # Tyrosine / Stop / Histidine / Glutamine
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    # Asparagine / Lysine / Aspartate / Glutamate
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    # Cysteine / Stop / Tryptophan
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    # Arginine / Glycine
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGA":"R","AGG":"R",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}
STOP_CODONS = {"TAA", "TAG", "TGA"}

# Explicitly illustrative preference weights. Replace with curated datasets
# (e.g. organism-specific CDS codon frequencies) for scientific use.
HOST_WEIGHTS = {
    "Escherichia coli": {
        "TTT":0.58,"TTC":0.42,"TTA":0.14,"TTG":0.13,"CTT":0.10,"CTC":0.10,"CTA":0.04,"CTG":0.49,
        "ATT":0.49,"ATC":0.39,"ATA":0.12,"ATG":1.0,
        "GTT":0.28,"GTC":0.20,"GTA":0.17,"GTG":0.35,
        "TCT":0.15,"TCC":0.15,"TCA":0.12,"TCG":0.14,"AGT":0.16,"AGC":0.28,
        "CCT":0.18,"CCC":0.13,"CCA":0.20,"CCG":0.49,
        "ACT":0.19,"ACC":0.40,"ACA":0.17,"ACG":0.24,
        "GCT":0.18,"GCC":0.26,"GCA":0.23,"GCG":0.33,
        "TAT":0.57,"TAC":0.43,"CAT":0.57,"CAC":0.43,"CAA":0.34,"CAG":0.66,
        "AAT":0.49,"AAC":0.51,"AAA":0.74,"AAG":0.26,
        "GAT":0.63,"GAC":0.37,"GAA":0.68,"GAG":0.32,
        "TGT":0.46,"TGC":0.54,"TGG":1.0,
        "CGT":0.36,"CGC":0.36,"CGA":0.07,"CGG":0.11,"AGA":0.07,"AGG":0.03,
        "GGT":0.35,"GGC":0.37,"GGA":0.13,"GGG":0.15,
    },
    "Yeast (S. cerevisiae)": {
        "TTT":0.59,"TTC":0.41,"TTA":0.28,"TTG":0.29,"CTT":0.13,"CTC":0.06,"CTA":0.14,"CTG":0.10,
        "ATT":0.46,"ATC":0.26,"ATA":0.28,"ATG":1.0,
        "GTT":0.39,"GTC":0.21,"GTA":0.21,"GTG":0.19,
        "TCT":0.26,"TCC":0.16,"TCA":0.21,"TCG":0.10,"AGT":0.16,"AGC":0.11,
        "CCT":0.31,"CCC":0.15,"CCA":0.42,"CCG":0.12,
        "ACT":0.35,"ACC":0.22,"ACA":0.30,"ACG":0.13,
        "GCT":0.38,"GCC":0.22,"GCA":0.29,"GCG":0.11,
        "TAT":0.56,"TAC":0.44,"CAT":0.64,"CAC":0.36,"CAA":0.69,"CAG":0.31,
        "AAT":0.59,"AAC":0.41,"AAA":0.58,"AAG":0.42,
        "GAT":0.65,"GAC":0.35,"GAA":0.70,"GAG":0.30,
        "TGT":0.63,"TGC":0.37,"TGG":1.0,
        "CGT":0.14,"CGC":0.06,"CGA":0.06,"CGG":0.04,"AGA":0.48,"AGG":0.22,
        "GGT":0.47,"GGC":0.19,"GGA":0.22,"GGG":0.12,
    },
    "Human": {
        "TTT":0.46,"TTC":0.54,"TTA":0.08,"TTG":0.13,"CTT":0.13,"CTC":0.20,"CTA":0.07,"CTG":0.39,
        "ATT":0.36,"ATC":0.48,"ATA":0.16,"ATG":1.0,
        "GTT":0.18,"GTC":0.24,"GTA":0.11,"GTG":0.47,
        "TCT":0.19,"TCC":0.22,"TCA":0.15,"TCG":0.06,"AGT":0.15,"AGC":0.23,
        "CCT":0.28,"CCC":0.33,"CCA":0.27,"CCG":0.12,
        "ACT":0.24,"ACC":0.36,"ACA":0.28,"ACG":0.12,
        "GCT":0.27,"GCC":0.40,"GCA":0.23,"GCG":0.10,
        "TAT":0.43,"TAC":0.57,"CAT":0.41,"CAC":0.59,"CAA":0.27,"CAG":0.73,
        "AAT":0.47,"AAC":0.53,"AAA":0.42,"AAG":0.58,
        "GAT":0.46,"GAC":0.54,"GAA":0.42,"GAG":0.58,
        "TGT":0.45,"TGC":0.55,"TGG":1.0,
        "CGT":0.08,"CGC":0.19,"CGA":0.11,"CGG":0.21,"AGA":0.20,"AGG":0.21,
        "GGT":0.16,"GGC":0.34,"GGA":0.25,"GGG":0.25,
    },
}

def parse_fasta(text: str) -> str:
    """Extract nucleotide sequence from a single FASTA record."""
    lines = [line.strip() for line in text.strip().splitlines() if line.strip()]
    if not lines:
        raise CodonOptimizationError("Please enter a DNA sequence.")
    if not lines[0].startswith(">"):
        raise CodonOptimizationError("FASTA input must begin with a header line starting with '>'.")
    if any(line.startswith(">") for line in lines[1:]):
        raise CodonOptimizationError("Please provide one FASTA record at a time.")
    return "".join(line for line in lines[1:] if not line.startswith(">"))

def validate_sequence(sequence: str) -> str:
    """Normalize and validate a complete coding DNA sequence with terminal stop."""
    seq = re.sub(r"\s+", "", sequence).upper()
    if not seq:
        raise CodonOptimizationError("DNA sequence is empty.")
    invalid = sorted(set(seq) - set("ATGC"))
    if invalid:
        raise CodonOptimizationError("Invalid character(s): " + ", ".join(invalid) + ". Use DNA alphabet A, T, G and C only.")
    if len(seq) % 3:
        raise CodonOptimizationError(f"Sequence length is {len(seq)} nt; a complete coding sequence must be divisible by 3.")
    if not seq.startswith("ATG"):
        raise CodonOptimizationError("For this version, the coding sequence must begin with the standard start codon ATG.")
    if seq[-3:] not in STOP_CODONS:
        raise CodonOptimizationError("For this version, the coding sequence must end with a stop codon (TAA, TAG or TGA).")
    internal_stops = [i + 1 for i in range(0, len(seq)-3, 3) if seq[i:i+3] in STOP_CODONS]
    if internal_stops:
        raise CodonOptimizationError(f"Internal stop codon detected at codon position {internal_stops[0]}.")
    return seq

def translate_dna(sequence: str) -> str:
    """Translate DNA in the provided reading frame, including '*' for stop."""
    return "".join(CODON_TABLE[sequence[i:i+3]] for i in range(0, len(sequence), 3))

def gc_content(sequence: str) -> float:
    """Return GC percentage."""
    return 100.0 * (sequence.count("G") + sequence.count("C")) / len(sequence) if sequence else 0.0

def _synonymous_groups():
    groups = defaultdict(list)
    for codon, amino_acid in CODON_TABLE.items():
        groups[amino_acid].append(codon)
    return groups

SYNONYMOUS = _synonymous_groups()

def optimize_sequence(sequence: str, host: str) -> dict:
    """Replace each sense codon by its highest-weight synonymous codon."""
    seq = validate_sequence(sequence)
    if host not in HOST_WEIGHTS:
        raise CodonOptimizationError(f"Unsupported host: {host}")
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    amino_acids = [CODON_TABLE[c] for c in codons]
    weights = HOST_WEIGHTS[host]
    optimized: List[str] = []
    for codon, aa in zip(codons, amino_acids):
        if aa == "*":
            optimized.append(codon)
            continue
        candidates = SYNONYMOUS[aa]
        # Stable deterministic tie-break by codon string.
        best = max(candidates, key=lambda c: (weights.get(c, 0.0), c == codon, c))
        optimized.append(best)
    optimized_seq = "".join(optimized)
    if translate_dna(seq) != translate_dna(optimized_seq):
        raise CodonOptimizationError("Internal verification failed: amino-acid sequence changed.")
    return {
        "optimized_sequence": optimized_seq,
        "original_codons": codons,
        "optimized_codons": optimized,
        "amino_acids": amino_acids,
        "codons_changed": sum(a != b for a, b in zip(codons, optimized)),
        "codons_total": len(codons),
    }

def calculate_cai(sequence: str, host: str) -> float:
    """Simplified CAI-like geometric mean of illustrative relative weights."""
    seq = validate_sequence(sequence)
    weights = HOST_WEIGHTS[host]
    values = []
    for i in range(0, len(seq)-3, 3):  # exclude terminal stop
        codon = seq[i:i+3]
        aa = CODON_TABLE[codon]
        if aa == "M" or aa == "W":
            continue
        group = SYNONYMOUS[aa]
        maximum = max(weights.get(c, 0.0) for c in group)
        relative = weights.get(codon, 0.0) / maximum if maximum else 0.0
        if relative > 0:
            values.append(relative)
    return math.exp(sum(math.log(v) for v in values) / len(values)) if values else 0.0
