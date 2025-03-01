import streamlit as st
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from Bio import pairwise2
import pandas as pd
import matplotlib.pyplot as plt

# Title of the Streamlit App
st.title("🔬 Phylogenetic Tree & Sequence Analysis")

# Input box for FASTA sequences
st.subheader("Enter FASTA Sequences")
fasta_input = st.text_area(
    "Paste sequences in FASTA format:",
    """>Seq1
ATGCGTACGTTAG
>Seq2
ATGCGTACGTGAG
>Seq3
ATGCGTTCGTTAG
>Seq4
ATGCGGACGTTAG""",
    height=200,
)


# Function to save input sequences to a FASTA file
def save_fasta(input_text, file_name="sequences.fasta"):
    with open(file_name, "w") as f:
        f.write(input_text.strip())


# Function to analyze sequences (alignment, SNPs, identity)
def analyze_sequences(input_fasta_file, reference_seq):
    alignment_results = []
    alignment = AlignIO.read(input_fasta_file, "fasta")

    for record in alignment:
        seq = str(record.seq)

        # Ensure sequence lengths are the same before alignment
        min_length = min(len(reference_seq), len(seq))
        trimmed_ref = reference_seq[:min_length]
        trimmed_seq = seq[:min_length]

        alignments = pairwise2.align.globalxx(trimmed_ref, trimmed_seq)
        best_alignment = alignments[0]
        similarity = sum(a == b for a, b in zip(best_alignment.seqA, best_alignment.seqB))
        percent_identity = (similarity / len(best_alignment.seqA)) * 100

        # Identify SNPs safely
        snps = [
            (i + 1, ref, alt)
            for i, (ref, alt) in enumerate(zip(trimmed_ref, trimmed_seq))
            if ref != alt
        ]
        snps_str = "; ".join([f"Pos {pos}: {ref} → {alt}" for pos, ref, alt in snps])

        alignment_results.append({
            "Sequence": record.id,
            "Percent Identity": percent_identity,
            "SNPs": snps_str,
            "Best Alignment": f"{best_alignment.seqA}\n{best_alignment.seqB}"
        })

    return pd.DataFrame(alignment_results)


# Function to build phylogenetic tree
def build_phylogenetic_tree(input_fasta_file):
    alignment = AlignIO.read(input_fasta_file, "fasta")
    calculator = DistanceCalculator("identity")
    constructor = DistanceTreeConstructor(calculator, "upgma")
    tree = constructor.build_tree(alignment)

    # Save tree as an image
    fig = plt.figure(figsize=(6, 4), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=ax)
    plt.savefig("tree.png")
    return "tree.png"


# Run analysis when the button is clicked
if st.button("Generate Phylogenetic Tree & Analyze Sequences"):
    if not fasta_input.strip():
        st.error("Please enter sequences in FASTA format!")
    else:
        # Save user input to FASTA file
        fasta_file_path = "sequences.fasta"
        save_fasta(fasta_input, fasta_file_path)
        st.success("FASTA file saved successfully!")

        # Extract reference sequence safely
        fasta_lines = fasta_input.strip().split("\n")
        first_seq = ""
        for line in fasta_lines:
            if not line.startswith(">"):
                first_seq += line.strip()

        if not first_seq:
            st.error("Invalid FASTA format! No sequences found.")
        else:
            # Generate phylogenetic tree
            st.subheader("📌 Phylogenetic Tree")
            try:
                tree_image = build_phylogenetic_tree(fasta_file_path)
                st.image(tree_image, caption="Phylogenetic Tree")
            except Exception as e:
                st.error(f"Error generating tree: {e}")

            # Analyze sequences
            st.subheader("📊 demio.pySequence Analysis")
            try:
                sequence_results = analyze_sequences(fasta_file_path, first_seq)
                st.dataframe(sequence_results)
            except Exception as e:
                st.error(f"Error analyzing sequences: {e}")
