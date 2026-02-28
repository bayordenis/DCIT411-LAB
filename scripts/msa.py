from Bio import AlignIO
import os
import subprocess

input_fasta = "data/sequences.fasta"
output_alignment = "results/alignments/aligned.fasta"

os.makedirs("results/alignments", exist_ok=True)

# Run MUSCLE alignment via subprocess
subprocess.run(["muscle", "-align", input_fasta, "-output", output_alignment], check=True)

alignment = AlignIO.read(output_alignment, "fasta")

print("\n=== MULTIPLE SEQUENCE ALIGNMENT ===")
print(alignment)

# Consensus sequence using modern Biopython API
new_alignment = alignment.alignment
consensus_array = []
for pos in range(new_alignment.shape[1]):
    col = new_alignment[:, pos]
    # Get most common residue
    from collections import Counter
    residue_counts = Counter(col)
    if residue_counts:
        consensus = residue_counts.most_common(1)[0][0]
        if consensus != '-':  # Skip gaps
            consensus_array.append(consensus)

consensus = ''.join(consensus_array)

print("\nConsensus Sequence:")
print(consensus)