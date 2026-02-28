from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices

# Load substitution matrix
matrix = substitution_matrices.load("BLOSUM62")

seq1 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
seq2 = "MKALYIAKQRQISFVKSHFSRQLDEELGLIEVQ"

def calculate_identity(alignment):
    seqA, seqB, score, start, end = alignment
    matches = sum(a == b for a, b in zip(seqA, seqB) if a != "-" and b != "-")
    length = min(len(seqA.replace("-", "")), len(seqB.replace("-", "")))
    return (matches / length) * 100

print("\n=== GLOBAL ALIGNMENT (Needleman–Wunsch) ===")
global_align = pairwise2.align.globalds(seq1, seq2, matrix, -10, -0.5)[0]
print(format_alignment(*global_align))
print(f"Identity: {calculate_identity(global_align):.2f}%")

print("\n=== LOCAL ALIGNMENT (Smith–Waterman) ===")
local_align = pairwise2.align.localds(seq1, seq2, matrix, -10, -0.5)[0]
print(format_alignment(*local_align))
print(f"Identity: {calculate_identity(local_align):.2f}%")