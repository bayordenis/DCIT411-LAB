from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices

matrix_BL = substitution_matrices.load("BLOSUM62")

# Load sequences
from Bio import SeqIO
seqs = list(SeqIO.parse("data/clean_sequences.fasta", "fasta"))

seq1 = str(seqs[0].seq)
seq2 = str(seqs[1].seq)


def identity(aln):
    a, b, score, start, end = aln
    match = sum(x == y for x, y in zip(a, b) if x != "-" and y != "-")
    length = min(len(a.replace("-", "")), len(b.replace("-", "")))
    return (match / length) * 100


print("\n-- Global Alignment --")
ga = pairwise2.align.globalds(seq1, seq2, matrix_BL, -10, -0.5)[0]
print(format_alignment(*ga))
print("Identity:", identity(ga))

print("\n-- Local Alignment --")
la = pairwise2.align.localds(seq1, seq2, matrix_BL, -10, -0.5)[0]
print(format_alignment(*la))
print("Identity:", identity(la))