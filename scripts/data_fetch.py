from Bio import Entrez, SeqIO

Entrez.email = "your.email@example.com"

accessions = [
    "NP_000257.1",  # example human protein
    "XP_011512.1"
]

with open("data/fetched_sequences.fasta", "w") as out:
    for acc in accessions:
        handle = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
        seq = SeqIO.read(handle, "fasta")
        SeqIO.write(seq, out, "fasta")
        handle.close()

print("Fetched sequences saved to data/fetched_sequences.fasta")