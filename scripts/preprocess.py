from Bio import SeqIO

input_file = "data/fetched_sequences.fasta"
output_file = "data/clean_sequences.fasta"

sequences = SeqIO.parse(input_file, "fasta")

with open(output_file, "w") as out:
    for record in sequences:
        seq = str(record.seq).replace("-", "").upper()
        record.seq = seq
        record.description = ""
        out.write(record.format("fasta"))

print("Preprocessing complete!")