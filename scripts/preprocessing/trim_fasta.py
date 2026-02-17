from Bio import SeqIO

input_fasta = "nextclade.aligned.fasta"
output_fasta = "nextclade.aligned.trimmed.fasta"

with open(output_fasta, "w") as out_f:
    for record in SeqIO.parse(input_fasta, "fasta"):
        # Trim 265 bp from start and 230 bp from end (corresponding to the 5' end and 3' end sections)
        record.seq = record.seq[265:-230]
        SeqIO.write(record, out_f, "fasta")