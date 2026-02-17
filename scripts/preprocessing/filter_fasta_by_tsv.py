import pandas as pd

# Step 1: Load the TSV file to get the sequence names to keep
tsv_path = "nextclade_filtered.tsv"
tsv_df = pd.read_csv(tsv_path, sep='\t')
sequence_names_to_keep = set(tsv_df['seqName'].astype(str))

# Step 2: Define input and output FASTA paths
fasta_path = "nextclade.aligned.fasta"
output_path = "filtered.nextclade.aligned.fasta"

# Step 3: Read and filter the FASTA
with open(fasta_path, 'r') as infile, open(output_path, 'w') as outfile:
    write = False
    for line in infile:
        if line.startswith('>'):
            seq_id = line[1:].strip().split()[0]
            write = seq_id in sequence_names_to_keep
        if write:
            outfile.write(line)
