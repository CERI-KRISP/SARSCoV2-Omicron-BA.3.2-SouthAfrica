import pandas as pd

# Input and output paths
input_path = "../nextclade_results/nextclade.tsv"
output_path = "lineage.tsv"

# Load the Nextclade TSV file
df = pd.read_csv(input_path, sep='\t')

# Ensure required columns exist
required_cols = ['seqName', 'Nextclade_pango']
missing = [col for col in required_cols if col not in df.columns]
if missing:
    raise ValueError(f"Missing columns in TSV file: {', '.join(missing)}")

# Extract relevant columns
output_df = df[required_cols]

# Save to TSV
output_df.to_csv(output_path, sep='\t', index=False)

print(f"Lineage data saved to: {output_path}")