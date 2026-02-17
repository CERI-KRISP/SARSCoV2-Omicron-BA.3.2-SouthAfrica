import pandas as pd

# Load Nextclade output file (TSV format)
input_file = "nextclade.tsv"
output_file = "nextclade_filtered.tsv"

# Read the TSV file
df = pd.read_csv(input_file, sep="\t", dtype=str)

# Convert numeric fields properly
df["coverage"] = pd.to_numeric(df["coverage"], errors="coerce") * 100

# Apply filters:
# 1. coverage > 90
# 2. qc.overallStatus != 'bad'
# 3. Nextclade_pango starts with "BA.3"
filtered_df = df[
    (df["coverage"] > 90) &
    (df["qc.overallStatus"] != "bad") &
    (df["Nextclade_pango"].str.startswith("BA.3", na=False))
]

# Write to new TSV
filtered_df.to_csv(output_file, sep="\t", index=False)

print(f"Filtered dataset saved: {len(filtered_df)} sequences retained.")
