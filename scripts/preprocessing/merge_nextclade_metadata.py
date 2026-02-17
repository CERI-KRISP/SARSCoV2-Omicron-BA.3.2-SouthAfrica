import pandas as pd

# Define file paths
nextclade_path = "nextclade_filtered.tsv"
metadata_path = "../../GISAID/Global metadata/metadata_tsv_2025_07_20/metadata_ba3.tsv"
output_path = "nextclade_filtered_gisaid_metadata.tsv"

# Load datasets
nextclade_df = pd.read_csv(nextclade_path, sep="\t", dtype=str)
metadata_df = pd.read_csv(metadata_path, sep="\t", dtype=str)

# Remove spaces from 'Virus name' for matching
metadata_df["Virus name"] = metadata_df["Virus name"].str.replace(" ", "", regex=False)

# Perform left join (nextclade is the base)
merged_df = nextclade_df.merge(metadata_df, how="left", left_on="seqName", right_on="Virus name")

# Reorder columns: metadata columns first, then nextclade columns
metadata_cols = metadata_df.columns.tolist()
nextclade_cols = [col for col in merged_df.columns if col not in metadata_cols]
ordered_cols = metadata_cols + nextclade_cols
merged_df = merged_df[ordered_cols]

# Save to file
merged_df.to_csv(output_path, sep="\t", index=False)

print(f"Merged dataset saved with metadata columns first: {len(merged_df)} records.")