import pandas as pd
from pathlib import Path

# Input and output files
input_file = next(Path("..").glob("gisaid_auspice_input_hcov-19*/*.metadata.tsv"))
df = pd.read_csv(input_file, sep="\t", dtype=str)
#input_file = "../gisaid_auspice_input_hcov-19_2025_07_23_14/1753280730234.metadata.tsv" #BA3.x
#input_file = "../gisaid_auspice_input_hcov-19_2025_07_23_14/1753280862428.metadata.tsv" #BA3
#input_file = "../gisaid_auspice_input_hcov-19_2025_08_23_20/1755981683734.metadata.tsv" #BA3
#input_file = "../gisaid_auspice_input_hcov-19_2025_08_08_09/1754645995182.metadata.tsv" #BA3 Southern Africa
#input_file = "../gisaid_auspice_input_hcov-19_2025_08_18_15/1755529391847.metadata.tsv" #global
#input_file = "../combined_metadata.tsv"
date_file = "dates.tsv"
dates_decimal_file = "dates_decimal.tsv"
location_file = "location_metadata.tsv"
country_file = "country.tsv"
division_file = "division.tsv"

# Read metadata
df = pd.read_csv(input_file, sep="\t", dtype=str)

# Save date metadata (strain | date)
df[["strain", "date"]].to_csv(date_file, sep="\t", index=False)

# Convert collection date to decimal date
def convert_to_decimal(date_str):
    try:
        year, month, day = map(int, date_str.split("-"))
        leap = year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)
        days_in_month = [31, 29 if leap else 28, 31, 30, 31, 30,
                         31, 31, 30, 31, 30, 31]
        doy = sum(days_in_month[:month - 1]) + day
        total_days = 366 if leap else 365
        return round(year + doy / total_days, 8)
    except:
        return "NA"

df["decimal_date"] = df["date"].apply(convert_to_decimal)
df[["strain", "decimal_date"]].to_csv(dates_decimal_file, sep="\t", index=False)

# Split location metadata (strain | continent | country | division | location)
# Assuming "location" column is not reliable (often empty), reconstruct from other fields
location_df = pd.DataFrame({
    "strain": df["strain"],
    "continent": df["region"],
    "country": df["country"],
    "division": df["division"],
    "location": df["location"]  # keep this, even if often blank
})
location_df.to_csv(location_file, sep="\t", header=False, index=False)

# Extract country for separate file
country_df = pd.DataFrame({
    "strain": df["strain"],
    "country": df["country"]
})
country_df.to_csv(country_file, sep="\t", index=False)

# Extract division for separate file
division_df = pd.DataFrame({
    "strain": df["strain"],
    "division": df["division"]
})
division_df.to_csv(division_file, sep="\t", index=False)

print("Date metadata saved to", date_file)
print("Decimal dates saved to", dates_decimal_file)
print("Location metadata saved to", location_file)
print("Country metadata saved to", country_file)
print("Division metadata saved to", division_file)