import os
import re
import pandas as pd

# Directory with raw sample files
raw_dir = "data/raw"

# Get all sample files (excluding the original tar archive)
all_files = [f for f in os.listdir(raw_dir) if f.endswith(".txt.gz") and ('_IMC' in f or '_IMH' in f)]
all_files.sort()

# Mapping from GSM number (e.g., "GSM5492071") to group label
# Provided by you: first six GSM are Healthy (HC), next five Asymptomatic (AM), last six Symptomatic (SM)
group_map = {
    "GSM5492071": "Healthy",
    "GSM5492072": "Healthy",
    "GSM5492073": "Healthy",
    "GSM5492074": "Healthy",
    "GSM5492075": "Healthy",
    "GSM5492076": "Healthy",
    "GSM5492077": "Asymptomatic",
    "GSM5492078": "Asymptomatic",
    "GSM5492079": "Asymptomatic",
    "GSM5492080": "Asymptomatic",
    "GSM5492081": "Asymptomatic",
    "GSM5492082": "Symptomatic",
    "GSM5492083": "Symptomatic",
    "GSM5492084": "Symptomatic",
    "GSM5492085": "Symptomatic",
    "GSM5492086": "Symptomatic",
    "GSM5492087": "Symptomatic",
}

# Prepare list to store metadata rows
metadata = []
for f in all_files:
    sample_id = f.replace(".txt.gz", "")  # e.g., "GSM5492071_IMC081"
    # Extract GSM number using regex (first part before underscore)
    gsm = sample_id.split("_")[0]         # e.g., "GSM5492071"
    if gsm not in group_map:
        raise ValueError(f"GSM number {gsm} not found in group_map for file {f}")
    group = group_map[gsm]
    metadata.append({"sample_id": sample_id, "group": group})

# Create a DataFrame and save
meta_df = pd.DataFrame(metadata)
# Ensure order matches the column order in the merged count matrix? It's okay, but we can sort by sample_id
meta_df = meta_df.sort_values("sample_id").reset_index(drop=True)

# Save to data/processed/sample_metadata.csv
output_path = "data/processed/sample_metadata.csv"
meta_df.to_csv(output_path, index=False)
print(meta_df)
print(f"Metadata saved to {output_path}")