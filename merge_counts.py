import pandas as pd
import gzip
import os

# Directory containing the per-sample count files, relative to repo root
raw_dir = "data/raw"

# Get all files ending with .txt.gz
all_files = [f for f in os.listdir(raw_dir) if f.endswith(".txt.gz")]

# Exclude the original tar archive by checking for '_IMC' or '_IMH' in filename
sample_files = [f for f in all_files if ('_IMC' in f) or ('_IMH' in f)]
# Alternative explicit exclusion: sample_files = [f for f in all_files if f != "GSE181179_raw_counts.txt.gz"]

sample_files.sort()
print(f"Found {len(sample_files)} sample files:")
for f in sample_files:
    print(f"  {f}")

dfs = []
for f in sample_files:
    sample_id = f.replace(".txt.gz", "")
    with gzip.open(os.path.join(raw_dir, f), 'rt') as fh:
        # Use whitespace as delimiter (space or tab)
        df = pd.read_csv(fh, sep='\s+')
    # Print columns for debugging
    print(f"\n{f} columns: {df.columns.tolist()}, shape: {df.shape}")

    # The count column should be the third column (index 2). If not, we will adjust.
    if df.shape[1] < 3:
        print(f"WARNING: {f} has fewer than 3 columns, skipping.")
        continue
    count_col = df.columns[2]
    df = df[['GeneID', count_col]].rename(columns={count_col: sample_id})
    dfs.append(df)

if not dfs:
    raise ValueError("No valid sample files were processed.")

merged = dfs[0]
for df in dfs[1:]:
    merged = merged.merge(df, on='GeneID', how='outer')

output_path = "data/processed/GSE181179_raw_counts_merged.txt"
merged.to_csv(output_path, sep='\t', index=False)
print(f"\nMerged matrix saved to {output_path}")
print(f"Shape: {merged.shape}")
print(merged.head())