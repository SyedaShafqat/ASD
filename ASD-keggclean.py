import pandas as pd

# Load file
df = pd.read_csv("ASD-kegg_pathways.csv")

# Remove everything after ' - ' from pathway_name
df["pathway_name"] = df["pathway_name"].str.split(" - ").str[0].str.strip()

# Save cleaned file
df.to_csv("ASD-kegg_pathways_cleaned.csv", index=False)

print("Done! Cleaned file saved as ASD-kegg_pathways_cleaned.csv")
