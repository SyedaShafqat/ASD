import pandas as pd

# -----------------------------
# File paths
# -----------------------------
DIRECT_FILE = "ASD-MyAutismDirectDrug.csv"
INDIRECT_FILE = "ASD-MyAutismIndirectDrug.csv"
OUTPUT_FILE = "ASD-MyAutismFullDrug.csv"

# -----------------------------
# Load files
# -----------------------------
direct_df = pd.read_csv(DIRECT_FILE, dtype=str)
indirect_df = pd.read_csv(INDIRECT_FILE, dtype=str)

# -----------------------------
# Combine both
# -----------------------------
combined = pd.concat([direct_df, indirect_df], ignore_index=True)

print(f"Total before removing duplicates: {len(combined)}")

# -----------------------------
# Remove duplicates by DrugBank ID
# -----------------------------
combined = combined.drop_duplicates(subset=["DrugBank ID"], keep="first")

print(f"Total after removing duplicates: {len(combined)}")

# -----------------------------
# Save final file
# -----------------------------
combined.to_csv(OUTPUT_FILE, index=False)

print(f"✅ Final merged file saved as {OUTPUT_FILE}")
print("Preview:")
print(combined[["DrugBank ID", "Name", "Relation"]].head(20))