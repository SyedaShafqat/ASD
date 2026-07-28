import os

import pandas as pd

# -----------------------------
# File paths
# -----------------------------
GENE_FILE = "ASD-genes_proteins_newmerged.csv"
DRUG_FILE = "ASD-drugbankcsv_with_umls.csv"
OUTPUT_FILE = "ASD-MyAutismIndirectDrug.csv"

# -----------------------------
# Load files
# -----------------------------
genes = pd.read_csv(GENE_FILE, dtype=str)
drugs = pd.read_csv(DRUG_FILE, dtype=str)

print(f"Genes file: {len(genes)} rows")
print(f"Drugs file: {len(drugs)} rows")

# -----------------------------
# Clean text (uppercase + strip)
# -----------------------------
for col in ["gene-symbol", "gene-name", "uniprot_id",
            "protein_gene_name", "protein_name"]:
    if col in genes.columns:
        genes[col] = genes[col].str.upper().str.strip()

if "Targets" in drugs.columns:
    drugs["Targets"] = drugs["Targets"].str.upper().str.strip()
else:
    print("WARNING: 'Targets' column not found in drugs file!")
    # Create empty Targets column if it doesn't exist
    drugs["Targets"] = ""

# -----------------------------
# Create a combined gene/protein set with mappings
# -----------------------------
# Create dictionaries to store identifiers and their preferred names
gene_to_name = {}
gene_to_symbol = {}

for _, row in genes.iterrows():
    # Store gene symbol mapping
    if pd.notna(row.get("gene-symbol")):
        gene_to_name[row["gene-symbol"]] = row.get("gene-name", row["gene-symbol"])
        gene_to_symbol[row["gene-symbol"]] = row["gene-symbol"]

    # Store uniprot mapping
    if pd.notna(row.get("uniprot_id")):
        gene_to_name[row["uniprot_id"]] = row.get("gene-name", row["uniprot_id"])
        gene_to_symbol[row["uniprot_id"]] = row.get("gene-symbol", row["uniprot_id"])

    # Store protein gene name mapping
    if pd.notna(row.get("protein_gene_name")):
        gene_to_name[row["protein_gene_name"]] = row.get("gene-name", row["protein_gene_name"])
        gene_to_symbol[row["protein_gene_name"]] = row.get("gene-symbol", row["protein_gene_name"])

    # Store protein name mapping
    if pd.notna(row.get("protein_name")):
        gene_to_name[row["protein_name"]] = row.get("gene-name", row["protein_name"])
        gene_to_symbol[row["protein_name"]] = row.get("gene-symbol", row["protein_name"])

# Create a set for quick lookup
gene_set = set(gene_to_name.keys())

print(f"Total unique gene/protein identifiers: {len(gene_set)}")


# -----------------------------
# Function to check target match and return matched genes
# -----------------------------
def find_matched_targets(target_string):
    if pd.isna(target_string) or target_string == "":
        return []

    # Split targets if multiple separated by comma/semicolon
    targets = [t.strip() for t in str(target_string).replace(";", ",").split(",")]

    matched = []
    for t in targets:
        t_clean = t.strip()
        if t_clean and t_clean in gene_set:
            matched.append({
                'identifier': t_clean,
                'gene_symbol': gene_to_symbol.get(t_clean, ''),
                'gene_name': gene_to_name.get(t_clean, '')
            })

    return matched


# -----------------------------
# Filter drugs that target ASD genes/proteins
# -----------------------------
matched_drugs = []
drugs_with_matches = 0

for idx, drug_row in drugs.iterrows():
    matched_targets = find_matched_targets(drug_row.get("Targets", ""))

    if matched_targets:  # If there's at least one match
        drugs_with_matches += 1
        # Create a row for each matched target
        for match in matched_targets:
            new_row = drug_row.copy()
            # Explicitly set Relation as INDIRECT
            new_row["Relation"] = "INDIRECT"
            new_row["Matched_Identifier"] = match['identifier']
            new_row["Matched_Gene_Symbol"] = match['gene_symbol']
            new_row["Matched_Gene_Name"] = match['gene_name']
            matched_drugs.append(new_row)

print(f"Drugs with matches: {drugs_with_matches}")

# Convert to DataFrame
if matched_drugs:
    matched_drugs_df = pd.DataFrame(matched_drugs)
    print(f"Created DataFrame with {len(matched_drugs_df)} rows")
else:
    print("No matches found - creating empty DataFrame")
    # Create empty DataFrame with expected columns
    all_cols = list(drugs.columns) + ["Relation", "Matched_Identifier",
                                      "Matched_Gene_Symbol", "Matched_Gene_Name"]
    matched_drugs_df = pd.DataFrame(columns=all_cols)

# -----------------------------
# Verify Relation column exists and has values
# -----------------------------
if "Relation" not in matched_drugs_df.columns:
    matched_drugs_df["Relation"] = "INDIRECT"
else:
    # Check if Relation column has any non-null values
    if matched_drugs_df["Relation"].isna().all():
        matched_drugs_df["Relation"] = "INDIRECT"

    # Fill any null values
    matched_drugs_df["Relation"] = matched_drugs_df["Relation"].fillna("INDIRECT")

print(
    f"Relation column values: {matched_drugs_df['Relation'].unique() if 'Relation' in matched_drugs_df.columns else 'Not found'}")

# -----------------------------
# Reorder columns to put new columns near Targets
# -----------------------------
# Get column list
cols = matched_drugs_df.columns.tolist()

# Find position of Targets column
if "Targets" in cols:
    targets_idx = cols.index("Targets")
    # Insert new columns after Targets
    new_cols_order = (cols[:targets_idx + 1] +
                      ["Matched_Identifier", "Matched_Gene_Symbol", "Matched_Gene_Name", "Relation"] +
                      [col for col in cols[targets_idx + 1:] if
                       col not in ["Relation", "Matched_Identifier", "Matched_Gene_Symbol", "Matched_Gene_Name"]])

    # Remove duplicates while preserving order
    seen = set()
    new_cols = []
    for col in new_cols_order:
        if col not in seen:
            seen.add(col)
            new_cols.append(col)

    matched_drugs_df = matched_drugs_df[new_cols]
else:
    # If Targets not found, put Relation near the front
    if "Relation" in cols:
        # Move Relation to front
        cols.insert(0, cols.pop(cols.index("Relation")))
        matched_drugs_df = matched_drugs_df[cols]

# -----------------------------
# Save
# -----------------------------
matched_drugs_df.to_csv(OUTPUT_FILE, index=False)

print(f"\n" + "=" * 60)
print("RESULTS")
print("=" * 60)
print(f"✅ Found {len(matched_drugs_df)} indirect drug-target pairs")
print(
    f"✅ Unique drugs: {matched_drugs_df['DrugBank ID'].nunique() if len(matched_drugs_df) > 0 and 'DrugBank ID' in matched_drugs_df.columns else 0}")
print(
    f"✅ Relation column: {matched_drugs_df['Relation'].iloc[0] if len(matched_drugs_df) > 0 and 'Relation' in matched_drugs_df.columns else 'No data'}")

print("\n👁️ Preview:")
if len(matched_drugs_df) > 0:
    preview_cols = ["DrugBank ID", "Name", "Targets", "Matched_Identifier",
                    "Matched_Gene_Symbol", "Matched_Gene_Name", "Relation"]
    existing_preview = [col for col in preview_cols if col in matched_drugs_df.columns]
    print(matched_drugs_df[existing_preview].head(10))

    # Show first few rows with focus on Relation
    print("\n🔍 Focusing on Relation column:")
    relation_check = matched_drugs_df[["DrugBank ID", "Name", "Relation"]].head(5)
    print(relation_check.to_string())
else:
    print("No matches found")

# -----------------------------
# Summary statistics
# -----------------------------
print("\n📊 Summary Statistics:")
if len(matched_drugs_df) > 0:
    print(f"Total drug-target pairs: {len(matched_drugs_df)}")
    print(
        f"Unique drugs: {matched_drugs_df['DrugBank ID'].nunique() if 'DrugBank ID' in matched_drugs_df.columns else 'N/A'}")
    print(
        f"Unique matched genes/proteins: {matched_drugs_df['Matched_Identifier'].nunique() if 'Matched_Identifier' in matched_drugs_df.columns else 'N/A'}")

    # Verify Relation values
    if 'Relation' in matched_drugs_df.columns:
        print(f"\n🔍 Relation column value counts:")
        print(matched_drugs_df['Relation'].value_counts(dropna=False))
    else:
        print("WARNING: Relation column missing from final DataFrame!")

    # Top 10 most targeted genes
    if 'Matched_Gene_Symbol' in matched_drugs_df.columns:
        print("\n🔬 Top 10 most targeted genes:")
        top_genes = matched_drugs_df['Matched_Gene_Symbol'].value_counts().head(10)
        for gene, count in top_genes.items():
            if gene and gene != '':  # Only show non-empty values
                print(f"   - {gene}: {count} drugs")

    # Top 10 drugs targeting most genes
    if 'Name' in matched_drugs_df.columns:
        print("\n💊 Top 10 drugs targeting most genes:")
        top_drugs = matched_drugs_df.groupby('Name').size().sort_values(ascending=False).head(10)
        for drug, count in top_drugs.items():
            print(f"   - {drug}: {count} targets")
else:
    print("No matches found in the data")

print(f"\n✅ Output saved to: {OUTPUT_FILE}")

# Final verification
print("\n" + "=" * 60)
print("FINAL VERIFICATION")
print("=" * 60)
if os.path.exists(OUTPUT_FILE):
    verification_df = pd.read_csv(OUTPUT_FILE)
    print(f"File exists: {OUTPUT_FILE}")
    print(f"Rows in saved file: {len(verification_df)}")
    if 'Relation' in verification_df.columns:
        print(f"Relation values in saved file: {verification_df['Relation'].unique()}")
        print(f"Sample of Relation column:")
        print(verification_df[['DrugBank ID', 'Name', 'Relation']].head(3))
    else:
        print("ERROR: Relation column missing from saved file!")
else:
    print(f"ERROR: File {OUTPUT_FILE} was not created!")