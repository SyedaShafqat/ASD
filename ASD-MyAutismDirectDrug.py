import pandas as pd
import re

# -----------------------------
# File paths
# -----------------------------
DRUGBANK_FILE = "ASD-drugbankcsv_with_umls.csv"
GENE_FILE = "ASD-genes_proteins_newmerged.csv"
OUTPUT_FILE = "ASD-MyAutismDirectDrug.csv"

# -----------------------------
# Your raw text containing drug names
# -----------------------------
raw_text = """
Fluoxetine
Donepezil
CX516
Citalopram
Buspirone
Valproic acid
Secretin human
Risperidone
Olanzapine
Methylphenidate
Guanfacine
Human immunoglobulin G
Aripiprazole
"""

# -----------------------------
# Load genes/proteins file
# -----------------------------
print("Loading genes/proteins file...")
genes = pd.read_csv(GENE_FILE, dtype=str)

# Clean gene data
for col in ["gene-symbol", "gene-name", "uniprot_id",
            "protein_gene_name", "protein_name"]:
    if col in genes.columns:
        genes[col] = genes[col].str.upper().str.strip()

# -----------------------------
# Create gene/protein mappings (same as indirect file)
# -----------------------------
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
# Function to find matched targets (same as indirect)
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
# Extract unique drug names
# -----------------------------
# Split by newline and clean
drug_names = [d.strip() for d in raw_text.split("\n") if d.strip()]
drug_names = list(set(drug_names))  # unique

print("\nDrugs extracted from text:")
print(drug_names)

# -----------------------------
# Load DrugBank CSV
# -----------------------------
drugbank = pd.read_csv(DRUGBANK_FILE, dtype=str)

# -----------------------------
# Clean for matching (ignore case)
# -----------------------------
drugbank["Name_clean"] = drugbank["Name"].str.upper().str.strip()
drugbank["Targets"] = drugbank["Targets"].str.upper().str.strip()
drug_names_clean = [d.upper() for d in drug_names]

# -----------------------------
# Filter DrugBank for drugs in text
# -----------------------------
filtered_base = drugbank[drugbank["Name_clean"].isin(drug_names_clean)].copy()

print(f"\nFound {len(filtered_base)} matching drugs in DrugBank")

# -----------------------------
# For each matched drug, find its targets (like indirect file)
# -----------------------------
matched_drugs = []

for _, drug_row in filtered_base.iterrows():
    matched_targets = find_matched_targets(drug_row.get("Targets", ""))

    if matched_targets:  # Drug has matching targets in gene list
        # Create a row for each matched target
        for match in matched_targets:
            new_row = drug_row.copy()
            new_row["Relation"] = "DIRECT"
            new_row["Matched_Identifier"] = match['identifier']
            new_row["Matched_Gene_Symbol"] = match['gene_symbol']
            new_row["Matched_Gene_Name"] = match['gene_name']
            matched_drugs.append(new_row)
    else:
        # Drug has no matching targets, still include it but with empty target fields
        new_row = drug_row.copy()
        new_row["Relation"] = "DIRECT"
        new_row["Matched_Identifier"] = ""
        new_row["Matched_Gene_Symbol"] = ""
        new_row["Matched_Gene_Name"] = ""
        matched_drugs.append(new_row)

# Convert to DataFrame
if matched_drugs:
    filtered = pd.DataFrame(matched_drugs)
    print(f"Created DataFrame with {len(filtered)} rows")
else:
    filtered = pd.DataFrame()

# Drop helper column if exists
if "Name_clean" in filtered.columns:
    filtered = filtered.drop(columns=["Name_clean"])

# -----------------------------
# Reorder columns to put new columns near Targets (like indirect)
# -----------------------------
if len(filtered) > 0:
    cols = filtered.columns.tolist()

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

        filtered = filtered[new_cols]

# -----------------------------
# Save to CSV
# -----------------------------
filtered.to_csv(OUTPUT_FILE, index=False)

print(f"\n" + "=" * 60)
print("RESULTS")
print("=" * 60)
print(f"✅ Saved {len(filtered)} rows to {OUTPUT_FILE}")
print(f"✅ Unique drugs: {filtered['DrugBank ID'].nunique() if len(filtered) > 0 else 0}")

print("\n📊 Summary Statistics:")
if len(filtered) > 0:
    drugs_with_targets = filtered[filtered['Matched_Identifier'] != '']['DrugBank ID'].nunique()
    drugs_without_targets = filtered[filtered['Matched_Identifier'] == '']['DrugBank ID'].nunique()

    print(f"Drugs with gene targets: {drugs_with_targets}")
    print(f"Drugs without gene targets: {drugs_without_targets}")
    print(f"Total drug-target pairs: {len(filtered[filtered['Matched_Identifier'] != ''])}")

    if 'Matched_Gene_Symbol' in filtered.columns:
        print("\n🔬 Genes targeted by these drugs:")
        top_genes = filtered[filtered['Matched_Gene_Symbol'] != '']['Matched_Gene_Symbol'].value_counts().head(10)
        for gene, count in top_genes.items():
            print(f"   - {gene}: {count} drugs")

print("\n👁️ Preview:")
preview_cols = ["DrugBank ID", "Name", "Targets", "Matched_Gene_Symbol", "Matched_Gene_Name", "Relation"]
existing_preview = [col for col in preview_cols if col in filtered.columns]
print(filtered[existing_preview].head(20))

# Show drugs without targets
print("\n📋 Drugs in your list that have no gene targets in database:")
drugs_no_targets = filtered[filtered['Matched_Identifier'] == ''][["DrugBank ID", "Name"]].drop_duplicates()
if len(drugs_no_targets) > 0:
    for _, row in drugs_no_targets.iterrows():
        print(f"   - {row['Name']} ({row['DrugBank ID']})")
else:
    print("   All drugs have gene targets!")