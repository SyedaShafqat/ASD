import pandas as pd
from collections import defaultdict

# ===============================
# LOAD FILES
# ===============================

genes_df = pd.read_csv("ASD-safariGenes_combined_umls.csv")
proteins_df = pd.read_csv("ASD-uniprotkb_with_umls_api.tsv", sep="\t")

# ===============================
# DROP UNWANTED COLUMNS
# ===============================

genes_df = genes_df.drop(columns=[
    "status", "ensembl-id", "chromosome", "eagle", "syndromic"
], errors="ignore")

proteins_df = proteins_df.drop(columns=[
    "Reviewed", "Length", "Pathway",
    "Pharmaceutical use", "Involvement in disease"
], errors="ignore")

# ===============================
# NORMALIZE COLUMN NAMES
# ===============================

proteins_df = proteins_df.rename(columns={
    "Entry": "uniprot_id",
    "Protein names": "protein_name",
    "Gene Names": "protein_gene_name",
    "Organism": "organism",
    "Function [CC]": "protein_function",
    "UMLS_ID": "UMLS_IDs"
})

# ===============================
# COMBINE GENE UMLS
# ===============================

genes_df["UMLS_IDs"] = (
    genes_df["UMLS_by_symbol"].fillna("") + "," +
    genes_df["UMLS_by_name"].fillna("")
)

genes_df = genes_df.drop(columns=["UMLS_by_symbol", "UMLS_by_name"])

# ===============================
# HELPER FUNCTION
# ===============================

def parse_umls(x):
    if pd.isna(x) or x == "":
        return []
    return [u.strip() for u in str(x).split(",") if u.strip()]

genes_df["UMLS_list"] = genes_df["UMLS_IDs"].apply(parse_umls)
proteins_df["UMLS_list"] = proteins_df["UMLS_IDs"].apply(parse_umls)

# ===============================
# 🔥 BUILD FAST LOOKUP INDEXES
# ===============================

protein_by_umls = defaultdict(list)
protein_by_gene = defaultdict(list)

for _, p in proteins_df.iterrows():

    # index by UMLS
    for umls in p["UMLS_list"]:
        protein_by_umls[umls].append(p)

    # index by gene symbol
    for g in str(p["protein_gene_name"]).lower().split():
        protein_by_gene[g].append(p)

# ===============================
# ⚡ FAST MERGE
# ===============================

merged_rows = []

for _, gene in genes_df.iterrows():
    gene_dict = gene.drop("UMLS_list").to_dict()
    gene_umls = set(gene["UMLS_list"])
    gene_symbol = str(gene["gene-symbol"]).lower()

    matched = False

    # --- Match by UMLS ---
    for umls in gene_umls:
        for p in protein_by_umls.get(umls, []):
            merged_rows.append({
                **gene_dict,
                "uniprot_id": p["uniprot_id"],
                "protein_gene_name": p["protein_gene_name"],
                "protein_name": p["protein_name"],
                "organism": p["organism"],
                "protein_function": p["protein_function"],
                "protein_UMLS_IDs": p["UMLS_IDs"],
                "match_type": "UMLS_ID"
            })
            matched = True

    # --- Fallback: gene symbol ---
    if not matched:
        for p in protein_by_gene.get(gene_symbol, []):
            merged_rows.append({
                **gene_dict,
                "uniprot_id": p["uniprot_id"],
                "protein_gene_name": p["protein_gene_name"],
                "protein_name": p["protein_name"],
                "organism": p["organism"],
                "protein_function": p["protein_function"],
                "protein_UMLS_IDs": p["UMLS_IDs"],
                "match_type": "gene_name"
            })
            matched = True

    # --- Force keep unmatched gene ---
    if not matched:
        merged_rows.append({
            **gene_dict,
            "uniprot_id": None,
            "protein_gene_name": None,
            "protein_name": None,
            "organism": None,
            "protein_function": None,
            "protein_UMLS_IDs": None,
            "match_type": None
        })

# ===============================
# FINAL DATAFRAME
# ===============================

merged_df = pd.DataFrame(merged_rows)

gene_cols = genes_df.drop("UMLS_list", axis=1).columns.tolist()
protein_cols = [
    "uniprot_id",
    "protein_gene_name",
    "protein_name",
    "organism",
    "protein_function",
    "protein_UMLS_IDs",
    "match_type"
]

merged_df = merged_df[gene_cols + protein_cols]

merged_df.to_csv("ASD-genes_proteins_newmerged.csv", index=False)

print("✅ Merge completed successfully")
print("Genes:", genes_df.shape[0])
print("Merged rows:", merged_df.shape[0])
print("Unique genes:", merged_df["gene-symbol"].nunique())
