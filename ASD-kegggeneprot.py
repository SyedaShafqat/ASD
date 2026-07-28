import pandas as pd
import re

# ===============================
# FILES
# ===============================

GENE_FILE = "ASD-genes_proteins_newmerged.csv"
PATHWAY_FILE = "ASD-kegg_pathways_with_single_umls.csv"
OUTPUT_FILE = "ASD-genes_proteins_pathways_matched.csv"

# ===============================
# LOAD DATA
# ===============================

genes = pd.read_csv(GENE_FILE)
pathways = pd.read_csv(PATHWAY_FILE)

# ===============================
# DROP UNWANTED COLUMNS
# ===============================

genes = genes.drop(columns=["match_type"], errors="ignore")

pathways = pathways.drop(
    columns=["description_tmp", "UMLS_score", "UMLS_mapping_method"],
    errors="ignore"
)

# ===============================
# HELPER FUNCTIONS
# ===============================

def split_umls(x):
    if pd.isna(x):
        return set()
    return set(u.strip() for u in str(x).split(";") if u.strip())

def normalize(text):
    if pd.isna(text):
        return ""
    text = text.lower()
    text = re.sub(r"[^a-z0-9\s]", " ", text)
    return re.sub(r"\s+", " ", text).strip()

# ===============================
# PREPARE UMLS SETS
# ===============================

genes["gene_umls_set"] = genes.apply(
    lambda row: split_umls(row["UMLS_IDs"]) | split_umls(row["protein_UMLS_IDs"]),
    axis=1
)


pathways["pathway_umls_set"] = pathways["UMLS_ID"].apply(split_umls)

# ===============================
# NORMALIZE TEXT FOR NAME MATCHING
# ===============================

genes["text_blob"] = (
    genes["gene-name"].fillna("") + " " +
    genes["protein_name"].fillna("") + " " +
    genes["protein_function"].fillna("")
).apply(normalize)

pathways["pathway_name_norm"] = pathways["pathway_name"].apply(normalize)

# ===============================
# MATCHING
# ===============================

matches = []

for _, g in genes.iterrows():
    gene_dict = g.drop(["gene_umls_set", "text_blob"]).to_dict()
    gene_matched = False

    for _, p in pathways.iterrows():

        # 1️⃣ UMLS match
        umls_overlap = g["gene_umls_set"] & p["pathway_umls_set"]
        if umls_overlap:
            matches.append({
                **gene_dict,
                "pathway_id": p["pathway_id"],
                "pathway_name": p["pathway_name"],
                "match_basis": "UMLS",
                "matched_umls": ";".join(umls_overlap)
            })
            gene_matched = True
            continue

        # 2️⃣ Name-based fallback
        if p["pathway_name_norm"] and p["pathway_name_norm"] in g["text_blob"]:
            matches.append({
                **gene_dict,
                "pathway_id": p["pathway_id"],
                "pathway_name": p["pathway_name"],
                "match_basis": "NAME",
                "matched_umls": None
            })
            gene_matched = True

    # 🚨 FORCE keep unmatched gene
    if not gene_matched:
        matches.append({
            **gene_dict,
            "pathway_id": None,
            "pathway_name": None,
            "match_basis": None,
            "matched_umls": None
        })

# ===============================
# SAVE RESULT
# ===============================

matched_df = pd.DataFrame(matches)
matched_df.to_csv(OUTPUT_FILE, index=False)

print("✅ Matching complete")
print("Total gene–pathway links:", len(matched_df))
