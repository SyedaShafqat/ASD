import pandas as pd

# ----------------------------
# 1️⃣ Load your ASD drugs file with UMLS IDs
# ----------------------------
asd_drugs = pd.read_csv("ASD_Drugs_Merged_with_UMLS_Name.csv")

# ----------------------------
# 2️⃣ Load DrugBank CSV
# ----------------------------
drugbank = pd.read_csv("ASD-drugbankcsv.csv")

# ----------------------------
# 3️⃣ Keep only needed columns from DrugBank
# ----------------------------
drugbank_subset = drugbank[[
    "DrugBank ID",
    "Indication",
    "Pharmacodynamics",
    "Mechanism of Action",
    "Pathways",
    "Targets"
]]

# ----------------------------
# 4️⃣ Merge on DrugBank ID
# ----------------------------
merged = asd_drugs.merge(
    drugbank_subset,
    left_on="DrugBank ID",
    right_on="DrugBank ID",
    how="left"
)

# ----------------------------
# 5️⃣ Save enriched file
# ----------------------------
merged.to_csv("ASD_drugs_with_targets.csv", index=False)

print(merged.head())
print("Total rows:", len(merged))
