import pandas as pd

# Load files
drug_df = pd.read_csv("ASD_Drugs_Merged_with_ATC.csv")
sider_df = pd.read_csv("ASD-iderCleaned.csv")

# Clean column names
drug_df.columns = drug_df.columns.str.strip()
sider_df.columns = sider_df.columns.str.strip()

# Rename sider ATC column
sider_df = sider_df.rename(columns={'atc_code': 'ATC_Codes'})

# Normalize ATC codes (uppercase + remove spaces)
drug_df['ATC_Codes'] = drug_df['ATC_Codes'].fillna('').astype(str).str.upper().str.replace(' ', '')
sider_df['ATC_Codes'] = sider_df['ATC_Codes'].fillna('').astype(str).str.upper().str.replace(' ', '')

# Split multiple ATC codes in drug file
drug_expanded = drug_df.assign(
    ATC_Codes=drug_df['ATC_Codes'].str.split(';')
).explode('ATC_Codes')

drug_expanded['ATC_Codes'] = drug_expanded['ATC_Codes'].str.strip()

# LEFT merge (keeps ALL drugs)
merged = drug_expanded.merge(
    sider_df,
    on='ATC_Codes',
    how='left'
)

# ---- IMPORTANT STEP ----
# Aggregate side effects back to original drug level
side_effect_cols = ['UMLS_concept_id', 'Frequency', 'Side_effect', 'Side_effect_name']

aggregated = merged.groupby(
    ['DrugBank ID', 'Name', 'Description', 'drug_type', 'ATC_Codes'],
    dropna=False
)[side_effect_cols].agg(lambda x: '; '.join(x.dropna().astype(str).unique()))

aggregated = aggregated.reset_index()

# Merge back to original drug file (so no drug is lost)
final_df = drug_df.merge(
    aggregated,
    on=['DrugBank ID', 'Name', 'Description', 'drug_type', 'ATC_Codes'],
    how='left'
)

# Save
final_df.to_csv("ASD_Drugs_with_SideEffects.csv", index=False)

print("Side effects added. No drugs were dropped.")
