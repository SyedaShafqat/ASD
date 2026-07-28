import pandas as pd
import numpy as np
import os

# ===============================
# FILE PATHS
# ===============================
CORE_FILE = "ASD-genes_proteins_pathways_matched.csv"
CLINVAR_FILE = "ASD-variant_summary.txt"
UNIPROT_FILE = "ASD-uniprotkb.tsv"
PPI_FILE = "ASD-protein.links.v12.0.txt"
ALIAS_FILE = "ASD-9606.protein.aliases.v12.0.txt"
EXPR_FILE = "ASD-GSE28521_Cerebellum_expression.csv"
PROBE_UNIPROT_FILE = "ASD-GPL570_probe_to_UNIPROT.csv"

OUTPUT_FILE = "ASD-genes_proteins_pathways_ENRICHED.csv"


# ===============================
# NORMALIZATION
# ===============================
def normalize(df):
    df.columns = (
        df.columns.astype(str)
        .str.strip()
        .str.lower()
        .str.replace(" ", "_")
        .str.replace("-", "_")
        .str.replace("#", "")
    )
    return df


# ===============================
# DATA CLEANING FUNCTIONS
# ===============================
def clean_uniprot_id(ids):
    """Clean UniProt IDs - handle multiple IDs separated by ; or ,"""
    if pd.isna(ids) or ids == "" or str(ids).strip().lower() in ['nan', 'none', '']:
        return ""

    id_str = str(ids).strip()
    # Remove version numbers if present (e.g., P12345.1 -> P12345)
    id_str = id_str.split('.')[0]

    # Take first ID if multiple separated by ; or ,
    if ';' in id_str:
        id_str = id_str.split(';')[0].strip()
    elif ',' in id_str:
        id_str = id_str.split(',')[0].strip()

    return id_str


def check_file_exists(file_path, file_name):
    """Check if file exists and has content"""
    if not os.path.exists(file_path):
        print(f"ERROR: {file_name} not found at {file_path}")
        return False

    file_size = os.path.getsize(file_path)
    if file_size == 0:
        print(f"ERROR: {file_name} is empty")
        return False

    print(f"✓ {file_name}: {file_size:,} bytes")
    return True


# ===============================
# MAIN SCRIPT
# ===============================
print("=" * 60)
print("GENOMIC DATA ENRICHMENT PIPELINE")
print("=" * 60)

# Check all files exist
print("\nChecking input files...")
files_to_check = [
    (CORE_FILE, "Core genes file"),
    (CLINVAR_FILE, "ClinVar file"),
    (UNIPROT_FILE, "UniProt file"),
    (PPI_FILE, "STRING PPI file"),
    (ALIAS_FILE, "STRING aliases file"),
    (EXPR_FILE, "Expression file"),
    (PROBE_UNIPROT_FILE, "Probe-UniProt mapping file")
]

all_files_ok = True
for file_path, file_name in files_to_check:
    if not check_file_exists(file_path, file_name):
        all_files_ok = False

if not all_files_ok:
    print("\n❌ Some input files are missing or empty. Please check file paths.")
    exit(1)

# ===============================
# LOAD FILES WITH ERROR HANDLING
# ===============================
print("\nLoading files...")
try:
    core = normalize(pd.read_csv(CORE_FILE))
    print(f"✓ Core genes: {core.shape[0]:,} genes, {core.shape[1]} columns")

    # Load ClinVar with appropriate chunking if too large
    clinvar_chunks = []
    chunk_size = 1000000
    print("Loading ClinVar (this may take a while)...")
    for chunk in pd.read_csv(CLINVAR_FILE, sep='\t', low_memory=False, chunksize=chunk_size):
        clinvar_chunks.append(chunk)
    clinvar = pd.concat(clinvar_chunks, ignore_index=True)
    clinvar = normalize(clinvar)
    print(f"✓ ClinVar: {clinvar.shape[0]:,} variants, {clinvar.shape[1]} columns")

    uniprot = normalize(pd.read_csv(UNIPROT_FILE, sep='\t'))
    print(f"✓ UniProt: {uniprot.shape[0]:,} entries, {uniprot.shape[1]} columns")

    # Load STRING files
    ppi = normalize(pd.read_csv(PPI_FILE, sep=' '))
    print(f"✓ STRING PPI: {ppi.shape[0]:,} interactions, {ppi.shape[1]} columns")

    aliases = normalize(pd.read_csv(ALIAS_FILE, sep='\t'))
    print(f"✓ STRING Aliases: {aliases.shape[0]:,} mappings, {aliases.shape[1]} columns")

    expr = normalize(pd.read_csv(EXPR_FILE))
    print(f"✓ Expression: {expr.shape[0]:,} probes, {expr.shape[1]} columns")

    probe_uni = normalize(pd.read_csv(PROBE_UNIPROT_FILE))
    print(f"✓ Probe-UniProt: {probe_uni.shape[0]:,} mappings, {probe_uni.shape[1]} columns")

except Exception as e:
    print(f"\n❌ Error loading files: {e}")
    exit(1)

# ===============================
# CLEAN AND PREPARE CORE DATA
# ===============================
print("\n" + "=" * 60)
print("PREPARING CORE DATA")
print("=" * 60)

# Clean UniProt IDs in core
core["uniprot_id_clean"] = core["uniprot_id"].apply(clean_uniprot_id)
unique_uniprot_ids = core["uniprot_id_clean"].nunique()
print(f"Unique clean UniProt IDs in core: {unique_uniprot_ids:,}")
print(f"Sample UniProt IDs: {core['uniprot_id_clean'].head().tolist()}")

# Clean gene symbols in core
core["gene_symbol_clean"] = core["gene_symbol"].astype(str).str.strip()
print(f"Unique gene symbols in core: {core['gene_symbol_clean'].nunique():,}")

# ===============================
# 1. UNIPROT ENRICHMENT
# ===============================
print("\n" + "=" * 60)
print("1. UNIPROT ENRICHMENT")
print("=" * 60)

# Clean UniProt IDs in uniprot data
if "entry" in uniprot.columns:
    uniprot = uniprot.rename(columns={"entry": "uniprot_id"})

uniprot["uniprot_id_clean"] = uniprot["uniprot_id"].apply(clean_uniprot_id)

# Prepare column renaming
column_mapping = {
    "entry_name": "uniprot_entry_name",
    "protein_names": "uniprot_protein_names",
    "gene_names": "uniprot_gene_names",
    "organism": "uniprot_organism",
    "length": "uniprot_protein_length",
    "pathway": "uniprot_pathway",
    "function_[cc]": "uniprot_function",
    "gene_ontology_(cellular_component)": "uniprot_go_cc",
    "involvement_in_disease": "uniprot_disease_association",
    "pharmaceutical_use": "uniprot_pharma_use",
    "pubmed_id": "uniprot_pubmed_ids",
    "doi_id": "uniprot_doi_ids"
}

# Apply renaming for columns that exist
for old_col, new_col in column_mapping.items():
    if old_col in uniprot.columns:
        uniprot = uniprot.rename(columns={old_col: new_col})

# Select columns to keep
keep_cols = ["uniprot_id_clean"] + [new_col for old_col, new_col in column_mapping.items() if
                                    new_col in uniprot.columns]
uniprot = uniprot[keep_cols].drop_duplicates(subset=["uniprot_id_clean"])

print(f"UniProt data after cleaning: {uniprot.shape}")
print(f"Columns available: {uniprot.columns.tolist()}")

# Merge with core
core = core.merge(
    uniprot,
    left_on="uniprot_id_clean",
    right_on="uniprot_id_clean",
    how="left",
    suffixes=("", "_uniprot")
)

print(f"After UniProt enrichment: {core.shape}")
print(f"Genes with UniProt data: {core['uniprot_entry_name'].notna().sum():,}")

# ===============================
# 2. CLINVAR ENRICHMENT
# ===============================
print("\n" + "=" * 60)
print("2. CLINVAR ENRICHMENT")
print("=" * 60)

# Find gene symbol column in ClinVar
gene_col = None
for col in ["genesymbol", "gene_symbol", "gene", "symbol"]:
    if col in clinvar.columns:
        gene_col = col
        break

if gene_col:
    clinvar = clinvar.rename(columns={gene_col: "gene_symbol"})
    print(f"Using '{gene_col}' as gene symbol column")
else:
    print("WARNING: No gene symbol column found in ClinVar")
    # Create empty ClinVar summary
    core["clinvar_variant_count"] = 0
    core["clinvar_clinical_significance"] = ""
    core["clinvar_phenotypes"] = ""
    core["clinvar_review_status"] = ""
    core["clinvar_submitter_count"] = 0

if gene_col:
    # Clean gene symbols
    clinvar["gene_symbol_clean"] = clinvar["gene_symbol"].astype(str).str.strip()

    # Filter to genes in core
    core_genes = set(core["gene_symbol_clean"].unique())
    clinvar_filtered = clinvar[clinvar["gene_symbol_clean"].isin(core_genes)].copy()

    print(f"ClinVar variants for core genes: {clinvar_filtered.shape[0]:,}")
    print(f"Core genes found in ClinVar: {clinvar_filtered['gene_symbol_clean'].nunique():,}")

    if not clinvar_filtered.empty:
        # Create summary statistics
        clinvar_summary = clinvar_filtered.groupby("gene_symbol_clean").agg(
            clinvar_variant_count=("variationid", "nunique"),
            clinvar_clinical_significance=("clinicalsignificance",
                                           lambda x: "; ".join(x.dropna().astype(str).unique()[:3])),
            clinvar_phenotypes=("phenotypelist",
                                lambda x: "; ".join(x.dropna().astype(str).str.split(";").explode().unique()[:5])),
            clinvar_review_status=("reviewstatus",
                                   lambda x: x.mode().iloc[0] if not x.mode().empty else ""),
            clinvar_submitter_count=("numbersubmitters", "max")
        ).reset_index()

        # Merge with core
        core = core.merge(
            clinvar_summary,
            left_on="gene_symbol_clean",
            right_on="gene_symbol_clean",
            how="left"
        )
    else:
        print("WARNING: No ClinVar variants found for core genes")
        core["clinvar_variant_count"] = 0
        core["clinvar_clinical_significance"] = ""
        core["clinvar_phenotypes"] = ""
        core["clinvar_review_status"] = ""
        core["clinvar_submitter_count"] = 0

print(f"After ClinVar enrichment: {core.shape}")
print(f"Genes with ClinVar variants: {(core['clinvar_variant_count'] > 0).sum():,}")

# ===============================
# 3. STRING PPI ENRICHMENT (SIMPLIFIED)
# ===============================
print("\n" + "=" * 60)
print("3. STRING PPI ENRICHMENT")
print("=" * 60)

# Try direct mapping from UniProt IDs
print("Attempting direct UniProt ID matching...")

# Get core UniProt IDs
core_uniprot_ids = set(core["uniprot_id_clean"].unique())
print(f"Core UniProt IDs to match: {len(core_uniprot_ids):,}")

# Check if we can extract UniProt IDs from STRING protein IDs
# STRING IDs are usually like "9606.ENSP..." - we need to map them

# First, let's check what's in the aliases file
print("\nChecking STRING aliases file...")
print(f"Aliases columns: {aliases.columns.tolist()}")
print(f"Aliases source values: {aliases['source'].unique()[:10]}")

# Look for UniProt mappings
if 'source' in aliases.columns:
    uniprot_sources = ['UniProt', 'UniProtKB', 'uniprotkb_ac', 'uniprotkb_id', 'UniProtKB-ID', 'UniProtKB-AC']
    found_sources = [s for s in uniprot_sources if s in aliases['source'].unique()]

    if found_sources:
        print(f"Found UniProt sources: {found_sources}")
        uniprot_aliases = aliases[aliases['source'].isin(found_sources)].copy()

        # Clean UniProt IDs in aliases
        uniprot_aliases['uniprot_id_clean'] = uniprot_aliases['alias'].apply(clean_uniprot_id)

        # Filter to core UniProt IDs
        uniprot_aliases_core = uniprot_aliases[uniprot_aliases['uniprot_id_clean'].isin(core_uniprot_ids)]

        print(f"STRING aliases for core UniProt IDs: {uniprot_aliases_core.shape[0]:,}")
        print(f"Core genes with STRING mappings: {uniprot_aliases_core['uniprot_id_clean'].nunique():,}")

        if not uniprot_aliases_core.empty:
            # Get STRING IDs for core proteins
            core_string_ids = set(uniprot_aliases_core['string_protein_id'].unique())

            # Filter PPI to interactions involving core proteins
            ppi_core = ppi[ppi['protein1'].isin(core_string_ids) | ppi['protein2'].isin(core_string_ids)]

            print(f"PPI interactions involving core proteins: {ppi_core.shape[0]:,}")

            if not ppi_core.empty:
                # Count partners per protein
                # First map STRING IDs back to UniProt
                string_to_uniprot = dict(zip(uniprot_aliases_core['string_protein_id'],
                                             uniprot_aliases_core['uniprot_id_clean']))

                # Count interactions
                interaction_counts = {}
                for _, row in ppi_core.iterrows():
                    prot1 = string_to_uniprot.get(row['protein1'])
                    prot2 = string_to_uniprot.get(row['protein2'])

                    if prot1 and prot1 in core_uniprot_ids:
                        interaction_counts.setdefault(prot1, 0)
                        interaction_counts[prot1] += 1

                    if prot2 and prot2 in core_uniprot_ids:
                        interaction_counts.setdefault(prot2, 0)
                        interaction_counts[prot2] += 1

                # Create DataFrame
                ppi_summary = pd.DataFrame({
                    'uniprot_id_clean': list(interaction_counts.keys()),
                    'num_STRING_partners': list(interaction_counts.values())
                })

                # Merge with core
                core = core.merge(
                    ppi_summary,
                    left_on='uniprot_id_clean',
                    right_on='uniprot_id_clean',
                    how='left'
                )

                print(f"Genes with STRING PPI data: {core['num_STRING_partners'].notna().sum():,}")
                print(f"Average partners per gene: {core['num_STRING_partners'].mean():.1f}")
            else:
                print("WARNING: No PPI interactions found for core proteins")
                core['num_STRING_partners'] = 0
        else:
            print("WARNING: No STRING mappings found for core UniProt IDs")
            core['num_STRING_partners'] = 0
    else:
        print("WARNING: No UniProt sources found in STRING aliases")
        core['num_STRING_partners'] = 0
else:
    print("WARNING: 'source' column not found in aliases file")
    core['num_STRING_partners'] = 0

# ===============================
# 4. EXPRESSION ENRICHMENT
# ===============================
print("\n" + "=" * 60)
print("4. EXPRESSION ENRICHMENT")
print("=" * 60)

# Clean UniProt IDs in probe mapping
probe_uni["uniprot_id_clean"] = probe_uni["uniprot_id"].apply(clean_uniprot_id)

# Remove empty IDs
probe_uni_clean = probe_uni[probe_uni["uniprot_id_clean"] != ""]

print(f"Probe-UniProt mappings after cleaning: {probe_uni_clean.shape[0]:,}")
print(f"Unique UniProt IDs in probe data: {probe_uni_clean['uniprot_id_clean'].nunique():,}")

if not probe_uni_clean.empty:
    # Check overlap with core
    probe_uniprot_ids = set(probe_uni_clean["uniprot_id_clean"].unique())
    common_ids = probe_uniprot_ids.intersection(core_uniprot_ids)

    print(f"UniProt IDs common with core: {len(common_ids):,}")

    if common_ids and "id_ref" in expr.columns and "id" in probe_uni_clean.columns:
        # Merge expression with probe mapping
        expr_merged = expr.merge(
            probe_uni_clean,
            left_on="id_ref",
            right_on="id",
            how="inner"
        )

        print(f"Expression data after merging: {expr_merged.shape}")

        # Get GSM columns
        gsm_cols = [c for c in expr_merged.columns if c.startswith("gsm")]
        print(f"GSM expression columns: {len(gsm_cols)}")

        if gsm_cols:
            # Group by uniprot_id and calculate mean expression
            expr_summary = expr_merged.groupby("uniprot_id_clean")[gsm_cols].mean().reset_index()

            # Calculate summary statistics
            expr_summary["mean_expression"] = expr_summary[gsm_cols].mean(axis=1)
            expr_summary["expression_variance"] = expr_summary[gsm_cols].var(axis=1)
            expr_summary["expression_range"] = expr_summary[gsm_cols].max(axis=1) - expr_summary[gsm_cols].min(axis=1)
            expr_summary["expression_samples"] = len(gsm_cols)

            print(f"Expression summary for {expr_summary.shape[0]:,} proteins")

            # Merge with core
            core = core.merge(
                expr_summary,
                left_on="uniprot_id_clean",
                right_on="uniprot_id_clean",
                how="left",
                suffixes=("", "_expr")
            )

            print(f"Genes with expression data: {core['mean_expression'].notna().sum():,}")
        else:
            print("WARNING: No GSM columns found in expression data")
            core["mean_expression"] = np.nan
            core["expression_variance"] = np.nan
            core["expression_range"] = np.nan
            core["expression_samples"] = 0
    else:
        print("WARNING: No common UniProt IDs or missing columns")
        core["mean_expression"] = np.nan
        core["expression_variance"] = np.nan
        core["expression_range"] = np.nan
        core["expression_samples"] = 0
else:
    print("WARNING: No valid probe-UniProt mappings after cleaning")
    core["mean_expression"] = np.nan
    core["expression_variance"] = np.nan
    core["expression_range"] = np.nan
    core["expression_samples"] = 0

# ===============================
# FINAL CLEANUP AND OUTPUT
# ===============================
print("\n" + "=" * 60)
print("FINAL PROCESSING")
print("=" * 60)

# Fill NaN values
fill_values = {
    "num_STRING_partners": 0,
    "clinvar_variant_count": 0,
    "clinvar_submitter_count": 0,
    "mean_expression": np.nan,
    "expression_variance": np.nan,
    "expression_range": np.nan,
    "expression_samples": 0
}

for col, default in fill_values.items():
    if col in core.columns:
        core[col] = core[col].fillna(default)

# Fill string columns
string_cols = ["clinvar_clinical_significance", "clinvar_phenotypes", "clinvar_review_status"]
for col in string_cols:
    if col in core.columns:
        core[col] = core[col].fillna("")

# Use clean UniProt ID as main ID
if "uniprot_id_clean" in core.columns:
    core["uniprot_id"] = core["uniprot_id_clean"]
    core = core.drop(columns=["uniprot_id_clean"])

if "gene_symbol_clean" in core.columns:
    core = core.drop(columns=["gene_symbol_clean"])

# Reorder columns for readability
priority_cols = [
    'gene_symbol', 'gene_name', 'uniprot_id', 'protein_name',
    'genetic_category', 'gene_score', 'number_of_reports',
    'clinvar_variant_count', 'num_STRING_partners', 'mean_expression'
]

existing_cols = core.columns.tolist()
ordered_cols = []
for col in priority_cols:
    if col in existing_cols:
        ordered_cols.append(col)
        existing_cols.remove(col)

# Add remaining columns sorted alphabetically
ordered_cols.extend(sorted(existing_cols))
core = core[ordered_cols]

# Save to CSV
print(f"\nSaving enriched data to {OUTPUT_FILE}...")
core.to_csv(OUTPUT_FILE, index=False)

# ===============================
# FINAL SUMMARY
# ===============================
print("\n" + "=" * 60)
print("ENRICHMENT SUMMARY")
print("=" * 60)

print(f"Total genes processed: {core.shape[0]:,}")
print(f"Total columns in output: {core.shape[1]}")

# Summary statistics
summary = {
    "UniProt data": core['uniprot_entry_name'].notna().sum(),
    "ClinVar variants": (core['clinvar_variant_count'] > 0).sum(),
    "STRING PPI partners": (core['num_STRING_partners'] > 0).sum(),
    "Expression data": core['mean_expression'].notna().sum()
}

print("\nEnrichment results:")
for item, count in summary.items():
    print(f"  - {item}: {count:,} genes")

# Show top genes by different metrics
print("\nTop 5 genes by ClinVar variants:")
if 'clinvar_variant_count' in core.columns:
    top_clinvar = core.nlargest(5, 'clinvar_variant_count')[['gene_symbol', 'clinvar_variant_count']]
    for _, row in top_clinvar.iterrows():
        print(f"  {row['gene_symbol']}: {int(row['clinvar_variant_count']):,} variants")

print("\nTop 5 genes by STRING partners:")
if 'num_STRING_partners' in core.columns:
    top_ppi = core.nlargest(5, 'num_STRING_partners')[['gene_symbol', 'num_STRING_partners']]
    for _, row in top_ppi.iterrows():
        print(f"  {row['gene_symbol']}: {int(row['num_STRING_partners']):,} partners")

print("\nSample of enriched data:")
sample = core.head(3)[['gene_symbol', 'uniprot_id', 'clinvar_variant_count',
                       'num_STRING_partners', 'mean_expression']]
print(sample.to_string())

print(f"\n✅ Enrichment complete! Output saved to {OUTPUT_FILE}")
print(f"File size: {os.path.getsize(OUTPUT_FILE):,} bytes")