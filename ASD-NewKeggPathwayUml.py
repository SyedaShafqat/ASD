import pandas as pd
import spacy
from scispacy.linking import EntityLinker

# ===============================
# FILE PATHS
# ===============================

INPUT_FILE = "ASD-kegg_pathways_cleaned.csv"
OUTPUT_FILE = "ASD-kegg_pathways_with_single_umls.csv"

# ===============================
# LOAD DATA
# ===============================

df = pd.read_csv(INPUT_FILE)

# ===============================
# AUTO-DETECT IMPORTANT COLUMNS
# ===============================

cols_lower = {c.lower(): c for c in df.columns}

def find_col(candidates):
    for c in candidates:
        if c in cols_lower:
            return cols_lower[c]
    return None

pathway_id_col = find_col(["pathway_id", "id", "kegg_id"])
pathway_name_col = find_col(["pathway_name", "name", "title"])
desc_col = find_col(["description", "definition", "summary"])

if pathway_name_col is None:
    raise ValueError(
        f"Pathway name column not found. Columns available: {df.columns.tolist()}"
    )

# Create empty description if missing
if desc_col is None:
    df["description_tmp"] = ""
    desc_col = "description_tmp"

# ===============================
# BUILD TEXT FOR NLP (NAME + DESC)
# ===============================

df["combined_text"] = (
    df[pathway_name_col].fillna("").astype(str).str.strip()
    + ". "
    + df[desc_col].fillna("").astype(str).str.strip()
).str.strip()

# ===============================
# LOAD SCISPACY + UMLS LINKER
# ===============================

print("🔹 Loading SciSpaCy model...")
nlp = spacy.load("en_core_sci_lg")

print("🔹 Adding UMLS linker...")
nlp.add_pipe(
    "scispacy_linker",
    config={
        "resolve_abbreviations": True,
        "name": "umls"
    }
)

# ===============================
# FUNCTION: GET SINGLE BEST UMLS ID
# ===============================

def extract_best_umls(doc):
    best = None

    for ent in doc.ents:
        ent_text = ent.text.lower()

        for cui, score in ent._.kb_ents:
            candidate = {
                "cui": cui,
                "score": score,
                "length": len(ent.text),
                "text": ent.text
            }

            # Prefer longer, more specific phrases
            if best is None:
                best = candidate
            else:
                if candidate["length"] > best["length"]:
                    best = candidate
                elif candidate["length"] == best["length"] and candidate["score"] > best["score"]:
                    best = candidate

    if best:
        return best["cui"], best["score"]
    return None, 0.0


# ===============================
# RUN LINKING
# ===============================

print("🔹 Linking KEGG pathways to ONE UMLS concept (this may take time)...")

best_umls = []
confidence_scores = []

for doc in nlp.pipe(df["combined_text"].tolist(), batch_size=64):
    cui, score = extract_best_umls(doc)
    best_umls.append(cui)
    confidence_scores.append(score)

df["UMLS_ID"] = best_umls
df["UMLS_score"] = confidence_scores
df["UMLS_mapping_method"] = "scispacy_best_score"

# ===============================
# SAVE OUTPUT
# ===============================

df.to_csv(OUTPUT_FILE, index=False)

# ===============================
# SUMMARY
# ===============================

print("✅ Finished successfully")
print("Total pathways:", len(df))
print("Pathways with UMLS ID:", df["UMLS_ID"].notna().sum())
print("\nPreview:")
cols_to_show = [c for c in [pathway_id_col, pathway_name_col, "UMLS_ID", "UMLS_score"] if c]
print(df[cols_to_show].head(10))
