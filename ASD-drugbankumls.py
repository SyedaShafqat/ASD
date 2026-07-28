import pandas as pd
import spacy
from scispacy.linking import EntityLinker

INPUT_FILE = "ASD-drugbankcsv.csv"
OUTPUT_FILE = "ASD-drugbankcsv_with_umls.csv"

# ----------------------------
# Load DrugBank CSV
# ----------------------------
df = pd.read_csv(INPUT_FILE)

if "Name" not in df.columns:
    raise ValueError(f"Column 'Name' not found. Columns: {df.columns.tolist()}")

# ----------------------------
# Load SciSpaCy + UMLS
# ----------------------------
print("Loading SciSpaCy model...")
nlp = spacy.load("en_core_sci_sm")

print("Adding UMLS linker...")
nlp.add_pipe("scispacy_linker", config={"name": "umls"})
linker = nlp.get_pipe("scispacy_linker")
kb = linker.kb

# ----------------------------
# Extract SINGLE BEST UMLS ID
# ----------------------------
def extract_single_drug_umls(drug_name: str):
    if not isinstance(drug_name, str) or not drug_name.strip():
        return None

    doc = nlp(drug_name)
    best = None

    for ent in doc.ents:
        for cui, score in ent._.kb_ents:
            if cui not in kb.cui_to_entity:
                continue

            entity = kb.cui_to_entity[cui]

            # Keep ONLY drug / chemical concepts
            semantic_types = entity.types
            if not any(t.startswith(("T109", "T121", "T195", "T200")) for t in semantic_types):
                continue

            candidate = {
                "cui": cui,
                "score": score,
                "name": entity.canonical_name
            }

            if best is None or candidate["score"] > best["score"]:
                best = candidate

    if best:
        return best["cui"]

    return None

# ----------------------------
# RUN EXTRACTION
# ----------------------------
print("Extracting ONE UMLS ID per drug...")
df["UMLS_ID"] = df["Name"].apply(extract_single_drug_umls)

# ----------------------------
# SAVE
# ----------------------------
df.to_csv(OUTPUT_FILE, index=False)

print("Saved:", OUTPUT_FILE)
print("Drugs with UMLS ID:", df["UMLS_ID"].notna().sum())
print(df[["Name", "UMLS_ID"]].head(10))
