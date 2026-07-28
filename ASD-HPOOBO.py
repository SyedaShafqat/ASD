import pandas as pd

# -----------------------------
# File paths
# -----------------------------
OBO_FILE = "ASD-hp.obo"  # your HPO OBO file
OUTPUT_FILE = "ASD-HPO_CLEAN.csv"

# -----------------------------
# Containers for parsed data
# -----------------------------
terms = []

# -----------------------------
# Read OBO file
# -----------------------------
with open(OBO_FILE, "r", encoding="utf-8") as f:
    current_term = {}
    for line in f:
        line = line.strip()
        if line == "[Term]":
            if current_term:  # save previous term
                terms.append(current_term)
            current_term = {
                "hpo_id": None,
                "hpo_name": None,
                "definition": None,
                "synonyms": [],
                "is_a": [],
                "alt_id": [],
                "xref": []
            }
        elif line.startswith("id: HP:"):
            current_term["hpo_id"] = line.split("id:")[1].strip()
        elif line.startswith("name:"):
            current_term["hpo_name"] = line.split("name:")[1].strip()
        elif line.startswith("def:"):
            # remove brackets for reference
            def_text = line.split("def:")[1].split("[")[0].strip().strip('"')
            current_term["definition"] = def_text
        elif line.startswith("synonym:"):
            syn_text = line.split('"')[1].strip()
            current_term["synonyms"].append(syn_text)
        elif line.startswith("is_a:"):
            parent = line.split("is_a:")[1].split("!")[0].strip()
            current_term["is_a"].append(parent)
        elif line.startswith("alt_id:"):
            alt = line.split("alt_id:")[1].strip()
            current_term["alt_id"].append(alt)
        elif line.startswith("xref:"):
            xref = line.split("xref:")[1].strip()
            current_term["xref"].append(xref)
    # Append last term
    if current_term:
        terms.append(current_term)

# -----------------------------
# Convert to DataFrame
# -----------------------------
df = pd.DataFrame(terms)

# Join list columns with '|'
for col in ["synonyms", "is_a", "alt_id", "xref"]:
    df[col] = df[col].apply(lambda x: "|".join(x) if x else None)

# -----------------------------
# Save CSV
# -----------------------------
df.to_csv(OUTPUT_FILE, index=False)
print(f"✅ HPO CSV saved: {OUTPUT_FILE}")
print(df.head())
