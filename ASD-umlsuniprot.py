import requests
import pandas as pd
import time

# =========================
# CONFIG
# =========================
API_KEY = "053eeced-a535-4242-b34a-78c1406430a2"

TGT_URL = "https://utslogin.nlm.nih.gov/cas/v1/api-key"
SEARCH_URL = "https://uts-ws.nlm.nih.gov/rest/search/current"

INPUT_FILE = "ASD-uniprotkb.tsv"
OUTPUT_FILE = "ASD-uniprotkb_with_umls_api.tsv"

# =========================
# STEP 1: Get TGT
# =========================
def get_tgt(api_key):
    r = requests.post(TGT_URL, data={"apikey": api_key})
    if r.status_code != 201:
        raise Exception("Failed to get TGT")
    return r.text.split('action="')[1].split('"')[0]

# =========================
# STEP 2: Get service ticket
# =========================
def get_service_ticket(tgt_url):
    r = requests.post(tgt_url, data={"service": "http://umlsks.nlm.nih.gov"})
    if r.status_code != 200:
        raise Exception("Failed to get service ticket")
    return r.text.strip()

# =========================
# STEP 3: Search UMLS
# =========================
def search_umls(term, tgt_url):
    if not term or str(term).strip() == "":
        return ""

    ticket = get_service_ticket(tgt_url)
    params = {
        "string": term,
        "pageSize": 1,
        "ticket": ticket
    }

    r = requests.get(SEARCH_URL, params=params)
    if r.status_code != 200:
        return ""

    results = r.json().get("result", {}).get("results", [])
    if results:
        return results[0].get("ui", "")

    return ""

# =========================
# STEP 4: Load UniProt
# =========================
print("Reading UniProt file...")
df = pd.read_csv(INPUT_FILE, sep="\t")

# Auto-detect columns
cols_lower = {c.lower(): c for c in df.columns}

entry_col = cols_lower.get("entry")
protein_col = cols_lower.get("protein names")
gene_col = cols_lower.get("gene names", None)

if entry_col is None or protein_col is None:
    raise ValueError("Required UniProt columns not found")

df["UMLS_ID"] = ""

# =========================
# STEP 5: Map proteins
# =========================
tgt_url = get_tgt(API_KEY)
print("Got TGT. Mapping UniProt proteins to UMLS...")

for i, row in df.iterrows():
    entry = row[entry_col]
    protein_name = str(row[protein_col]).split("(")[0].strip()

    try:
        cui = search_umls(protein_name, tgt_url)
        df.at[i, "UMLS_ID"] = cui
        print(f"{i+1}/{len(df)}: {entry} | {protein_name} → {cui}")
        time.sleep(0.2)  # avoid throttling
    except Exception as e:
        print(f"Error for {entry}: {e}")

# =========================
# STEP 6: Save
# =========================
df.to_csv(OUTPUT_FILE, sep="\t", index=False)
print("Saved:", OUTPUT_FILE)
print("Proteins with UMLS:", (df["UMLS_ID"].str.len() > 0).sum())
