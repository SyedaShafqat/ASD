import pandas as pd
import requests
from lxml import html
import time

# ==========================
# UMLS API CONFIG
# ==========================
API_KEY = "053eeced-a535-4242-b34a-78c1406430a2"
TGT_URL = "https://utslogin.nlm.nih.gov/cas/v1/api-key"
SEARCH_URL = "https://uts-ws.nlm.nih.gov/rest/search/current"

# ==========================
# Get Ticket Granting Ticket
# ==========================
def get_tgt(api_key):
    response = requests.post(TGT_URL, data={'apikey': api_key})
    tree = html.fromstring(response.text)
    tgt = tree.xpath('//form/@action')[0]
    return tgt

# ==========================
# Get Service Ticket
# ==========================
def get_service_ticket(tgt):
    response = requests.post(tgt, data={'service': 'http://umlsks.nlm.nih.gov'})
    return response.text

# ==========================
# Search UMLS CUI
# ==========================
def get_cui(drug_name, tgt):
    try:
        ticket = get_service_ticket(tgt)
        params = {
            'string': drug_name,
            'ticket': ticket,
            'searchType': 'exact'
        }
        r = requests.get(SEARCH_URL, params=params)
        data = r.json()

        results = data.get('result', {}).get('results', [])

        if results and results[0]['ui'] != "NONE":
            return results[0]['ui']
        else:
            return None
    except:
        return None

# ==========================
# Load your file
# ==========================
df = pd.read_csv("Myside_effects_with_drugname_atc.csv", dtype=str)

# Clean names
df["Drug_Name"] = df["Drug_Name"].astype(str).str.strip()

# Remove duplicates to avoid repeated API calls
unique_drugs = df["Drug_Name"].dropna().unique()

# ==========================
# Get TGT once
# ==========================
tgt = get_tgt(API_KEY)

# ==========================
# Fetch CUIs
# ==========================
drug_cui_dict = {}

for drug in unique_drugs:
    print(f"Searching UMLS for: {drug}")
    cui = get_cui(drug, tgt)
    drug_cui_dict[drug] = cui
    time.sleep(0.5)  # prevent API overload

# ==========================
# Map back to dataframe
# ==========================
df["Drug_UMLS_CUI"] = df["Drug_Name"].map(drug_cui_dict)

# ==========================
# Save result
# ==========================
df.to_csv("ASD-SIDER_with_Drug_UMLS.csv", index=False)

print("✅ UMLS CUI added successfully!")
print(df.head())