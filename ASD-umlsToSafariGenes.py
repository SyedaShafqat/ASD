import requests
import pandas as pd
import time

API_KEY = "053eeced-a535-4242-b34a-78c1406430a2"
TGT_URL = "https://utslogin.nlm.nih.gov/cas/v1/api-key"
SEARCH_URL = "https://uts-ws.nlm.nih.gov/rest/search/current"

# -------------------------
# Step 1: Get Ticket Granting Ticket (TGT)
# -------------------------
def get_tgt(api_key):
    response = requests.post(TGT_URL, data={"apikey": api_key})
    if response.status_code != 201:
        raise Exception(f"Failed to get TGT: {response.text}")
    # TGT URL is in the 'action' of the returned HTML form
    tgt_url = response.text.split('action="')[1].split('"')[0]
    return tgt_url

# -------------------------
# Step 2: Get single-use service ticket
# -------------------------
def get_service_ticket(tgt_url):
    response = requests.post(tgt_url, data={"service": "http://umlsks.nlm.nih.gov"})
    if response.status_code != 200:
        raise Exception(f"Failed to get service ticket: {response.text}")
    return response.text.strip()

# -------------------------
# Step 3: Search UMLS for a term
# -------------------------
def search_umls(term, tgt_url):
    ticket = get_service_ticket(tgt_url)
    params = {"string": term, "pageSize": 1, "ticket": ticket}
    response = requests.get(SEARCH_URL, params=params)
    if response.status_code != 200:
        print(f"Error for term '{term}': {response.status_code}")
        return ""
    results = response.json().get("result", {}).get("results", [])
    if results:
        return results[0]["ui"]
    return ""

# -------------------------
# Step 4: Map gene-symbols
# -------------------------
df = pd.read_csv("safariGenes.csv")
df["UMLS_ID"] = ""

tgt_url = get_tgt(API_KEY)
print("Got TGT. Mapping gene symbols...")

for i, symbol in enumerate(df["gene-symbol"]):
    if pd.isna(symbol) or str(symbol).strip() == "":
        continue
    try:
        cui = search_umls(symbol, tgt_url)
        df.at[i, "UMLS_ID"] = cui
        print(f"{i+1}/{len(df)}: {symbol} -> {cui}")
        time.sleep(0.2)  # avoid throttling
    except Exception as e:
        print(f"Error mapping {symbol}: {e}")

df.to_csv("ASD-safariGenes_with_umls_api.csv", index=False)
print("Done!")
