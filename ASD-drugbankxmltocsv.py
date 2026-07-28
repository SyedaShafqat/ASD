import xml.etree.ElementTree as ET
import pandas as pd

# ==========================
# File paths
# ==========================
xml_file = "ASD-drugbank.xml"
output_csv = "ASD-drugbankcsv.csv"

# ==========================
# Load XML
# ==========================
tree = ET.parse(xml_file)
root = tree.getroot()

# DrugBank namespace
ns = {"db": "http://www.drugbank.ca"}

rows = []

for drug in root.findall("db:drug", ns):

    def get_text(path):
        element = drug.find(path, ns)
        return element.text.strip() if element is not None and element.text else ""

    # -------------------------
    # Basic Information
    # -------------------------
    drug_id = ""
    for i in drug.findall("db:drugbank-id", ns):
        if i.attrib.get("primary") == "true":
            drug_id = i.text
            break

    name = get_text("db:name")
    description = get_text("db:description")
    drug_type = drug.attrib.get("type", "")

    # -------------------------
    # Groups
    # -------------------------
    groups = "; ".join(
        g.text for g in drug.findall("db:groups/db:group", ns)
        if g.text
    )

    indication = get_text("db:indication")
    pharmacodynamics = get_text("db:pharmacodynamics")
    mechanism = get_text("db:mechanism-of-action")

    # -------------------------
    # Targets
    # -------------------------
    targets = []

    for t in drug.findall("db:targets/db:target", ns):
        pname = t.find("db:polypeptide/db:name", ns)
        gene = t.find("db:polypeptide/db:gene-name", ns)

        if gene is not None and gene.text:
            targets.append(gene.text)

        elif pname is not None and pname.text:
            targets.append(pname.text)

    targets = "; ".join(sorted(set(targets)))

    # -------------------------
    # Enzymes
    # -------------------------
    enzymes = []

    for e in drug.findall("db:enzymes/db:enzyme", ns):
        gene = e.find("db:gene-name", ns)
        if gene is not None and gene.text:
            enzymes.append(gene.text)

    enzymes = "; ".join(sorted(set(enzymes)))

    # -------------------------
    # Transporters
    # -------------------------
    transporters = []

    for t in drug.findall("db:transporters/db:transporter", ns):
        gene = t.find("db:gene-name", ns)
        if gene is not None and gene.text:
            transporters.append(gene.text)

    transporters = "; ".join(sorted(set(transporters)))

    # -------------------------
    # Pathways
    # -------------------------
    pathways = "; ".join(
        p.find("db:name", ns).text
        for p in drug.findall("db:pathways/db:pathway", ns)
        if p.find("db:name", ns) is not None
    )

    # -------------------------
    # ATC Codes
    # -------------------------
    atc = "; ".join(
        a.attrib.get("code", "")
        for a in drug.findall("db:atc-codes/db:atc-code", ns)
    )

    # -------------------------
    # External IDs
    # -------------------------
    ids = []

    for ext in drug.findall("db:external-identifiers/db:external-identifier", ns):

        resource = ext.find("db:resource", ns)
        identifier = ext.find("db:identifier", ns)

        if resource is not None and identifier is not None:
            if resource.text in [
                "SNOMED-CT",
                "MeSH",
                "CAS"
            ]:
                ids.append(f"{resource.text}: {identifier.text}")

    ids = "; ".join(ids)

    rows.append({
        "DrugBank ID": drug_id,
        "Name": name,
        "Description": description,
        "Type": drug_type,
        "Groups": groups,
        "Indication": indication,
        "Pharmacodynamics": pharmacodynamics,
        "Mechanism of Action": mechanism,
        "Targets": targets,
        "Enzymes": enzymes,
        "Transporters": transporters,
        "Pathways": pathways,
        "ATC Codes": atc,
        "SNOMED/MeSH/CAS IDs": ids
    })

# ==========================
# Save CSV
# ==========================
df = pd.DataFrame(rows)
df.to_csv(output_csv, index=False, encoding="utf-8-sig")

print(f"Done! Extracted {len(df)} drugs.")
print(f"Saved as {output_csv}")