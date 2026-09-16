"""12a_prefill_manual_scores.py - suggested values for the manual criteria used by 12_candidate_prioritization.R

Writes docs/annotation/candidate_manual_scores.csv (one row per candidate gene) with:
  protein_detectability  suggested from Human Protein Atlas IHC reliability
                         (Enhanced 3, Supported 2, Approved 1, Uncertain/none 0)
  druggability           suggested from DGIdb interactions
                         (approved-drug interaction 3, any drug interaction 2, druggable-genome category only 1, none 0)
  literature             left blank: must be scored manually (0-3) by the team
Existing rows are preserved; only empty cells are filled. Sources are recorded per gene.
Usage: python scripts/12a_prefill_manual_scores.py [extra_gene_list.txt]
"""
import csv, json, os, sys, time, urllib.request, urllib.parse

OUT = "docs/annotation/candidate_manual_scores.csv"
genes = [r["gene"] for r in csv.DictReader(open("results/bulk/intersect_genes.csv", encoding="utf-8"))]
if len(sys.argv) > 1:
    genes += [g.strip() for g in open(sys.argv[1], encoding="utf-8") if g.strip()]
genes = list(dict.fromkeys(genes))

cols = ["gene", "protein_detectability", "literature", "druggability",
        "hpa_ihc_reliability", "dgidb_n_interactions", "dgidb_approved_drugs", "dgidb_categories", "notes"]
rows = {}
if os.path.exists(OUT):
    for r in csv.DictReader(open(OUT, encoding="utf-8")):
        rows[r["gene"]] = r
for g in genes:
    rows.setdefault(g, {c: "" for c in cols} | {"gene": g})

def get(url, data=None, headers=None, timeout=60):
    req = urllib.request.Request(url, data=data, headers=headers or {"User-Agent": "invivo-target-validation"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8")

# ---- Human Protein Atlas: IHC reliability ----
hpa_map = {"enhanced": 3, "supported": 2, "approved": 1, "uncertain": 0}
for g in genes:
    r = rows[g]
    if r.get("hpa_ihc_reliability"):
        continue
    try:
        url = ("https://www.proteinatlas.org/api/search_download.php?search=" + urllib.parse.quote(g) +
               "&format=json&columns=g,relih&compress=no")
        res = json.loads(get(url))
        hit = next((x for x in res if x.get("Gene") == g), None)
        rel = (hit or {}).get("Reliability (IH)") or "none"
        r["hpa_ihc_reliability"] = rel
        if r.get("protein_detectability", "") == "":
            r["protein_detectability"] = str(hpa_map.get(str(rel).lower(), 0))
    except Exception as e:
        r["hpa_ihc_reliability"] = r.get("hpa_ihc_reliability") or f"query_failed"
    time.sleep(0.2)

# ---- DGIdb (GraphQL v5) ----
q = """query($names:[String!]!){ genes(names:$names){ nodes{ name
  geneCategories{ name }
  interactions{ drug{ name approved } } } } }"""
try:
    payload = json.dumps({"query": q, "variables": {"names": genes}}).encode("utf-8")
    res = json.loads(get("https://dgidb.org/api/graphql", data=payload,
                         headers={"Content-Type": "application/json", "User-Agent": "invivo-target-validation"}, timeout=120))
    nodes = {n["name"]: n for n in res.get("data", {}).get("genes", {}).get("nodes", [])}
    for g in genes:
        r = rows[g]; n = nodes.get(g)
        inter = (n or {}).get("interactions") or []
        approved = sorted({i["drug"]["name"] for i in inter if i.get("drug") and i["drug"].get("approved")})
        cats = sorted({c["name"] for c in (n or {}).get("geneCategories") or []})
        r["dgidb_n_interactions"] = str(len(inter))
        r["dgidb_approved_drugs"] = ";".join(approved[:10])
        r["dgidb_categories"] = ";".join(cats)
        if r.get("druggability", "") == "":
            r["druggability"] = "3" if approved else ("2" if inter else ("1" if "DRUGGABLE GENOME" in cats else "0"))
except Exception as e:
    print("DGIdb query failed:", e)

with open(OUT, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
    w.writeheader()
    for g in sorted(rows):
        w.writerow({c: rows[g].get(c, "") for c in cols})
print(f"wrote {OUT}: {len(rows)} genes")
