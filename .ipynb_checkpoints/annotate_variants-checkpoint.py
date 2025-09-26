#!/usr/bin/env python3
"""
Given a list of rsIDs or genomic coordinates in a CSV, fetch:
- Ensembl VEP (via Ensembl REST) consequence + mapping to gene
- dbSNP presence via Ensembl or NCBI
- ClinVar aggregate status via NCBI E-utilities (esearch/efetch)
- GTEx eQTLs via GTEx REST (best-effort)

Input CSV format: one column "variant" with entries like:
  rs558171988
  19:45410250:G:C
  chr19:45410250:G:C

Usage:
  pip install requests pandas
  python annotate_variants.py --in variants.csv --out annotated.csv

Note: This script uses public REST APIs and is rate-limited; 
it’s suitable for tens to hundreds of variants;
for very large batches use Ensembl VEP offline or bulk ClinVar dumps.
"""
import argparse
import requests
import pandas as pd
import time

ENSEMBL_REST = "https://rest.ensembl.org"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GTEX_BASE = "https://gtexportal.org/rest/v1"

HEADERS = {"Content-Type":"application/json", "Accept":"application/json", "User-Agent":"miRNA-seed-annotator/1.0"}

def query_ensembl_vep(variant):
    # Accept variant as rsID or CHR:POS:REF:ALT (VEP accepts rsIDs or HGVS)
    url = f"{ENSEMBL_REST}/vep/human/hgvs/{variant}"
    resp = requests.get(url, headers=HEADERS, timeout=30)
    if resp.status_code == 200:
        return resp.json()
    # try as rsID endpoint
    url2 = f"{ENSEMBL_REST}/variation/human/{variant}"
    resp2 = requests.get(url2, headers=HEADERS, timeout=30)
    if resp2.status_code == 200:
        return resp2.json()
    return None

def query_ncbi_clinvar_by_rs(rsid):
    # Use esearch db=clinvar to find record for rs#
    params = {"db":"clinvar", "term": rsid, "retmode":"json"}
    r = requests.get(NCBI_EUTILS + "/esearch.fcgi", params=params, timeout=30)
    if r.status_code != 200:
        return None
    js = r.json()
    ids = js.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return None
    # fetch summary (efetch or esummary)
    params2 = {"db":"clinvar", "id": ids[0], "retmode":"xml"}
    r2 = requests.get(NCBI_EUTILS + "/efetch.fcgi", params=params2, timeout=30)
    if r2.status_code == 200:
        return r2.text
    return None

def query_gtex_eqtl_by_variant(rsid_or_variant):
    # GTEx doesn't expose all eQTLs for arbitrary rs easily via the public REST without paging;
    # a practical approach is to use the variant endpoint by variantId (rsid) if present
    # This endpoint returns eQTLs (tissue-level) but may require exact formatting.
    url = f"{GTEX_BASE}/association/singleVariantEqtl?variantId={rsid_or_variant}"
    r = requests.get(url, timeout=30)
    if r.status_code == 200:
        return r.json()
    # fallback: GTEx variant search
    url2 = f"{GTEX_BASE}/variants/search?variantId={rsid_or_variant}"
    r2 = requests.get(url2, timeout=30)
    if r2.status_code == 200:
        return r2.json()
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="input_csv", required=True)
    ap.add_argument("--out", dest="out_csv", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.input_csv)
    if 'variant' not in df.columns:
        raise SystemExit("Input CSV must contain a 'variant' column with rsIDs or CHR:POS:REF:ALT entries")

    rows = []
    for idx, r in df.iterrows():
        variant = str(r['variant']).strip()
        entry = {"variant": variant}
        # Ensembl VEP
        try:
            vep = query_ensembl_vep(variant)
            entry['ensembl_vep'] = vep
        except Exception as e:
            entry['ensembl_vep'] = None
            entry['vep_error'] = str(e)
        # ClinVar via NCBI
        try:
            if variant.startswith("rs"):
                clin = query_ncbi_clinvar_by_rs(variant)
                entry['clinvar_xml'] = clin
            else:
                entry['clinvar_xml'] = None
        except Exception as e:
            entry['clinvar_xml'] = None
            entry['clinvar_error'] = str(e)
        # GTEx
        try:
            gtex = query_gtex_eqtl_by_variant(variant)
            entry['gtex'] = gtex
        except Exception as e:
            entry['gtex'] = None
            entry['gtex_error'] = str(e)

        rows.append(entry)
        time.sleep(0.34)  # be polite ( ~ 3 requests/sec )

    outdf = pd.DataFrame(rows)
    outdf.to_json(args.out_csv, orient='records', indent=2)
    print("Wrote", args.out_csv)

if __name__ == "__main__":
    main()
