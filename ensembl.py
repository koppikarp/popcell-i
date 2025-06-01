#!/usr/bin/env python3
"""
ensembl.py  –  fetch CDS and **coding-only** exons for a list of genes.

Input : genes.csv   (must have a column 'gene_name')
Output:
    output_cds_summary.csv       (gene_name,gene_id,transcript_id,cds_sequence)
    output_transcript_exons.json (transcript_id → list[coding exon blocks])
"""

import json, time, pandas as pd, requests

ENSEMBL = "https://rest.ensembl.org"


# -------------------------------------------------------------------- helpers
def get_gene_id(symbol):
    url = f"{ENSEMBL}/lookup/symbol/homo_sapiens/{symbol}?content-type=application/json"
    r = requests.get(url)
    return r.json().get("id") if r.ok else None


def get_canonical_transcript(gene_id):
    url = f"{ENSEMBL}/lookup/id/{gene_id}?expand=1;content-type=application/json"
    r = requests.get(url)
    if not r.ok:
        return None, None
    for t in r.json().get("Transcript", []):
        if t.get("is_canonical"):
            return t["id"], t
    return None, None


def fetch_cds(transcript_id):
    url = f"{ENSEMBL}/sequence/id/{transcript_id}?type=cds"
    r = requests.get(url, headers={"Content-Type": "text/x-fasta"})
    if not r.ok:
        return None
    return "".join(r.text.splitlines()[1:])  # drop FASTA header


# ----------------------------------------------------------------- new core
def coding_blocks_from_exons(exons, translation):
    """
    Intersect each exon (genomic coords) with the CDS span defined by
    translation['start']..['end'].

    Handles both + and – strands transparently.
    """
    cds_low  = min(translation["start"], translation["end"])
    cds_high = max(translation["start"], translation["end"])

    blocks = []
    for ex in exons:
        ov_start = max(ex["start"], cds_low)
        ov_end   = min(ex["end"],   cds_high)
        if ov_start <= ov_end:      # exon contributes coding bases
            blocks.append({
                "chrom":  ex["seq_region_name"],
                "strand": ex["strand"],
                "start":  int(ov_start),
                "end":    int(ov_end)
            })
    return blocks


def fetch_gene_data(symbol):
    gid = get_gene_id(symbol)
    if not gid:
        print(f"[WARN] {symbol}: gene ID not found")
        return None

    tid, tjson = get_canonical_transcript(gid)
    if not tid:
        print(f"[WARN] {symbol}: canonical transcript missing")
        return None

    cds = fetch_cds(tid)
    if not cds:
        print(f"[WARN] {symbol}: CDS fetch failed")
        return None

    cds_exons = coding_blocks_from_exons(tjson["Exon"], tjson["Translation"])

    return {
        "gene_name": symbol,
        "gene_id": gid,
        "transcript_id": tid,
        "cds_sequence": cds,
        "exon_structure": cds_exons,
    }


# -------------------------------------------------------------------- batch
def process_gene_list(csv_file, prefix="output"):
    genes = pd.read_csv(csv_file)["gene_name"]
    rows, exon_map = [], {}
    for g in genes:
        print("Processing", g)
        d = fetch_gene_data(g)
        if d:
            rows.append({k: d[k] for k in
                         ("gene_name", "gene_id", "transcript_id", "cds_sequence")})
            exon_map[d["transcript_id"]] = d["exon_structure"]
        time.sleep(0.2)

    pd.DataFrame(rows).to_csv(f"{prefix}_cds_summary.csv", index=False)
    with open(f"{prefix}_transcript_exons.json", "w") as fh:
        json.dump(exon_map, fh, indent=2)

    print(f"✓ Wrote {prefix}_cds_summary.csv and {prefix}_transcript_exons.json")


# -------------------------------------------------------------------- entry
if __name__ == "__main__":
    process_gene_list("genes.csv")
