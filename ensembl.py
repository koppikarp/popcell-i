import pandas as pd
import requests
import time
import json

def get_ensembl_gene_id(gene_name):
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_name}?content-type=application/json"
    response = requests.get(url)
    if response.ok:
        return response.json().get("id")
    else:
        print(f"[WARN] Could not get gene ID for {gene_name}")
        return None

def get_canonical_transcript_info(gene_id):
    url = f"https://rest.ensembl.org/lookup/id/{gene_id}?expand=1;content-type=application/json"
    response = requests.get(url)
    if not response.ok:
        return None, None
    transcripts = response.json().get("Transcript", [])
    for t in transcripts:
        if t.get("is_canonical"):
            return t["id"], t  # return transcript ID and full transcript metadata
    return None, None

def get_cds_sequence(transcript_id):
    url = f"https://rest.ensembl.org/sequence/id/{transcript_id}?type=cds"
    headers = {"Content-Type": "text/x-fasta"}
    response = requests.get(url, headers=headers)
    if not response.ok:
        return None
    fasta = response.text.split('\n')
    return ''.join(fasta[1:])  # remove header, join sequence lines

def fetch_gene_data(gene_name):
    gene_id = get_ensembl_gene_id(gene_name)
    if not gene_id:
        return None
    transcript_id, transcript_data = get_canonical_transcript_info(gene_id)
    if not transcript_id or not transcript_data:
        return None
    cds_seq = get_cds_sequence(transcript_id)
    if not cds_seq:
        return None
    # extract exon structure
    exons = transcript_data.get("Exon", [])
    exon_list = [
        {
            "start": e["start"],
            "end": e["end"],
            "strand": transcript_data["strand"],
            "chrom": transcript_data["seq_region_name"]
        }
        for e in exons
    ]
    return {
        "gene_name": gene_name,
        "gene_id": gene_id,
        "transcript_id": transcript_id,
        "cds_sequence": cds_seq,
        "exon_structure": exon_list
    }

def process_gene_list(input_csv, output_prefix="output"):
    df = pd.read_csv(input_csv)
    results = []
    exon_map = {}
    for gene in df['gene_name']:
        print(f"Processing: {gene}")
        data = fetch_gene_data(gene)
        if data:
            results.append({
                "gene_name": data["gene_name"],
                "gene_id": data["gene_id"],
                "transcript_id": data["transcript_id"],
                "cds_sequence": data["cds_sequence"]
            })
            exon_map[data["transcript_id"]] = data["exon_structure"]
        time.sleep(0.2)  # Respect Ensembl rate limits

    # Write CDS info to CSV
    pd.DataFrame(results).to_csv(f"{output_prefix}_cds_summary.csv", index=False)

    # Write exon structure to JSON
    with open(f"{output_prefix}_transcript_exons.json", "w") as f:
        json.dump(exon_map, f, indent=2)

# Example usage
if __name__ == "__main__":
    process_gene_list("genes.csv")  # Input CSV with column "gene_name"
