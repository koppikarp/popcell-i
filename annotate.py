#!/usr/bin/env python3
"""
annotate.py

Combine a scores table with Ensembl CDS + exon data to
add codon + genomic context (±100 bp) for each amino-acid position.
"""

import json
import requests
import time
import pandas as pd
from pathlib import Path
from textwrap import wrap
import urllib3

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

ENSEMBL_REST = "https://rest.ensembl.org"
GRCH38       = "GRCh38"        # coord-system tag for hg38

# --- a minimal codon table (DNA) ---
CODON2AA = {
    # Phe / Leu
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    # Ser / Tyr / STOP / Cys / Trp
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    # Leu / Pro
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    # His / Gln / Arg
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    # Ile / Met / Thr
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    # Asn / Lys / Ser / Arg
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    # Val / Ala / Asp / Glu / Gly
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def translate_codon(codon: str) -> str:
    """Return one-letter amino acid for an uppercase DNA codon."""
    return CODON2AA.get(codon.upper(), 'X')

def load_inputs(scores_csv: str, cds_csv: str, exon_json: str):
    scores = pd.read_csv(scores_csv)
    cds_df = pd.read_csv(cds_csv)
    with open(exon_json) as fh:
        exon_map = json.load(fh)
    cds_dict = cds_df.set_index("gene_name").to_dict(orient="index")
    return scores, cds_dict, exon_map

def offset_to_genomic(transcript_exons, strand, cds_offset):
    """
    Given 0-based CDS offset, return chrom, strand (+1/-1),
    and genomic start, end (inclusive) of the first codon base.
    """
    # Sort exons 5'→3' (transcript direction)
    exons = sorted(
        transcript_exons,
        key=lambda e: e["start"],
        reverse=(strand == -1)
    )

    consumed = 0
    for ex in exons:
        ex_len = ex["end"] - ex["start"] + 1
        if consumed + ex_len > cds_offset:
            within = cds_offset - consumed
            if strand == 1:
                g_start = ex["start"] + within
            else:  # negative strand: coords decrease 5'→3'
                g_start = ex["end"] - within
            g_end = g_start + 2 if strand == 1 else g_start - 2
            if strand == -1 and g_end > g_start:
                g_start, g_end = g_end, g_start
            return ex["chrom"], strand, g_start, g_end
        consumed += ex_len
    raise ValueError("CDS offset out of bounds")

def fetch_sequence(chrom, start, end, strand):
    """
    Fetch sequence for chrom:start..end:strand (1-based, inclusive) using Ensembl REST.
    Returns uppercase DNA string.
    """
    # make sure start <= end for the URL
    region = f"{chrom}:{start}..{end}:{strand}"
    url = f"{ENSEMBL_REST}/sequence/region/homo_sapiens/{region}?coord_system_version={GRCH38}"
    r = requests.get(url, headers={"Content-Type": "text/plain"}, verify=False)
    r.raise_for_status()
    return r.text.strip()

def annotate_scores(scores_df, cds_dict, exon_map, pause=0.2):
    out_rows = []
    for idx, row in scores_df.iterrows():
        gene = row["Gene"]
        data = cds_dict.get(gene)
        if data is None:
            print(f"[WARN] Gene {gene} not found in cds_summary.csv – skipping.")
            continue

        gene_id       = data["gene_id"]
        transcript_id = data["transcript_id"]
        cds_seq       = data["cds_sequence"].upper()
        exons         = exon_map[transcript_id]
        strand        = exons[0]["strand"]   # all exons share strand

        aa_pos = int(row["Position"])
        cds_offset = (aa_pos - 1) * 3
        codon_seq = cds_seq[cds_offset : cds_offset + 3]

        if len(codon_seq) != 3:
            print(f"[WARN] {gene} position {aa_pos} beyond CDS – skipping.")
            continue

        aa_expected = row["Residue"]
        aa_observed = translate_codon(codon_seq)
        if aa_observed != aa_expected:
            print(f"[WARN] residue mismatch {gene} pos {aa_pos} "
                  f"(csv {aa_expected} vs codon {aa_observed}) – skipping.")
            continue

        chrom, strand_sign, g_start, g_end = offset_to_genomic(exons, strand, cds_offset)

        # ±100 bp context
        ctx_start = g_start - 100 if strand_sign == 1 else g_end - 100
        ctx_end   = g_end + 100   if strand_sign == 1 else g_start + 100
        if ctx_start < 1:
            ctx_start = 1
        # Fetch sequence (max 201 bp) - note: for negative strand Ensembl
        # returns the reverse-complement automatically when strand = -1.
        context_seq = fetch_sequence(chrom, ctx_start, ctx_end, strand_sign)
        # sanity: might be shorter at chromosome ends
        context_seq = context_seq.upper()

        out_row = row.to_dict()
        out_row.update({
            "gene_id": gene_id,
            "transcript_id": transcript_id,
            "codon_sequence": codon_seq,
            "chrom": chrom,
            "strand": strand_sign,
            "codon_genomic_start": g_start,
            "codon_genomic_end": g_end,
            "context_sequence": context_seq,
        })
        out_rows.append(out_row)
        time.sleep(pause)       # keep REST API happy

    return pd.DataFrame(out_rows)

def main():
    import argparse, sys
    p = argparse.ArgumentParser(
        description="Annotate score rows with codon + genomic context")
    p.add_argument("--scores", required=True,
                   help="CSV with Position,Residue,…,Gene columns")
    p.add_argument("--cds", required=True,
                   help="cds_summary.csv generated earlier")
    p.add_argument("--exons", default="transcript_exons.json",
                   help="JSON file with transcript→exon structure")
    p.add_argument("--out", default="annotated_scores.csv",
                   help="output CSV file")
    args = p.parse_args()

    scores, cds_dict, exon_map = load_inputs(args.scores, args.cds, args.exons)
    annotated = annotate_scores(scores, cds_dict, exon_map)
    if annotated.empty:
        sys.exit("No rows annotated (see warnings above).")
    annotated.to_csv(args.out, index=False)
    print(f"Wrote {len(annotated)} rows to {args.out}")

if __name__ == "__main__":
    main()
