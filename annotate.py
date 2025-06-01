#!/usr/bin/env python3
"""
annotate.py

Combine a scores table with Ensembl CDS + exon data to
add codon + genomic context (±200 bp) for each amino-acid position.

Anchor for the window depends on NeighborFlag:
  - 'a' → first codon base
  - 'z' → third codon base
"""

import json
import time
import urllib3
from pathlib import Path

import pandas as pd
import requests

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

ENSEMBL_REST = "https://rest.ensembl.org"
GRCH38 = "GRCh38"  # coord-system tag for hg38

# --- minimal DNA codon table -------------------------------------------------
CODON2AA = {
    # Phe / Leu
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    # Ser / Tyr / STOP / Cys / Trp
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    # Leu / Pro
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # His / Gln / Arg
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    # Ile / Met / Thr
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # Asn / Lys / Ser / Arg
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    # Val / Ala / Asp / Glu / Gly
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


# -----------------------------------------------------------------------------


def translate_codon(codon: str) -> str:
    """Return one-letter amino-acid for a DNA codon."""
    return CODON2AA.get(codon.upper(), "X")


def load_inputs(scores_csv: str, cds_csv: str, exon_json: str):
    scores = pd.read_csv(scores_csv)
    cds_df = pd.read_csv(cds_csv)
    with open(exon_json) as fh:
        exon_map = json.load(fh)
    cds_dict = cds_df.set_index("gene_name").to_dict(orient="index")
    return scores, cds_dict, exon_map


# -----------------------------------------------------------------------------


def cds_base_to_genomic(exons, strand, cds_base_offset):
    """
    Map an arbitrary CDS base (0-based) to a genomic coordinate.

    Returns (chrom, strand, genomic_pos_1based)
    """
    # order exons 5′→3′ along the transcript
    exons_sorted = sorted(
        exons, key=lambda e: e["start"], reverse=(strand == -1)
    )
    consumed = 0
    for ex in exons_sorted:
        ex_len = ex["end"] - ex["start"] + 1
        if consumed + ex_len > cds_base_offset:
            within = cds_base_offset - consumed
            if strand == 1:
                g_pos = ex["start"] + within
            else:
                g_pos = ex["end"] - within
            return ex["chrom"], strand, g_pos
        consumed += ex_len
    raise ValueError("CDS offset out of bounds")


# -----------------------------------------------------------------------------


def fetch_sequence(chrom, start, end, strand):
    """
    Fetch sequence for chrom:start..end:strand (1-based, inclusive) via Ensembl REST.
    Returns uppercase DNA string (already reverse-complemented when strand = -1).
    """
    region = f"{chrom}:{start}..{end}:{strand}"
    url = (
        f"{ENSEMBL_REST}/sequence/region/homo_sapiens/{region}"
        f"?coord_system_version={GRCH38}"
    )
    r = requests.get(url, headers={"Content-Type": "text/plain"}, verify=False)
    r.raise_for_status()
    return r.text.strip().upper()


# -----------------------------------------------------------------------------


def annotate_scores(scores_df, cds_dict, exon_map, pause=0.2):
    out_rows = []
    for idx, row in scores_df.iterrows():
        gene = row["Gene"]
        data = cds_dict.get(gene)
        if data is None:
            print(f"[WARN] Gene {gene} not in cds_summary.csv – skipping.")
            continue

        gene_id, transcript_id = data["gene_id"], data["transcript_id"]
        cds_seq = data["cds_sequence"].upper()
        exons = exon_map[transcript_id]
        strand = exons[0]["strand"]

        aa_pos = int(row["Position"])
        cds_offset = (aa_pos - 1) * 3  # first base of this codon in CDS space

        # get the codon in CDS space
        codon_seq = cds_seq[cds_offset: cds_offset + 3]
        if len(codon_seq) != 3:
            print(f"[WARN] {gene} position {aa_pos} beyond CDS – skipping.")
            continue

        if translate_codon(codon_seq) != row["Residue"]:
            print(
                f"[WARN] residue mismatch {gene} pos {aa_pos} "
                f"(CSV {row['Residue']} vs codon {translate_codon(codon_seq)}) – skipping."
            )
            continue

        # genomic coordinate of first & third codon bases
        chrom, _, base1 = cds_base_to_genomic(exons, strand, cds_offset)
        _, _, base3 = cds_base_to_genomic(exons, strand, cds_offset + 2)

        # choose anchor according to NeighborFlag
        flag = str(row["NeighborFlag"]).strip().lower()
        anchor = base3 if flag == "z" else base1

        # ±200 bp window around the chosen anchor
        ctx_start = max(anchor - 200, 1)
        ctx_end = anchor + 200
        context_seq = fetch_sequence(chrom, ctx_start, ctx_end, strand)

        out_rows.append(
            {
                **row.to_dict(),
                "gene_id": gene_id,
                "transcript_id": transcript_id,
                "codon_sequence": codon_seq,
                "chrom": chrom,
                "strand": strand,
                "codon_genomic_base1": base1,
                "codon_genomic_base3": base3,
                "context_sequence": context_seq,
            }
        )
        time.sleep(pause)  # polite pause for Ensembl
    return pd.DataFrame(out_rows)


# -----------------------------------------------------------------------------


def main():
    import argparse
    import sys

    p = argparse.ArgumentParser(
        description="Annotate score rows with codon + genomic context"
    )
    p.add_argument(
        "--scores",
        required=True,
        help="CSV with Position,Residue,NeighborFlag,…,Gene columns",
    )
    p.add_argument(
        "--cds",
        required=True,
        help="cds_summary.csv generated earlier",
    )
    p.add_argument(
        "--exons",
        default="transcript_exons.json",
        help="JSON file mapping transcript → exon structure",
    )
    p.add_argument(
        "--out",
        default="annotated_scores.csv",
        help="output CSV file",
    )
    args = p.parse_args()

    scores, cds_dict, exon_map = load_inputs(args.scores, args.cds, args.exons)
    annotated = annotate_scores(scores, cds_dict, exon_map)
    if annotated.empty:
        sys.exit("No rows annotated (see warnings above).")
    annotated.to_csv(args.out, index=False)
    print(f"Wrote {len(annotated)} rows to {args.out}")


# -----------------------------------------------------------------------------


if __name__ == "__main__":
    main()
