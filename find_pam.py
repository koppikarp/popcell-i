#!/usr/bin/env python3
"""
find_pam.py – locate all NGG / CCN PAMs in coding-oriented sequences.

context_sequence is already 5'→3' in CDS orientation no matter what the
genomic strand column says; we therefore scan that string directly.

windows (offset from anchor in coding coords)
---------------------------------------------
Flag  Motif   min   max
----  ------  ----  ----
 a    NGG     -30   +3
 a    CCN     -3    +30
 z    NGG     -30   +4
 z    CCN     -5    +30
"""

import re, sys, pandas as pd, urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# overlap-aware regexes
PAT_NGG = re.compile(r"(?=([ACGT]GG))")   # look-ahead → overlapping hits
PAT_CCN = re.compile(r"(?=(CC[ACGT]))")

# -------------------------------------------------------------------- #
def genomic_pos_from_index(idx: int,
                           anchor_genomic: int,
                           strand: int,
                           upstream_len: int) -> int:
    """map coding-string index → 1-based hg38 coordinate."""
    if strand == 1:                    # plus-strand gene
        return anchor_genomic - upstream_len + idx
    else:                              # minus-strand gene
        return anchor_genomic + upstream_len - idx

# -------------------------------------------------------------------- #
def find_pams_in_row(row):
    seq    = row["context_sequence"].upper()         # coding orientation
    strand = int(row["strand"])                      # only for coords
    codon  = row["codon_sequence"].upper()
    flag   = str(row["NeighborFlag"]).strip().lower()

    anchor_genomic = (int(row["codon_genomic_base3"])
                      if flag == "z"
                      else int(row["codon_genomic_base1"]))

    # --- pick the codon occ. nearest centre --------------------------------
    hits = [m.start() for m in re.finditer(re.escape(codon), seq)]
    if not hits:
        raise ValueError("codon not found in context_sequence")
    mid          = len(seq) // 2
    idx_codon    = min(hits, key=lambda i: abs(i - mid))
    idx_anchor   = idx_codon if flag != "z" else idx_codon + 2
    upstream_len = idx_anchor

    # --- windows -----------------------------------------------------------
    if flag == "z":
        NGG_min, NGG_max = -30, 4
        CCN_min, CCN_max = -5, 30
    else:                              # 'a'
        NGG_min, NGG_max = -30, 3
        CCN_min, CCN_max = -3, 30

    NGG_hits, CCN_hits = [], []

    # --- scan NGG (overlapping) --------------------------------------------
    for m in PAT_NGG.finditer(seq):
        idx = m.start()                           # index of N
        off = idx - idx_anchor
        if NGG_min <= off <= NGG_max:
            NGG_hits.append(str(genomic_pos_from_index(
                idx, anchor_genomic, strand, upstream_len)))

    # --- scan CCN (overlapping) --------------------------------------------
    for m in PAT_CCN.finditer(seq):
        idx = m.start()                           # first C
        off = idx - idx_anchor
        if CCN_min <= off <= CCN_max:
            CCN_hits.append(str(genomic_pos_from_index(
                idx, anchor_genomic, strand, upstream_len)))

    return ";".join(NGG_hits), ";".join(CCN_hits)

# -------------------------------------------------------------------- #
def main():
    import argparse
    p = argparse.ArgumentParser(description="find NGG / CCN PAMs around codons")
    p.add_argument("--infile",  default="annotated_scores.csv",
                   help="input CSV from annotate.py")
    p.add_argument("--outfile", default="annotated_scores_with_PAMs.csv",
                   help="output CSV with PAM columns added")
    args = p.parse_args()

    df = pd.read_csv(args.infile)
    need = {"NeighborFlag", "context_sequence",
            "codon_genomic_base1", "codon_genomic_base3",
            "codon_sequence", "strand"}
    miss = need - set(df.columns)
    if miss:
        sys.exit(f"missing columns: {', '.join(miss)}")

    pam = df.apply(find_pams_in_row, axis=1, result_type="expand")
    pam.columns = ["NGG_PAMs", "CCN_PAMs"]

    out = pd.concat([df, pam], axis=1)
    out.to_csv(args.outfile, index=False)
    print(f"✓ wrote {len(out)} rows to {args.outfile}")

if __name__ == "__main__":
    main()
