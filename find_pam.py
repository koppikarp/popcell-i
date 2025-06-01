#!/usr/bin/env python3
"""
find_pam.py  –  search CRISPR-Cas9 NGG / CCN PAMs around each codon.

Input   : annotated_scores.csv  (from annotate.py)
Output  : annotated_scores_with_PAMs.csv  (adds NGG_PAMs, CCN_PAMs columns)

Anchor base rules
-----------------
NeighborFlag == 'a'  →  first codon base  (codon_genomic_base1)
NeighborFlag == 'z'  →  third codon base  (codon_genomic_base3)

Windows (bp offset from anchor, genomic coords)
-----------------------------------------------
Flag  Motif   min   max
----  ------  ----  ----
 a    NGG     -30   +3
 a    CCN     -3    +30
 z    NGG     -30   +4
 z    CCN     -5    +30

context_sequence is already in coding orientation on both strands.
"""

import re
import sys
import pandas as pd
import urllib3

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


# --------------------------------------------------------------------- #
# coordinate conversion helper
# --------------------------------------------------------------------- #
def genomic_pos_from_index(idx, anchor_genomic, strand, upstream_len):
    """
    Convert a 0-based index inside context_sequence to an absolute hg38
    coordinate (1-based).

    idx           : index of base in the string
    anchor_genomic: genomic coordinate of the anchor base
    strand        : +1 or –1
    upstream_len  : bases to the *left* of anchor inside the string
    """
    if strand == 1:
        ctx_start = anchor_genomic - upstream_len
        return ctx_start + idx
    else:  # strand == -1
        ctx_end = anchor_genomic + upstream_len
        return ctx_end - idx


# --------------------------------------------------------------------- #
# per-row PAM finder
# --------------------------------------------------------------------- #
def find_pams_in_row(row):
    seq    = row["context_sequence"].upper()
    strand = int(row["strand"])
    codon  = row["codon_sequence"].upper()
    flag   = str(row["NeighborFlag"]).strip().lower()

    anchor_genomic = (int(row["codon_genomic_base3"])
                      if flag == "z"
                      else int(row["codon_genomic_base1"]))

    # ------------------------------------------------------------------
    # 1. locate the *correct* occurrence of the codon
    # ------------------------------------------------------------------
    best_idx_anchor = None
    for m in re.finditer(re.escape(codon), seq):
        idx_codon = m.start()                      # first base of this match
        idx_anchor = idx_codon if flag != "z" else idx_codon + 2
        upstream_len = idx_anchor                 # bases left of anchor
        pos = genomic_pos_from_index(
            idx_anchor, anchor_genomic, strand, upstream_len
        )
        if pos == anchor_genomic:
            best_idx_anchor = idx_anchor
            upstream_len_correct = upstream_len
            break

    if best_idx_anchor is None:
        # Fallback: use occurrence closest to middle (should never happen)
        matches = [m.start() for m in re.finditer(re.escape(codon), seq)]
        if not matches:
            raise ValueError("Codon not found in context_sequence.")
        mid    = len(seq) // 2
        idx_c  = min(matches, key=lambda i: abs(i - mid))
        best_idx_anchor = idx_c if flag != "z" else idx_c + 2
        upstream_len_correct = best_idx_anchor

    # ------------------------------------------------------------------
    # 2. define PAM search windows
    # ------------------------------------------------------------------
    if flag == "z":
        NGG_min, NGG_max = -30, 4
        CCN_min, CCN_max = -5, 30
    else:  # 'a'
        NGG_min, NGG_max = -30, 3
        CCN_min, CCN_max = -3, 30

    NGG_hits, CCN_hits = [], []

    # ------------------------------------------------------------------
    # 3. scan NGG motifs
    # ------------------------------------------------------------------
    for m in re.finditer(r"[ACGT]GG", seq):
        idx = m.start()  # index of the N in NGG
        pos = genomic_pos_from_index(
            idx, anchor_genomic, strand, upstream_len_correct
        )
        offset = pos - anchor_genomic
        if NGG_min <= offset <= NGG_max:
            NGG_hits.append(str(pos))

    # ------------------------------------------------------------------
    # 4. scan CCN motifs
    # ------------------------------------------------------------------
    for m in re.finditer(r"CC[ACGT]", seq):
        idx = m.start()  # first C
        pos = genomic_pos_from_index(
            idx, anchor_genomic, strand, upstream_len_correct
        )
        offset = pos - anchor_genomic
        if CCN_min <= offset <= CCN_max:
            CCN_hits.append(str(pos))

    return ";".join(NGG_hits), ";".join(CCN_hits)


# --------------------------------------------------------------------- #
# main driver
# --------------------------------------------------------------------- #
def main():
    import argparse

    ap = argparse.ArgumentParser(description="Find NGG / CCN PAMs around codons")
    ap.add_argument("--infile",  default="annotated_scores.csv",
                    help="input CSV from annotate.py")
    ap.add_argument("--outfile", default="annotated_scores_with_PAMs.csv",
                    help="output CSV with PAM columns added")
    args = ap.parse_args()

    df = pd.read_csv(args.infile)

    # required columns check
    needed = {"NeighborFlag", "context_sequence",
              "codon_genomic_base1", "codon_genomic_base3",
              "codon_sequence", "strand"}
    miss = needed - set(df.columns)
    if miss:
        sys.exit(f"ERROR: missing columns: {', '.join(miss)}")

    pam_df = df.apply(find_pams_in_row, axis=1, result_type="expand")
    pam_df.columns = ["NGG_PAMs", "CCN_PAMs"]

    out_df = pd.concat([df, pam_df], axis=1)
    out_df.to_csv(args.outfile, index=False)
    print(f"Wrote {len(out_df)} rows to {args.outfile}")


if __name__ == "__main__":
    main()
