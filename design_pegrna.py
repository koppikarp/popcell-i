#!/usr/bin/env python3
"""
design_pegrna.py

Choose the closest valid PAM for each row of *annotated_scores_with_PAMs.csv*
and build pegRNA components (protospacer, PBS, RTT).

Implemented cases
-----------------
1. NGG PAM   & NeighborFlag 'a'
2. CCN PAM   & NeighborFlag 'z'
3. NGG PAM   & NeighborFlag 'z'
4. CCN PAM   & NeighborFlag 'a'

Filters
-------
* No component may contain BsaI, BbsI or BsmBI sites (or their RCs).
* Protospacer GC content must be 15–85 %.
"""

import csv
import json
import re
import sys
from pathlib import Path
from pam_surgery import patch_pam
import Bio.Seq as Seq


import pandas as pd

RESTRICTION_6MERS = {
    "GGTCTC", "GAGACC",  # BsaI
    "GAAGAC", "GTCTTC",  # BbsI
    "CGTCTC", "GAGACG",  # BsmBI
}

# add reverse complements
RESTRICTION_6MERS |= {s[::-1].translate(str.maketrans("ACGT", "TGCA"))
                      for s in RESTRICTION_6MERS}


def rc(seq: str) -> str:
    """Reverse-complement (DNA)."""
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def gc_fraction(seq: str) -> float:
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq)


def genomic_to_index(pos: int,
                     anchor_genomic: int,
                     upstream_len: int,
                     strand: int) -> int:
    """
    Convert absolute hg38 coordinate → index in context_sequence.
    anchor_genomic = genomic coord of first (flag a) or third (flag z) codon base.
    upstream_len   = bases left of anchor inside context string.
    """
    if strand == 1:
        ctx_start = anchor_genomic - upstream_len
        return pos - ctx_start
    else:
        ctx_end = anchor_genomic + upstream_len
        return ctx_end - pos


def locate_anchor(row) -> tuple[int, int]:
    """
    Return (anchor_index_in_seq, upstream_len)
    """
    seq = row["context_sequence"].upper()
    codon = row["codon_sequence"].upper()
    flag = row["NeighborFlag"].lower()

    # find the codon match closest to centre
    matches = [m.start() for m in re.finditer(re.escape(codon), seq)]
    if not matches:
        raise ValueError("Codon not found in context_sequence")

    mid = len(seq) // 2
    codon_idx = min(matches, key=lambda i: abs(i - mid))

    anchor_idx = codon_idx if flag == "a" else codon_idx + 2
    return anchor_idx, anchor_idx  # upstream_len == anchor_idx


def build_components(row) -> dict:
    """
    Return dict with keys:
      chosen_PAM_pos, PAM_type, protospacer, PBS, RTT
    or empty strings if no valid PAM found.
    """
    seq = row["context_sequence"].upper()
    flag = row["NeighborFlag"].lower()
    strand = int(row["strand"])

    # anchor genomic coordinate
    anchor_genomic = int(row["codon_genomic_base1"] if flag == "a"
                         else row["codon_genomic_base3"])

    # locate anchor inside seq
    anchor_idx, upstream_len = locate_anchor(row)

    # build dict of genomic_pos → index_in_seq for all PAMs
    pam_positions = []
    for pam_type, col in (("NGG", "NGG_PAMs"), ("CCN", "CCN_PAMs")):
        for pos_str in filter(None, row[col].split(";")):
            gpos = int(pos_str)
            idx = genomic_to_index(gpos, anchor_genomic, upstream_len, strand)
            pam_positions.append((pam_type, gpos, idx))

    # sort by |distance|, tie-break upstream on + strand
    pam_positions.sort(key=lambda tup: (abs(tup[1] - anchor_genomic),
                                        tup[1] if strand == 1 else -tup[1]))

    for pam_type, gpos, idx in pam_positions:
        # pick only PAMs compatible with this flag
        if   flag == "a" and pam_type != "NGG":
            continue
        elif flag == "z" and pam_type != "CCN":
            continue
        elif flag not in ("a", "z"):
            continue  # other flags: not handled

        if pam_type == "NGG":     # ---------- NGG protospacer / PBS
            N_idx = idx
            protospacer = seq[N_idx - 20:N_idx]
            pbs = rc(seq[N_idx - 18:N_idx - 3])

            # RTT and PAM mod logic
            if flag == "a":
                cod_start = anchor_idx
                distance = abs(N_idx - cod_start)
                sign = '-' if (N_idx - cod_start) < 0 else '+'
                phase = distance % 3
                if phase == 0:
                    window = seq[N_idx:N_idx + 3]
                    seq = seq[:N_idx] + patch_pam(window, pam_type) + seq[N_idx+3:]
                elif phase == 1:
                    if sign == '-':
                        window = seq[N_idx - 2:N_idx + 4]
                        seq = seq[:N_idx-2] + patch_pam(window, pam_type) + seq[N_idx+4:]
                    else: # sign == '+'
                        window = seq[N_idx - 1:N_idx + 5]
                        seq = seq[:N_idx-1] + patch_pam(window, pam_type) + seq[N_idx+5:]
                elif phase == 2:
                    if sign == '-':
                        window = seq[N_idx - 1:N_idx + 5]
                        seq = seq[:N_idx-1] + patch_pam(window, pam_type) + seq[N_idx+5:]
                    else: # sign == '+'
                        window = seq[N_idx - 2:N_idx + 4]
                        seq = seq[:N_idx-2] + patch_pam(window, pam_type) + seq[N_idx+4:]

                rtt1 = seq[N_idx - 3: cod_start]
                rtt2 = "NNN"
                rtt3 = seq[cod_start: cod_start + 27]
                rtt = rc(rtt1 + rtt2 + rtt3)
            else:  # NGG / z
                cod_end = anchor_idx
                cod_start = cod_end - 2
                distance = abs(N_idx - cod_start)
                sign = '-' if (N_idx - cod_start) < 0 else '+'
                phase = distance % 3
                if phase == 0:
                    window = seq[N_idx:N_idx + 3]
                    seq = seq[:N_idx] + patch_pam(window, pam_type) + seq[N_idx+3:]
                elif phase == 1:
                    if sign == '-':
                        window = seq[N_idx - 2:N_idx + 4]
                        seq = seq[:N_idx-2] + patch_pam(window, pam_type) + seq[N_idx+4:]
                    else: # sign == '+'
                        window = seq[N_idx - 1:N_idx + 5]
                        seq = seq[:N_idx-1] + patch_pam(window, pam_type) + seq[N_idx+5:]
                elif phase == 2:
                    if sign == '-':
                        window = seq[N_idx - 1:N_idx + 5]
                        seq = seq[:N_idx-1] + patch_pam(window, pam_type) + seq[N_idx+5:]
                    else: # sign == '+'
                        window = seq[N_idx - 2:N_idx + 4]
                        seq = seq[:N_idx-2] + patch_pam(window, pam_type) + seq[N_idx+4:]
                rtt1 = seq[N_idx - 3:cod_end + 1]
                rtt2 = "NNN"
                rtt3 = seq[cod_end + 1:cod_end + 28]
                rtt = rc(rtt1 + rtt2 + rtt3)

        else:                     # ---------- CCN protospacer / PBS
            C_idx = idx
            protospacer = rc(seq[C_idx + 3:C_idx + 23])    # 20 nt
            pbs = rc(rc(seq[C_idx + 6:C_idx + 21]))         # 15 nt

            # RTT and PAM mod logic
            if flag == "z":
                cod_end = anchor_idx
                cod_start = cod_end - 2
                distance = abs(C_idx - cod_start)
                sign = '-' if (C_idx - cod_start) < 0 else '+'
                phase = distance % 3
                if phase == 0:
                    window = seq[C_idx:C_idx + 3]
                    seq = seq[:C_idx] + patch_pam(window, pam_type) + seq[C_idx+3:]
                elif phase == 1:
                    if sign == '-':
                        window = seq[C_idx - 2:C_idx + 4]
                        seq = seq[:C_idx-2] + patch_pam(window, pam_type) + seq[C_idx+4:]
                    else: # sign == '+'
                        window = seq[C_idx - 1:C_idx + 5]
                        seq = seq[:C_idx-1] + patch_pam(window, pam_type) + seq[C_idx+5:]
                elif phase == 2:
                    if sign == '-':
                        window = seq[C_idx - 1:C_idx + 5]
                        seq = seq[:C_idx-1] + patch_pam(window, pam_type) + seq[C_idx+5:]
                    else: # sign == '+'
                        window = seq[C_idx - 2:C_idx + 4]
                        seq = seq[:C_idx-2] + patch_pam(window, pam_type) + seq[C_idx+4:]
                rtt1 = rc(seq[cod_end + 1:C_idx + 6])
                rtt2 = "NNN"
                rtt3 = rc(seq[cod_end - 26: cod_end + 1])
                rtt = rc(rtt1 + rtt2 + rtt3)
            else:  # CCN / a
                cod_start = anchor_idx
                distance = abs(C_idx - cod_start)
                sign = '-' if (C_idx - cod_start) < 0 else '+'
                phase = distance % 3
                if phase == 0:
                    window = seq[C_idx:C_idx + 3]
                    seq = seq[:C_idx] + patch_pam(window, pam_type) + seq[C_idx+3:]
                elif phase == 1:
                    if sign == '-':
                        window = seq[C_idx - 2:C_idx + 4]
                        seq = seq[:C_idx-2] + patch_pam(window, pam_type) + seq[C_idx+4:]
                    else: # sign == '+'
                        window = seq[C_idx - 1:C_idx + 5]
                        seq = seq[:C_idx-1] + patch_pam(window, pam_type) + seq[C_idx+5:]
                elif phase == 2:
                    if sign == '-':
                        window = seq[C_idx - 1:C_idx + 5]
                        seq = seq[:C_idx-1] + patch_pam(window, pam_type) + seq[C_idx+5:]
                    else: # sign == '+'
                        window = seq[C_idx - 2:C_idx + 4]
                        seq = seq[:C_idx-2] + patch_pam(window, pam_type) + seq[C_idx+4:]
                rtt1 = rc(seq[cod_start:C_idx + 6])
                rtt2 = "NNN"
                rtt3 = rc(seq[cod_start - 27:cod_start])
                rtt = rc(rtt1 + rtt2 + rtt3)

        # length sanity
        if len(protospacer) != 20 or len(pbs) != 15:
            continue  # skip malformed slices

        # GC filter
        gc = gc_fraction(protospacer)
        if not (0.15 <= gc <= 0.85):
            continue

        # restriction-site filter
        bad = any(site in part for part in (protospacer, pbs, rtt)
                  for site in RESTRICTION_6MERS if site)
        if bad:
            continue

        # passed all filters
        return {
            "chosen_PAM_pos": str(gpos),
            "PAM_type": pam_type,
            "protospacer": protospacer,
            "PBS": pbs,
            "RTT": rtt,
        }

    # no valid PAM found
    return {
        "chosen_PAM_pos": "",
        "PAM_type": "",
        "protospacer": "",
        "PBS": "",
        "RTT": "",
    }


def main():
    import argparse

    ap = argparse.ArgumentParser(description="Design pegRNA components")
    ap.add_argument("--infile", default="annotated_scores_with_PAMs.csv",
                    help="input CSV (default: annotated_scores_with_PAMs.csv)")
    ap.add_argument("--outfile", default="pegRNA_designs.csv",
                    help="output CSV (default: pegRNA_designs.csv)")
    args = ap.parse_args()

    df = pd.read_csv(args.infile)
    comps = df.apply(build_components, axis=1, result_type="expand")
    result = pd.concat([df, comps], axis=1)
    result.to_csv(args.outfile, index=False)
    print(f"✓ Wrote {len(result)} rows to {args.outfile}")


if __name__ == "__main__":
    main()
