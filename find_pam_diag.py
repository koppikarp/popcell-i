#!/usr/bin/env python3
"""
find_pam_diag.py

Diagnostic script:
  • Reads select_positions CSV -> gets Gene and the largest Position
  • Reads gene-annotated FASTA (uppercase = CDS) -> translates each record
  • Prints CDS length vs. expected length
  • Shows a pairwise global alignment (identity %) with the canonical
    UniProt sequence if available.

Requires Biopython (pairwise2) and requests.
"""

import argparse
import sys
import re
from pathlib import Path
from typing import List, Dict

import pandas as pd
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
import requests

# ─────────────────────────── regex helpers ────────────────────────────
REFSEQ_RE = re.compile(r'(?:NM|NR|XM|XR|NP|XP)_[0-9]+(?:\.[0-9]+)?')
GENE_RE   = re.compile(r'^>(\S+)')               # first token of header (gene symbol)
CODING_RE = re.compile(r'[A-Z]')                 # uppercase = CDS base

# ─────────────────────────── FASTA parsing ────────────────────────────
def parse_fasta(path: Path) -> List[Dict]:
    recs = []
    with path.open() as fh:
        header, chunks = None, []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    recs.append(_make_rec(header, chunks))
                header, chunks = line[1:], []
            else:
                chunks.append(line)
        if header:
            recs.append(_make_rec(header, chunks))
    return recs


def _make_rec(header: str, chunks: List[str]) -> Dict:
    seq = "".join(chunks)
    gene  = GENE_RE.match(">" + header).group(1)
    ref   = REFSEQ_RE.search(header).group(0)
    # collect uppercase (CDS) bases
    cds_nt = "".join(b for b in seq if b.isupper())
    cds_aa = str(Seq(cds_nt).translate(to_stop=False))
    return {"gene": gene,
            "refseq": ref,
            "header": header,
            "cds_nt": cds_nt,
            "cds_aa": cds_aa}


# ─────────────────────────── UniProt helper ────────────────────────────
def fetch_uniprot_sequence(accession: str) -> str | None:
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        lines = r.text.splitlines()
        return "".join(lines[1:])  # drop header
    except Exception:
        return None


# ─────────────────────────── main ────────────────────────────
def main() -> None:
    ap = argparse.ArgumentParser("Diagnostic length/alignment checker")
    ap.add_argument("-f", "--fasta", required=True, type=Path,
                    help="FASTA with gene-prefixed headers (uppercase=CDS)")
    ap.add_argument("-p", "--positions", required=True, type=Path,
                    help="CSV from select_positions.py")
    args = ap.parse_args()

    # 1. CSV → expected gene + length
    df = pd.read_csv(args.positions)
    gene = df["Gene"].iloc[0]
    expected_len = int(df["Position"].max())          # biggest position

    # 2. fetch canonical UniProt protein (one accession per CSV)
    protein_acc = df["Protein"].iloc[0]
    canonical_seq = fetch_uniprot_sequence(protein_acc)
    if canonical_seq:
        sys.stderr.write(f"[INFO] fetched UniProt sequence {protein_acc} "
                         f"({len(canonical_seq)} aa)\n")
    else:
        sys.stderr.write(f"[WARN] could not fetch UniProt seq {protein_acc}; "
                         f"using {expected_len} × 'X'\n")
        canonical_seq = "X" * expected_len

    # 3. iterate fasta records
    for rec in parse_fasta(args.fasta):
        if rec["gene"] != gene:
            continue

        cds_len = len(rec["cds_aa"])
        print(f"\n=== {rec['refseq']}  Gene={gene} ===")
        print(f"  CDS length      : {cds_len}")
        print(f"  Expected length : {expected_len}")

        # pairwise alignment (global, identity scoring = xx)
        align = pairwise2.align.globalxx(canonical_seq, rec["cds_aa"], one_alignment_only=True)[0]
        identity = align.score / max(len(canonical_seq), len(rec["cds_aa"])) * 100
        print(f"  Identity        : {identity:.1f}%")
        print(format_alignment(*align, full_sequences=True))

    print("\n[DONE]")


if __name__ == "__main__":
    main()
