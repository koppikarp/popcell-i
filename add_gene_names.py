#!/usr/bin/env python3
"""
add_gene_names.py

Usage
-----
python add_gene_names.py input.fasta [-o output.fasta]

Function
--------
For every FASTA header that contains a RefSeq accession (NM_, NR_, XM_, XR_,
NP_, XP_), look up the human gene symbol via MyGene.info and prepend it:

    >GENE_SYM original header ...

If the lookup fails, the accession itself is used in place of the gene symbol.
All sequence lines are passed through unchanged.

The script caches look-ups so each accession is queried only once.
"""

import argparse
import re
import sys
from pathlib import Path
from functools import lru_cache
import requests


# Regex to capture common RefSeq accessions with optional version (e.g. .2)
REFSEQ_RE = re.compile(
    r'(?:NM|NR|XM|XR|NP|XP)_[0-9]+(?:\.[0-9]+)?',
    re.IGNORECASE,
)

MYGENE_URL = "https://mygene.info/v3/query"


@lru_cache(maxsize=None)
def refseq_to_gene_symbol(refseq: str) -> str:
    """
    Query MyGene.info for a human gene symbol given a RefSeq accession.
    Falls back to the accession itself if nothing is found or on error.
    """
    params = {
    "q": refseq,
    "scopes": "refseq",   # <— add this
    "fields": "symbol",
    "species": "human",
    "size": 1,
}
    try:
        r = requests.get(MYGENE_URL, params=params, timeout=10)
        r.raise_for_status()
        data = r.json()
        if data.get("total", 0):
            return data["hits"][0].get("symbol", refseq)
    except Exception:
        pass
    return refseq


def process_fasta(in_path: Path, out_handle) -> None:
    """
    Read FASTA from *in_path* and write modified headers and sequences to *out_handle*.
    """
    with in_path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")  # strip leading '>' and newline
                match = REFSEQ_RE.search(header)
                if match:
                    refseq = match.group(0)
                    gene   = refseq_to_gene_symbol(refseq)
                    out_handle.write(f">{gene} {header}\n")
                else:
                    sys.stderr.write(f"[WARN] no RefSeq ID found in header: {header}\n")
                    out_handle.write(line)

            else:
                out_handle.write(line)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Add gene symbols to FASTA headers using RefSeq accessions."
    )
    parser.add_argument("fasta", type=Path, help="Input FASTA file")
    parser.add_argument(
        "-o", "--output", type=Path, help="Output file (default: stdout)", default=None
    )
    args = parser.parse_args()

    # Validate input file exists
    if not args.fasta.exists():
        sys.exit(f"Error: {args.fasta} not found.")

    # Decide output destination
    if args.output:
        with args.output.open("w") as out_fh:
            process_fasta(args.fasta, out_fh)
    else:
        process_fasta(args.fasta, sys.stdout)


if __name__ == "__main__":
    main()
