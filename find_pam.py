#!/usr/bin/env python3
"""
find_pam.py

Link positions chosen in the scoring CSV to their genomic coordinates in a
gene-annotated FASTA file whose coding regions are denoted by uppercase.

Outputs (TSV to stdout):
    Gene  RefSeq  StartBaseIdx  Position  Residue  TotalPenalty
    NeighborFlag  Protein  GeneCopy

Notes
-----
* Handles multiple FASTA records per gene; each valid match is emitted.
* Assumes '+' strand.  Handling the '-' strand would require reverse-complement
  logic that is *not* implemented here.
"""

import argparse
import re
import sys
from pathlib import Path
from typing import List, Dict, Tuple

import pandas as pd

# Translation dictionary (standard genetic code)
CODON_TABLE = {
    **{c: aa for c, aa in zip(
        # 1st position row
        ["TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TAA","TAG","TGT","TGC","TGA","TGG"],
        ["F","F","L","L","S","S","S","S","Y","Y","*","*","C","C","*","W"])},
    **{c: aa for c, aa in zip(
        # 2nd position row
        ["CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG"],
        ["L","L","L","L","P","P","P","P","H","H","Q","Q","R","R","R","R"])},
    **{c: aa for c, aa in zip(
        # 3rd position row
        ["ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG"],
        ["I","I","I","M","T","T","T","T","N","N","K","K","S","S","R","R"])},
    **{c: aa for c, aa in zip(
        # 4th position row
        ["GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG"],
        ["V","V","V","V","A","A","A","A","D","D","E","E","G","G","G","G"])}
}

# Regex tools
REFSEQ_RE = re.compile(r'(?:NM|NR|XM|XR|NP|XP)_[0-9]+(?:\.[0-9]+)?')
GENE_RE   = re.compile(r'^>(\S+)')               # first token of header (gene symbol)


# ---------------------------------------------------------------------------#
# Helpers
# ---------------------------------------------------------------------------#
def translate_codon(codon: str) -> str:
    """Translate an uppercase 3-nt codon using the standard code."""
    return CODON_TABLE.get(codon.upper(), "X")   # X for any ambiguous/unknown


def parse_fasta(path: Path) -> List[Dict]:
    """
    Read a FASTA with gene-prefixed headers.
    Returns a list of dicts, each containing:
        gene       – gene symbol (1st token)
        refseq     – RefSeq accession (first NM_ / NR_ … in header)
        header     – full header line (without '>')
        seq        – full sequence (newlines removed)
        coding_idx – list[int] indices (1-based) of every uppercase base in seq
    """
    records: List[Dict] = []
    with path.open() as fh:
        header = None
        seq_chunks: List[str] = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append(_build_record(header, seq_chunks))
                header  = line[1:]          # drop '>'
                seq_chunks = []
            else:
                seq_chunks.append(line)
        # final record
        if header is not None:
            records.append(_build_record(header, seq_chunks))
    return records


def _build_record(header: str, seq_chunks: List[str]) -> Dict:
    seq = "".join(seq_chunks)
    gene_m = GENE_RE.match(">" + header)      # re-add '>' so regex works
    gene = gene_m.group(1) if gene_m else "UNKNOWN"

    ref_m = REFSEQ_RE.search(header)
    refseq = ref_m.group(0) if ref_m else "NA"

    coding_idx: List[int] = []
    for i, base in enumerate(seq, start=1):   # 1-based
        if base.isupper():
            coding_idx.append(i)

    # quick sanity: length multiple of 3?
    if len(coding_idx) % 3 != 0:
        sys.stderr.write(f"[WARN] coding length in {header} not multiple of 3 "
                         f"({len(coding_idx)} nt)\n")

    return {"gene": gene,
            "refseq": refseq,
            "header": header,
            "seq": seq,
            "coding_idx": coding_idx}


def find_codon_start(record: Dict, aa_pos: int) -> Tuple[int, str]:
    """
    Given a FASTA record and a 1-based amino-acid position, return
    (start_idx_in_original_seq, codon_string).

    Raises IndexError if aa_pos exceeds CDS length.
    """
    idx = record["coding_idx"]
    if aa_pos * 3 > len(idx):
        raise IndexError("AA position beyond coding length")

    start_idx = idx[(aa_pos - 1) * 3]          # first nt of codon
    codon_nts = [record["seq"][i - 1] for i in idx[(aa_pos - 1) * 3 : aa_pos * 3]]
    codon = "".join(codon_nts)
    return start_idx, codon


# ---------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------#
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Map scored amino-acid positions onto genomic FASTA coordinates."
    )
    parser.add_argument("-f", "--fasta",   required=True, type=Path, help="FASTA with gene-prefixed headers")
    parser.add_argument("-p", "--positions", required=True, type=Path, help="select_positions_with_gene CSV")
    args = parser.parse_args()

    # ------------------------------------------------------------------- #
    # 1. Load FASTA records
    # ------------------------------------------------------------------- #
    sys.stderr.write("[INFO] parsing FASTA …\n")
    fasta_records = parse_fasta(args.fasta)
    # Group by gene for quick lookup
    by_gene: Dict[str, List[Dict]] = {}
    for rec in fasta_records:
        by_gene.setdefault(rec["gene"], []).append(rec)

    # ------------------------------------------------------------------- #
    # 2. Read positions CSV
    # ------------------------------------------------------------------- #
    df = pd.read_csv(args.positions)
    required_cols = {"Position", "Residue", "TotalPenalty",
                     "NeighborFlag", "Protein", "Gene"}
    if not required_cols.issubset(df.columns):
        sys.exit(f"Error: CSV must contain {', '.join(required_cols)}")

    # ------------------------------------------------------------------- #
    # 3. Iterate rows and map to each matching FASTA record
    # ------------------------------------------------------------------- #
    out_rows = []
    for _, row in df.iterrows():
        gene = row["Gene"]
        aa_pos = int(row["Position"])
        residue_expected = row["Residue"]

        if gene not in by_gene:
            sys.stderr.write(f"[WARN] gene {gene} not found in FASTA – skipping\n")
            continue

        for rec in by_gene[gene]:
            try:
                cds_start_idx, codon = find_codon_start(rec, aa_pos)
            except IndexError:
                sys.stderr.write(f"[WARN] {gene}/{rec['refseq']} too short "
                                 f"for AA position {aa_pos}\n")
                continue

            residue_actual = translate_codon(codon)
            if residue_actual != residue_expected:
                sys.stderr.write(
                    f"[WARN] {gene}/{rec['refseq']} AA mismatch at pos {aa_pos}: "
                    f"{residue_actual} (DNA {codon}) ≠ expected {residue_expected}\n"
                )
                continue      # skip mismatching transcripts

            out_rows.append({
                "Gene": gene,
                "RefSeq": rec["refseq"],
                "StartBaseIdx": cds_start_idx,
                "Position": row["Position"],
                "Residue": row["Residue"],
                "TotalPenalty": row["TotalPenalty"],
                "NeighborFlag": row["NeighborFlag"],
                "Protein": row["Protein"],
                "GeneCopy": gene           # redundant, per request
            })

    # ------------------------------------------------------------------- #
    # 4. Emit TSV
    # ------------------------------------------------------------------- #
    if not out_rows:
        sys.stderr.write("[INFO] no matches produced any output\n")
        return

    out_df = pd.DataFrame(out_rows, columns=[
        "Gene", "RefSeq", "StartBaseIdx",
        "Position", "Residue", "TotalPenalty",
        "NeighborFlag", "Protein", "GeneCopy",
    ])
    out_df.to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    main()
