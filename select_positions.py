#!/usr/bin/env python3
"""
select_positions.py

Given a CSV ranked by TotalPenalty, extract the protein accession from the
filename (e.g. "P53367_scores.csv" → protein "P53367"), look up its corresponding
gene name via UniProt, and then return:

  • The row whose Position is the minimum (first residue)   → NeighborFlag 'z'
  • The row whose Position is the maximum (last  residue)   → NeighborFlag 'a'
  • The n best-scoring rows (lowest TotalPenalty) excluding those two

Output columns:
  Position, Residue, TotalPenalty, NeighborFlag, Protein, Gene

If UniProt lookup fails, the gene name will default to the accession itself.
"""

import argparse
import os
import sys
import pandas as pd
import requests


def get_accession_and_gene(csv_path: str) -> (str, str):
    """
    Derive the protein accession from the CSV filename (before the first underscore),
    then query UniProt to fetch the gene name. If the lookup fails, return the accession
    as the gene name.

    Example:
      "/full/path/to/P53367_scores.csv"  → protein_acc = "P53367"
      UniProt lookup → gene_name = "ABC1" (if found), else gene_name = "P53367"
    """
    fname = os.path.basename(csv_path)
    # Strip extension and then split on underscore
    name_root, _ = os.path.splitext(fname)
    protein_acc = name_root.split("_")[0]

    # UniProt REST API (newer endpoint)
    url = f"https://rest.uniprot.org/uniprotkb/{protein_acc}.json"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        data = r.json()
        # Under "genes", pick the first geneName.value if present
        gene_entry = data.get("genes", [{}])[0]
        gene_name = gene_entry.get("geneName", {}).get("value", protein_acc)
    except Exception:
        gene_name = protein_acc

    return protein_acc, gene_name


def compute_neighbor_flags(df: pd.DataFrame) -> pd.Series:
    """
    For every row (sorted by Position), assign:
      • 'z' if the only neighbor with lower TotalPenalty is behind (larger Position),
      • 'a' if the only neighbor with lower TotalPenalty is in front (smaller Position),
      • If both exist, compare PrevPenalty vs. NextPenalty:
          – 'a' if PrevPenalty <= NextPenalty (i.e. in‐front neighbor is lower or tie)
          – 'z' otherwise
    """
    ordered = df.sort_values("Position").reset_index(drop=True)
    ordered["PrevPenalty"] = ordered["TotalPenalty"].shift(1)
    ordered["NextPenalty"] = ordered["TotalPenalty"].shift(-1)

    def flag(row):
        prev_pen = row["PrevPenalty"]
        next_pen = row["NextPenalty"]

        if pd.isna(prev_pen):
            # No previous row → only next exists → flag 'z'
            return "z"
        if pd.isna(next_pen):
            # No next row → only previous exists → flag 'a'
            return "a"
        # Both neighbors exist
        return "a" if prev_pen <= next_pen else "z"

    ordered["NeighborFlag"] = ordered.apply(flag, axis=1)
    return ordered.set_index("Position")["NeighborFlag"]


def select_positions(df: pd.DataFrame, n: int, protein: str, gene: str) -> pd.DataFrame:
    """
    – Identify the first (min Position) and last (max Position) rows.
    – Take the top‐n best‐scoring rows (lowest TotalPenalty) excluding those two.
    – Compute NeighborFlag for each returned row.
    – Add Protein and Gene columns (same value on every row).
    """
    min_pos = df["Position"].min()
    max_pos = df["Position"].max()

    first_row = df.loc[df["Position"] == min_pos].iloc[0]
    last_row  = df.loc[df["Position"] == max_pos].iloc[0]

    ranked = (
        df.sort_values("TotalPenalty", ascending=True)
          .query("Position not in (@min_pos, @max_pos)")
          .head(n)
    )

    result = pd.concat(
        [first_row.to_frame().T, last_row.to_frame().T, ranked],
        ignore_index=True,
    )

    # Compute flags for *all* positions, then map onto result
    neighbor_flags = compute_neighbor_flags(df)
    result["NeighborFlag"] = result["Position"].map(neighbor_flags)
    # Add Protein and Gene columns
    result["Protein"] = protein
    result["Gene"] = gene

    # Reorder columns for clarity
    return result[
        ["Position", "Residue", "TotalPenalty", "NeighborFlag", "Protein", "Gene"]
    ]


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Select the first & last positions plus the top-n best-scoring positions from "
            "a ranked CSV, and label each with 'a'/'z' neighbor flags. "
            "Also prepend Protein/Gene (derived from filename via UniProt)."
        )
    )
    parser.add_argument("csv_file", help="Path to the ranked CSV file (e.g. P53367_scores.csv)")
    parser.add_argument(
        "-n",
        type=int,
        required=True,
        help="Number of additional top positions to return (excluding first & last)",
    )
    args = parser.parse_args()

    # 1) Read CSV
    try:
        df = pd.read_csv(args.csv_file)
    except Exception as e:
        sys.exit(f"Error reading CSV: {e!s}")

    # 2) Validate required columns
    required_cols = {"Position", "Residue", "TotalPenalty"}
    if not required_cols.issubset(df.columns):
        sys.exit(f"CSV must contain columns: {', '.join(required_cols)}")

    # 3) Check n is valid
    if args.n <= 0:
        sys.exit("Argument -n must be a positive integer.")
    if args.n > (len(df) - 2):
        sys.exit(
            f"-n = {args.n} is too large: only {len(df)-2} positions remain once first/last are excluded."
        )

    # 4) Derive protein accession & gene name
    protein_acc, gene_name = get_accession_and_gene(args.csv_file)

    # 5) Select rows + compute flags + append Protein/Gene
    out_df = select_positions(df, args.n, protein_acc, gene_name)

    # 6) Write to stdout as CSV
    out_df.to_csv(sys.stdout, index=False)


if __name__ == "__main__":
    main()
