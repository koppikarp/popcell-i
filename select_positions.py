#!/usr/bin/env python3
"""
select_positions.py

Select first- and last-position rows plus the top-n best-scoring rows from a
ranked penalty CSV, add NeighborFlag/Protein/Gene, and write the result to a
new CSV file that can be consumed by find_pam.py
"""

import argparse
import os
import sys
import pandas as pd
import requests

# ─────────────────────────── helpers ────────────────────────────
def get_accession_and_gene(csv_path: str) -> tuple[str, str]:
    fname = os.path.basename(csv_path)
    protein_acc = os.path.splitext(fname)[0].split("_")[0]

    url = f"https://rest.uniprot.org/uniprotkb/{protein_acc}.json"
    try:
        data = requests.get(url, timeout=10).json()
        gene_name = data.get("genes", [{}])[0].get("geneName", {}).get("value", protein_acc)
    except Exception:
        gene_name = protein_acc
    return protein_acc, gene_name


def compute_neighbor_flags(df: pd.DataFrame) -> pd.Series:
    ordered = df.sort_values("Position").reset_index(drop=True)
    ordered["PrevPenalty"] = ordered["TotalPenalty"].shift(1)
    ordered["NextPenalty"] = ordered["TotalPenalty"].shift(-1)

    def flag(row):
        if pd.isna(row["PrevPenalty"]):   # first row
            return "z"
        if pd.isna(row["NextPenalty"]):   # last row
            return "a"
        return "a" if row["PrevPenalty"] <= row["NextPenalty"] else "z"

    return ordered.apply(flag, axis=1).set_axis(ordered["Position"])


def select_positions(df: pd.DataFrame, n: int, protein: str, gene: str) -> pd.DataFrame:
    min_pos, max_pos = df["Position"].min(), df["Position"].max()
    first_row = df.loc[df["Position"] == min_pos].iloc[0]
    last_row  = df.loc[df["Position"] == max_pos].iloc[0]

    ranked = (
        df.sort_values("TotalPenalty")
          .query("Position not in (@min_pos, @max_pos)")
          .head(n)
    )

    # if you want the first and last positions too (but just use normal popcell for this)
    result = pd.concat([first_row.to_frame().T,
                        last_row.to_frame().T,
                        ranked], ignore_index=True)
    
    # for only internal positions)
    result = pd.concat([ranked], ignore_index=True)

    result["NeighborFlag"] = result["Position"].map(compute_neighbor_flags(df))
    result["Protein"] = protein
    result["Gene"] = gene
    return result[["Position", "Residue", "TotalPenalty",
                   "NeighborFlag", "Protein", "Gene"]]

# ──────────────────────────── CLI ───────────────────────────────
def main() -> None:
    p = argparse.ArgumentParser(
        description="Create a select-positions CSV ready for find_pam.py")
    p.add_argument("csv_file",
                   help="Input CSV ranked by TotalPenalty (e.g. P53367_scores.csv)")
    p.add_argument("-n", type=int, required=True,
                   help="How many best-scoring rows to keep (besides first & last)")
    p.add_argument("-o", "--out", default=None,
                   help="Output CSV file (default: <protein>_selected_positions.csv)")
    args = p.parse_args()

    try:
        df = pd.read_csv(args.csv_file)
    except Exception as exc:
        sys.exit(f"Error reading CSV: {exc}")

    needed = {"Position", "Residue", "TotalPenalty"}
    if not needed.issubset(df.columns):
        sys.exit(f"CSV must contain columns: {', '.join(needed)}")
    if args.n <= 0 or args.n > len(df) - 2:
        sys.exit("Argument -n must be >0 and leave at least two rows for first/last.")

    protein_acc, gene_name = get_accession_and_gene(args.csv_file)
    out_df = select_positions(df, args.n, protein_acc, gene_name)

    out_path = args.out or f"{protein_acc}_selected_positions.csv"
    out_df.to_csv(out_path, index=False)
    print(f"[DONE] wrote {len(out_df)} rows to {out_path}")

if __name__ == "__main__":
    main()
