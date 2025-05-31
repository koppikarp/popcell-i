#!/usr/bin/env python3
"""
select_positions.py

Return (n + 2) rows from a ranked-by-TotalPenalty CSV:

 • Row with the smallest Position  (first residue)   → flag 'z'
 • Row with the largest  Position  (last  residue)   → flag 'a'
 • The n best-scoring rows *excluding* those two

Add column 'NeighborFlag' = 'a' or 'z' depending on whether the
closest-scoring lower-penalty neighbor is in front ('a') or behind ('z').
"""

import argparse
import sys
import pandas as pd


def compute_neighbor_flags(df: pd.DataFrame) -> pd.Series:
    """
    For every row, label 'a' if the closest neighbor with lower
    TotalPenalty is earlier (smaller Position); else 'z'.
    """
    # Ensure rows are ordered by Position so .shift() yields neighbors
    ordered = df.sort_values("Position").reset_index(drop=True)

    # Penalty of previous / next positions in sequence
    ordered["PrevPenalty"] = ordered["TotalPenalty"].shift(1)
    ordered["NextPenalty"] = ordered["TotalPenalty"].shift(-1)

    def flag(row):
        prev_pen = row["PrevPenalty"]
        next_pen = row["NextPenalty"]

        # Decide when only one neighbor exists
        if pd.isna(prev_pen):
            return "z"   # only next neighbor
        if pd.isna(next_pen):
            return "a"   # only previous neighbor

        # Both neighbors exist – choose the lower-penalty one
        return "a" if prev_pen <= next_pen else "z"

    ordered["NeighborFlag"] = ordered.apply(flag, axis=1)

    # Map from Position to flag for fast lookup later
    return ordered.set_index("Position")["NeighborFlag"]


def select_positions(df: pd.DataFrame, n: int) -> pd.DataFrame:
    # Identify first and last positions in the sequence
    min_pos = df["Position"].min()
    max_pos = df["Position"].max()

    first_row = df.loc[df["Position"] == min_pos].iloc[0]
    last_row  = df.loc[df["Position"] == max_pos].iloc[0]

    # Best n rows excluding first/last positions
    ranked = (
        df.sort_values("TotalPenalty", ascending=True)
          .query("Position not in (@min_pos, @max_pos)")
          .head(n)
    )

    # Assemble result DataFrame
    result = pd.concat(
        [first_row.to_frame().T, last_row.to_frame().T, ranked],
        ignore_index=True,
    )

    # Add NeighborFlag column
    neighbor_flags = compute_neighbor_flags(df)
    result["NeighborFlag"] = result["Position"].map(neighbor_flags)

    return result[["Position", "Residue", "TotalPenalty", "NeighborFlag"]]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Select first & last positions plus the top-n best-scoring positions"
                    " and label each with an 'a' or 'z' neighbor flag."
    )
    parser.add_argument("csv_file", help="Path to ranked CSV file")
    parser.add_argument(
        "-n", type=int, required=True,
        help="Number of additional top positions to return (excluding first & last)"
    )
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.csv_file)
    except Exception as exc:
        sys.exit(f"Error reading CSV: {exc}")

    needed = {"Position", "Residue", "TotalPenalty"}
    if not needed.issubset(df.columns):
        sys.exit(f"CSV must contain columns: {', '.join(needed)}")

    if args.n <= 0:
        sys.exit("n must be a positive integer")

    available_core = len(df) - 2
    if args.n > available_core:
        sys.exit(
            f"n={args.n} exceeds the available positions ({available_core}) "
            "after excluding first and last."
        )

    out = select_positions(df, args.n)
    out.to_csv(sys.stdout, index=False)


if __name__ == "__main__":
    main()
