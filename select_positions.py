#!/usr/bin/env python3
"""
select_positions.py

Given a CSV ranked from the smallest to largest TotalPenalty,
return n + 2 rows:

  • The row whose Position is the minimum value in the file
  • The row whose Position is the maximum value in the file
  • The n best-scoring rows (lowest TotalPenalty) *excluding*
    those two positions

Only the columns Position, Residue and TotalPenalty are printed.
"""

import argparse
import sys
import pandas as pd


def select_positions(df: pd.DataFrame, n: int) -> pd.DataFrame:
    """
    Parameters
    ----------
    df : pandas.DataFrame
        Must contain columns 'Position', 'Residue', 'TotalPenalty'.
    n : int
        Number of top-ranked rows to return in addition to the
        first and last position.

    Returns
    -------
    pandas.DataFrame
        n + 2 rows with columns Position, Residue, TotalPenalty.
        Ordered as [min-position row, max-position row, top-n rows].
    """
    # Identify first and last positions in the sequence
    min_pos = df["Position"].min()
    max_pos = df["Position"].max()

    first_row = df.loc[df["Position"] == min_pos].iloc[0]
    last_row  = df.loc[df["Position"] == max_pos].iloc[0]

    # Start from best (lowest) TotalPenalty and exclude first/last positions
    ranked = (
        df.sort_values("TotalPenalty", ascending=True)
          .query("Position not in (@min_pos, @max_pos)")
          .head(n)
    )

    # Assemble result: first, last, then the n best
    result = pd.concat(
        [first_row.to_frame().T, last_row.to_frame().T, ranked],
        ignore_index=True,
    )

    return result[["Position", "Residue", "TotalPenalty"]]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Select the first & last positions plus the top-n best-scoring positions."
    )
    parser.add_argument("csv_file", help="Path to the ranked CSV file")
    parser.add_argument(
        "-n",
        type=int,
        required=True,
        help="Number of additional top positions to return (excluding first & last)",
    )
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.csv_file)
    except Exception as exc:
        sys.exit(f"Error reading CSV: {exc}")

    needed_cols = {"Position", "Residue", "TotalPenalty"}
    if not needed_cols.issubset(df.columns):
        sys.exit(f"CSV must contain columns: {', '.join(needed_cols)}")

    if args.n <= 0:
        sys.exit("n must be a positive integer")

    available_core = len(df) - 2  # after dropping first/last
    if args.n > available_core:
        sys.exit(
            f"n={args.n} exceeds the available positions ({available_core}) "
            "when first and last are excluded."
        )

    out = select_positions(df, args.n)
    out.to_csv(sys.stdout, index=False)


if __name__ == "__main__":
    main()
