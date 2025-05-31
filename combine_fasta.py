#!/usr/bin/env python3
"""
combine_sequences.py

Merge an “exon-only” FASTA with a “full-sequence” FASTA by
re-inserting 5'/3' UTR and intronic bases.

The two FASTA files must share identical headers (or at least the same
RefSeq ID and range fragment), e.g.

    >hg38_ncbiRefSeqCurated_NM_001113491.2 range=chr17:77281469-77500623 ...

Usage
-----
python combine_sequences.py -e exon.fa -f full.fa [-o combined.fa]

Positional / optional arguments
-------------------------------
  -e / --exon   path to the FASTA containing exon-only (coding-uppercase) seqs
  -f / --full   path to the FASTA containing full genomic-style seqs
  -o / --out    output FASTA file (default: combined.fa)
"""

import argparse
import re
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# --------------------------------------------------------------------------- #
# Regex helpers
# --------------------------------------------------------------------------- #
REFSEQ_RE  = re.compile(r'([NX][MR]_\d+\.\d+)')      # NM_ / NR_ / XM_ / XR_ IDs
RANGE_RE   = re.compile(r'range=(\S+)')
CODING_RE  = re.compile(r'[A-Z]')                    # coding bases are uppercase


def parse_fasta_header(header: str):
    """Return (RefSeq accession, range) or (None, None) if not found."""
    refseq_match = REFSEQ_RE.search(header)
    range_match  = RANGE_RE.search(header)

    refseq = refseq_match.group(1) if refseq_match else None
    rng    = range_match.group(1)  if range_match  else None

    if not (refseq and rng):
        sys.stderr.write(
            f"[WARN] header missing token(s): {header}\n"
        )
    return refseq, rng


def find_coding_region(sequence: str):
    """
    Return (start_idx, end_idx) — 0-based half-open — of the
    contiguous coding stretch (uppercase) within *sequence*.
    If none, return (0, len(sequence)).
    """
    start = CODING_RE.search(sequence)
    end   = CODING_RE.search(sequence[::-1])  # search from right

    if start and end:
        return start.start(), len(sequence) - end.start()
    return 0, len(sequence)


def combine_sequences(exon_file: Path, full_seq_file: Path, output_path: Path) -> Path:
    """Core logic (unchanged except for minor logging cosmetics)."""
    intron_records   = list(SeqIO.parse(exon_file, "fasta"))
    full_seq_records = list(SeqIO.parse(full_seq_file, "fasta"))

    sys.stderr.write(f"[INFO] loaded {len(intron_records)} exon records\n")
    sys.stderr.write(f"[INFO] loaded {len(full_seq_records)} full records\n")

    # Build lookup from (RefSeq, range) -> sequence string
    full_seq_dict = {}
    for record in full_seq_records:
        refseq_id, range_info = parse_fasta_header(record.description)
        if refseq_id and range_info:
            full_seq_dict[(refseq_id, range_info)] = str(record.seq)

    sys.stderr.write(f"[INFO] dictionary built with {len(full_seq_dict)} keys\n")

    combined_records = []
    for intron_record in intron_records:
        refseq_id, range_info = parse_fasta_header(intron_record.description)

        key = (refseq_id, range_info)
        if key not in full_seq_dict:
            sys.stderr.write(f"[WARN] no matching full sequence for {key}\n")
            continue

        intron_seq = str(intron_record.seq)
        full_seq   = full_seq_dict[key]

        # 1. coding region inside exon-only record
        coding_start, coding_end = find_coding_region(intron_seq)

        # 2. coding region inside full sequence
        coding_start_full = CODING_RE.search(full_seq).start()
        coding_end_full   = len(full_seq) - CODING_RE.search(full_seq[::-1]).start()

        upstream   = intron_seq[:coding_start].lower()
        downstream = intron_seq[coding_end:].lower()

        utr_5 = full_seq[len(upstream):coding_start_full].lower()
        utr_3 = full_seq[coding_end_full:len(full_seq) - len(downstream)].lower()

        combined_seq = (
            upstream +
            utr_5 +
            intron_seq[coding_start:coding_end] +
            utr_3 +
            downstream
        )

        combined_records.append(
            SeqRecord(
                Seq(combined_seq),
                id=intron_record.id,
                description=intron_record.description
            )
        )

    sys.stderr.write(f"[INFO] writing {len(combined_records)} combined records\n")
    SeqIO.write(combined_records, output_path, "fasta")
    return output_path


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Re-insert UTR/introns into exon-only sequences."
    )
    parser.add_argument("-e", "--exon", required=True, type=Path,
                        help="FASTA with exon-only sequences (uppercase = CDS)")
    parser.add_argument("-f", "--full", required=True, type=Path,
                        help="FASTA with full genomic sequences")
    parser.add_argument("-o", "--out",  default="combined.fa", type=Path,
                        help="Output FASTA (default: combined.fa)")
    args = parser.parse_args()

    if not args.exon.exists() or not args.full.exists():
        sys.exit("Error: One or both input FASTA files do not exist.")

    out_path = combine_sequences(args.exon, args.full, args.out)
    sys.stderr.write(f"[DONE] combined FASTA written to {out_path}\n")


if __name__ == "__main__":
    main()
