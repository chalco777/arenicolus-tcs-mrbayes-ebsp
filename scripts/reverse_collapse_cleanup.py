#!/usr/bin/env python3
"""
reverse_collapse_cleanup.py

This script reverses the effect of collapsing identical sequences in a
haplotype alignment by expanding each unique sequence back to the number
of individuals that share that haplotype.  It reads a NEXUS alignment along with a tab‐separated file listing each haplotype identifier and the number of individuals represented by
that identifier (the “count” column).  It writes a new NEXUS file in
which each unique sequence is repeated according to its count and the
taxon names are suffixed with an index (e.g. `_1`, `_2`, …) to
distinguish the replicated records.

Usage (example):
    python reverse_collapse_cleanup.py \
        -i ../results/phylogenetic_analysis/alignment/mtDNA_concat.nex \
        --counts ../data/sample_to_region_mtDNA.tsv \
        --output ../results/haplotypes/diversity_stats/mtDNA_concat_reversed.nex

The counts file should have at least two columns: a header row with
"sample" and "count" (tab separated), where each subsequent line gives
the taxon name and the integer count.  Missing or empty counts are
treated as `1`.
"""

import argparse
import csv
import sys
from pathlib import Path

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


def load_counts(tsv_file: Path) -> dict:
    """Load counts per sample from a TSV file.

    The TSV is expected to have at least two columns: `sample` and `count`.
    If the count field is blank or cannot be parsed as an integer, a
    default of 1 is used.

    Args:
        tsv_file: Path to the tab-separated counts file.

    Returns:
        A dictionary mapping sample names (str) to integer counts.
    """
    counts: dict[str, int] = {}
    with tsv_file.open() as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            name = row.get('sample') or row.get('Sample') or row.get('taxon')
            if not name:
                continue
            count_str = (row.get('count') or '').strip()
            try:
                count = int(count_str) if count_str else 1
            except ValueError:
                count = 1
            counts[name] = count
    return counts


def reverse_collapse_nexus(input_path: Path, counts: dict[str, int], output_path: Path) -> None:
    """Expand a collapsed NEXUS alignment according to counts.

    Reads the alignment at `input_path`, looks up each record's ID in
    `counts` (defaulting to 1 if not found), and writes a new NEXUS
    alignment to `output_path` in which each record appears as many
    times as its count.  Replicated records have their ID suffixed
    with `_1`, `_2`, etc., to remain unique.

    Args:
        input_path: Path to the input NEXUS file containing unique haplotypes.
        counts: Dictionary of taxon counts.
        output_path: Path where the expanded NEXUS should be written.
    """
    # Read alignment using Biopython.  AlignIO returns a
    # MultipleSeqAlignment object whose records have .id and .seq
    alignment = AlignIO.read(input_path, 'nexus')

    expanded_records = []
    for record in alignment:
        # Determine how many copies to make for this record
        n = counts.get(record.id, 1)
        for i in range(1, n + 1):
            # Create a copy of the SeqRecord (shallow copy retains .seq)
            new_record = record[:]
            # Assign a new ID and name with suffix
            new_id = f"{record.id}_{i}"
            new_record.id = new_id
            new_record.name = new_id
            expanded_records.append(new_record)

    # Create a new MultipleSeqAlignment from the expanded records
    expanded_alignment = MultipleSeqAlignment(expanded_records)
    # Write the expanded alignment in NEXUS format
    AlignIO.write(expanded_alignment, output_path, 'nexus')


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Expand collapsed NEXUS alignment using counts from TSV.")
    parser.add_argument("--input", "-i", type=Path, required=True, help="Path to input NEXUS file (collapsed alignment)")
    parser.add_argument("--counts", "-c", type=Path, required=True, help="TSV file with sample and count columns")
    parser.add_argument("--output", "-o", type=Path, required=True, help="Path to output NEXUS file (expanded alignment)")

    args = parser.parse_args(argv)

    counts = load_counts(args.counts)
    reverse_collapse_nexus(args.input, counts, args.output)


if __name__ == "__main__":
    main()