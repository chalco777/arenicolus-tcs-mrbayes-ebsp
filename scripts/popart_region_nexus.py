#!/usr/bin/env python3
"""
This script reads a sample‑to‑region table and a Nexus alignment and
appends a PopART compatible traits block to the alignment.  
Inputs:
    * A tab‑separated file where each line has a sample identifier,
      the region code and a count column. Sample with a region of “NA” is
      is excluded.
    * A Nexus alignment. Taxon names
      repeated with suffixes like `_2`, `_3` to indicate
      multiple haplotypes of the same sample.
# Important:
    * Regions beginning with E retain their subdivision (e.g. “EA”,
      “EB”, “EC” are kept distinct).
    * If a sample appears more than once in the table (e.g. ESP9196
      with entries in two regions) the replicate counts are added
      sequentially.

Usage should be like:
    python popart_region_nexus.py <sample_to_region.tsv> <input.nex> <output.nex>
"""

import sys
import csv
from typing import Dict, List, Tuple


def parse_sample_table(path: str) -> Tuple[Dict[str, str], List[str]]:
    """Read the sample to region table and build a mapping from each
    taxon label in the Nexus file to its region.

    Returns a dictionary mapping Nexus taxon names (including
    replicate suffixes) to region codes, and a list of all region
    codes encountered (with duplicates).  The list of regions is used
    later to determine the set of unique trait labels.
    """
    region_map: Dict[str, str] = {}
    regions: List[str] = []
    replicate_counter: Dict[str, int] = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample = row["sample"].strip()
            region = row["region"].strip()
            try:
                count = int(row["count"])
            except ValueError:
                # Skip rows where count is not a number
                continue
            # Ignore outgroup samples
            if region == "NA":
                continue
            # Normalise the region code
            if region.startswith("E") and len(region) > 1:
                norm_region = region  # keep EA/EB/EC as is
            else:
                norm_region = region[0]  # take the first letter for A/B/C/D
            # Determine the starting index for replicates of this sample
            start_idx = replicate_counter.get(sample, 0) + 1
            replicate_counter[sample] = replicate_counter.get(sample, 0) + count
            # Assign each replicate a name and region
            for i in range(count):
                idx = start_idx + i
                # The first haplotype keeps the original sample name,
                # additional haplotypes get a suffix
                taxon_name = sample if idx == 1 else f"{sample}_{idx}"
                region_map[taxon_name] = norm_region
                regions.append(norm_region)
    return region_map, regions


def extract_taxa_from_nexus(path: str) -> List[str]:
    """Parse a Nexus file and return the list of taxon names in the order
    they appear in the MATRIX.  This function only reads enough of the
    file to find the MATRIX section and stops at the semicolon
    terminating the matrix.
    """
    taxa: List[str] = []
    reading_matrix = False
    with open(path, "r") as f:
        for line in f:
            stripped = line.strip()
            # Look for the start of the matrix
            if not reading_matrix and stripped.lower().startswith("matrix"):
                reading_matrix = True
                continue
            if reading_matrix:
                if stripped.startswith(";"):
                    # End of the matrix
                    break
                if not stripped:
                    continue
                parts = stripped.split()
                if parts:
                    taxa.append(parts[0])
    return taxa


def build_traits_block(taxa: List[str], region_map: Dict[str, str], region_order: List[str]) -> str:
    """Create a traits block for the given taxa.

    Each taxon will have a row with as many comma‑separated values as
    there are region categories.  The entry is '1' for the taxon's
    region and '0' for all others.
    """
    lines = []
    # Header lines follow the PopART example【741038014716694†L98-L119】
    lines.append("BEGIN TRAITS;")
    lines.append(f"  Dimensions NTRAITS={len(region_order)};")
    lines.append("  Format labels=yes missing=? separator=Comma;")
    lines.append(f"  TraitLabels {' '.join(region_order)};")
    lines.append("  Matrix")
    # Build a mapping from region code to index for quick lookup
    region_index = {r: i for i, r in enumerate(region_order)}
    for taxon in taxa:
        region = region_map.get(taxon)
        if region is None:
            # Skip taxa without a region (outgroups)
            continue
        # Start with all zeros
        row = ["0"] * len(region_order)
        # Set 1 at the appropriate position
        row[region_index[region]] = "1"
        # Join numbers with commas to match PopART format
        values = ",".join(row)
        lines.append(f"{taxon} {values}")
    lines.append("  ;")
    lines.append("END;")
    # Join with newlines
    return "\n".join(lines)


def main(sample_path: str, nexus_input: str, nexus_output: str) -> None:
    
    # first, let's read the sample table and build region assignments
    region_map, regions = parse_sample_table(sample_path)


    # We now determine the unique regions in a stable order.  Sort first by the
    # length of the region code and then alphabetically so that
    # single‑letter regions (A, B, C, D) come before multi‑letter
    # subdivisions (EA, EB, EC).
    unique_regions = sorted(set(regions), key=lambda r: (len(r), r))
    # Extract taxa names from the Nexus file
    taxa = extract_taxa_from_nexus(nexus_input)
    # Build the traits block
    traits_block = build_traits_block(taxa, region_map, unique_regions)
    
    
    
    # Write out the original Nexus content followed by the traits block
    with open(nexus_input, "r") as fin, open(nexus_output, "w") as fout:
        fout.write(fin.read().rstrip())
        fout.write("\n\n")
        fout.write(traits_block)
        fout.write("\n")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage: python3 popart_region_script.py <sample_to_region.tsv> <input.nex> <output.nex>",
            file=sys.stderr,
        )
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])