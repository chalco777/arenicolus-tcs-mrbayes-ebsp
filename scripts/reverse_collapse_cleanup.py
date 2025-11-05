#!/usr/bin/env python3
# Regenerates previous (to deduplication) NEXUS by duplicating rows from the "counts" file.
# - For duplicates, keep the first label as-is and add suffixes _2, _3, ...

import argparse, csv, re
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, type=Path)
ap.add_argument("-c", "--counts", required=True, type=Path)
ap.add_argument("-o", "--output", required=True, type=Path)
args = ap.parse_args()

# Read counts TSV (expects a header with sample/taxon/name and count)
counts = {}
with args.counts.open(newline="") as fh:
    r = csv.DictReader(fh, delimiter="\t")
    hdr = [h.strip().lower() for h in (r.fieldnames or [])]
    if not hdr:
        raise SystemExit("Counts TSV needs a header")
    # pick name column
    for cand in ("sample", "taxon", "name"):
        if cand in hdr:
            name_col = cand
            break
    else:
        raise SystemExit("Counts TSV: missing 'sample' (or 'taxon'/'name') column")
    count_col = "count" if "count" in hdr else None
    for row in r:
        name = (row.get(name_col) or "").strip()
        if not name:
            continue
        c = row.get(count_col) if count_col else ""
        try:
            c = int(str(c).strip()) if str(c).strip() else 1
        except ValueError:
            c = 1
        counts[name] = max(1, c)

text = args.input.read_text()

# Find DATA block
m_begin = re.search(r'(?is)\bbegin\s+data\s*;', text)
if not m_begin:
    raise SystemExit("BEGIN DATA; not found")
# find the next END; after BEGIN DATA;
m_end = re.search(r'(?is)\bend\s*;', text[m_begin.end():])
if not m_end:
    raise SystemExit("END; for DATA not found")
data_start = m_begin.start()
data_end   = m_begin.end() + m_end.end()
data_block = text[data_start:data_end]

# Find MATRIX block inside DATA
m_mat_line = re.search(r'(?im)^\s*matrix\s*\r?\n', data_block)
if not m_mat_line:
    raise SystemExit("MATRIX line not found")
mat_start = m_mat_line.end()
m_semicolon = re.search(r'(?im)^\s*;\s*$', data_block[mat_start:])
if not m_semicolon:
    raise SystemExit("Closing ';' for MATRIX not found")
mat_end = mat_start + m_semicolon.start()
matrix_text = data_block[mat_start:mat_end]

# Build new matrix
new_lines = []
ntax = 0
pat = re.compile(r'^(\s*)(\S+)(\s+)(\S.*)$')  # leading, label, space(s), rest
for raw in matrix_text.splitlines():
    if not raw.strip() or raw.lstrip().startswith('['):
        new_lines.append(raw)
        continue
    m = pat.match(raw)
    if not m:
        new_lines.append(raw)  # keep as-is
        continue
    lead, label, mid, rest = m.groups()
    n = counts.get(label, 1)
    # first occurrence: original label
    new_lines.append(f"{lead}{label}{mid}{rest}")
    ntax += 1
    # duplicates with suffixes
    for k in range(2, n + 1):
        new_label = f"{label}_{k}"
        new_lines.append(f"{lead}{new_label}{mid}{rest}")
        ntax += 1

new_matrix_text = "\n".join(new_lines)

# Replace MATRIX content
new_data_block = data_block[:mat_start] + new_matrix_text + data_block[mat_end:]

# Update NTAX before MATRIX (only the last occurrence)
head = new_data_block[:m_mat_line.start()]
tail = new_data_block[m_mat_line.start():]
occ = list(re.finditer(r'(?i)(NTAX\s*=\s*)(\d+)', head))
if occ:
    a = occ[-1]
    s, e = a.span(2)
    head = head[:s] + str(ntax) + head[e:]
new_data_block = head + tail

# Stitch back to full text
out_text = text[:data_start] + new_data_block + text[data_end:]
args.output.write_text(out_text)
print(f"Done. NTAX={ntax}")
