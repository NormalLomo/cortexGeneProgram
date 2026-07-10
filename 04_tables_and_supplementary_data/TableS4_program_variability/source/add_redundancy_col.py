#!/usr/bin/env python3
# Add a `redundancy_cluster` column to Supp Table 2 (supp_program_validity_eta2.tsv),
# tagging each program with its near-duplicate cluster at gene-cosine > 0.85
# (from redundancy_summary.txt, intermediate threshold). Idempotent.
import csv
import pathlib

TSV = pathlib.Path(
    "CORTEX_PROGRAM_ROOT/results/crossregion_v1/supp_program_validity_eta2.tsv"
)

# gene-cosine > 0.85 near-duplicate clusters (verified from redundancy_summary.txt)
CLUSTERS = {
    "P13": "near-dup C1 (oligo/myelin: P13/P26/P37)",
    "P26": "near-dup C1 (oligo/myelin: P13/P26/P37)",
    "P37": "near-dup C1 (oligo/myelin: P13/P26/P37)",
    "P2":  "near-dup C2 (P2/P5)",
    "P5":  "near-dup C2 (P2/P5)",
    "P18": "near-dup C3 (P18/P52)",
    "P52": "near-dup C3 (P18/P52)",
    "P38": "near-dup C4 (OPC/ECM: P38/P43)",
    "P43": "near-dup C4 (OPC/ECM: P38/P43)",
    "P40": "near-dup C5 (microglia: P40/P49)",
    "P49": "near-dup C5 (microglia: P40/P49)",
}
DEFAULT = "distinct"

rows = []
with TSV.open(newline="", encoding="utf-8") as fh:
    reader = csv.reader(fh, delimiter="\t")
    header = next(reader)
    rows = list(reader)

if "redundancy_cluster" in header:
    print("SKIP: column already present")
    raise SystemExit(0)

prog_idx = header.index("program")
header.append("redundancy_cluster")
n_tagged = 0
for r in rows:
    prog = r[prog_idx]
    tag = CLUSTERS.get(prog, DEFAULT)
    if tag != DEFAULT:
        n_tagged += 1
    r.append(tag)

# backup then write
bak = TSV.with_suffix(TSV.suffix + ".bak_redundancy")
TSV.replace(bak)
with TSV.open("w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(header)
    w.writerows(rows)

print(f"ADDED redundancy_cluster column; {n_tagged} programs tagged as near-dup, "
      f"{len(rows) - n_tagged} distinct. Backup: {bak.name}")
