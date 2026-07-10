#!/usr/bin/env python3
# NAMEFIX 20260614: rewrite sp1_long.csv ptitle column to authoritative
# program_names.tsv name_short. Only the 'ptitle' field is touched; all other
# fields/rows unchanged. Streaming row-by-row (file is ~150MB).
import csv, sys, os, shutil

F = "CORTEX_PROGRAM_ROOT/results/crossregion_v1/markcorr_v2/spatial_for_R/sp1_long.csv"
TMP = F + ".namefix_tmp"

# old ptitle -> new ptitle (authoritative name_short + (Pn) keep)
MAP = {
    "Mid-deep IT neuropil/cytoskeleton (P8)": "Neurofilament cytoskeleton (pan-neuronal) (P8)",
    "Oligodendrocyte dev. (P13)":            "Oligodendrocyte development (P13)",
    "Myelination (P45)":                     "Myelination (P45)",  # already authoritative
}

n_total = 0
n_changed = 0
seen = {}
with open(F, newline="") as fin, open(TMP, "w", newline="") as fout:
    rd = csv.reader(fin)
    wr = csv.writer(fout)
    header = next(rd)
    wr.writerow(header)
    pi = header.index("ptitle")
    for row in rd:
        n_total += 1
        old = row[pi]
        seen[old] = seen.get(old, 0) + 1
        if old in MAP and MAP[old] != old:
            row[pi] = MAP[old]
            n_changed += 1
        wr.writerow(row)

# sanity: every distinct old ptitle must be a known one (else abort, don't replace)
unknown = [k for k in seen if k not in MAP]
if unknown:
    os.remove(TMP)
    print("ABORT: unexpected ptitle values present, no change made:")
    for u in unknown:
        print("   ", repr(u))
    sys.exit(1)

shutil.move(TMP, F)
print(f"rows={n_total} changed={n_changed}")
print("distinct ptitle seen (old):")
for k, v in sorted(seen.items()):
    print(f"   {v:>8}  {k!r} -> {MAP[k]!r}")
