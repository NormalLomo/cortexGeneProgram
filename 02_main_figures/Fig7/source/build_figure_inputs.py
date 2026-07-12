#!/usr/bin/env python3
"""Build the Fig. 7 display-name table from program_names.tsv.

  - display name = name_short, formatted "P{n} {name_short}"
  - confidence == brain-weak  -> trailing " *" appended (weak-evidence marker)
  - confidence == brain-sig   -> no marker
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPOSITORY_ROOT))
from workflow.root_contract import add_canonical_root_argument, resolve_canonical_root

parser = argparse.ArgumentParser(description=__doc__)
add_canonical_root_argument(parser)
args = parser.parse_args()
canonical_root = resolve_canonical_root(args.canonical_root)
CR = canonical_root / "results/crossregion_v1"
FIG = canonical_root / "results/xspecies_humanmap_v1/spatial_xspecies/figures/Fig_spatial_univ"

nm = pd.read_csv(f"{CR}/program_names.tsv", sep="\t")

rows = []
for _, r in nm.iterrows():
    n = int(r["program"])
    P = f"P{n}"
    program = f"program_{n}"
    short = str(r["name_short"]).strip()
    # Display names use the authority name_short verbatim.
    conf = str(r["confidence"]).strip()
    star = " *" if conf == "brain-weak" else ""
    curated = f"{P} {short}{star}"
    rows.append({"program": program, "P": P, "curated_name": curated})

out = pd.DataFrame(rows)
out.to_csv(FIG / "program_names_curated.tsv", sep="\t", index=False)
print("WROTE", FIG / "program_names_curated.tsv", "n=", len(out))

# --- blacklist grep on the 60 final labels ---
BLACK = ["kidney","eye","ear","taste","skeletal","cardiac","muscle","heart","liver","renal",
         "osteoclast","endoderm","nodal","retinoic","epithelial","sensory organ",
         "substantia nigra"]
print("\n=== BLACKLIST GREP (must be 0) ===")
hits = []
for _, r in out.iterrows():
    low = r["curated_name"].lower()
    for b in BLACK:
        if b in low:
            hits.append((r["P"], b, r["curated_name"]))
if hits:
    for h in hits:
        print("  HIT:", h)
else:
    print("  0 hits — clean")
print("total blacklist hits:", len(hits))

# --- consistency spot-check against public program_names.tsv labels ---
print("\n=== CONSISTENCY SPOT-CHECK (forfig vs program_names.tsv name_short) ===")
nm_short = dict(zip([f"P{int(p)}" for p in nm["program"]], nm["name_short"]))
for P in ["P40","P54","P15","P53","P51","P34"]:
    fig_name = out[out["P"]==P].iloc[0]["curated_name"]
    src_short = nm_short[P]
    ok = src_short in fig_name
    print(f"  {P}: program_names.tsv name_short='{src_short}'  forfig='{fig_name}'  match={ok}")

# --- CJK check ---
print("\n=== CJK CHECK (must be 0) ===")
cjk = [c for s in out["curated_name"] for c in s if ord(c) > 0x2E7F]
print("CJK chars in 60 labels:", len(cjk))
