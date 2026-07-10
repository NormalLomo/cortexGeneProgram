#!/usr/bin/env python
"""v18: rebuild program_names_curated_forfig from the CORRECT authoritative table.

P0#4 consistency fix (羅老師): v17's program_names_curated_forfig.tsv was copied from the STALE
results/crossregion_v1/program_labels.tsv, which mis-labelled P15/P34/P51/P53 (and others) with
endothelial/chemokine/TNF/vascular identities. The figure_release main text uses
results/crossregion_v1/program_names.tsv (clean, biologically correct, E2-validated).

v18 rebuilds the forfig table from program_names.tsv:
  - display name = name_short, formatted "P{n} {name_short}"
  - confidence == brain-weak  -> trailing " *" appended (weak-evidence marker)
  - confidence == brain-sig   -> no marker
The table keeps the SAME columns/keys as v17 (program / P / curated_name) so panel a R and
panel d py read it unchanged. ONLY the display strings change; matrix/clustering/order/colours
are untouched downstream.
"""
import pandas as pd

CR = "CORTEX_PROGRAM_ROOT/results/crossregion_v1"
FIG = "CORTEX_PROGRAM_ROOT/results/xspecies_humanmap_v1/spatial_xspecies/figures/Fig_spatial_univ"

# correct authoritative table
nm = pd.read_csv(f"{CR}/program_names.tsv", sep="\t")
# stale v17 forfig (for the remap diff)
v17 = pd.read_csv(f"{FIG}/program_names_curated_forfig.tsv", sep="\t")
v17_map = dict(zip(v17["P"], v17["curated_name"]))

rows = []
remap_rows = []
for _, r in nm.iterrows():
    n = int(r["program"])
    P = f"P{n}"
    program = f"program_{n}"
    short = str(r["name_short"]).strip()
    # single-source rule (羅老師 20260615): display name MUST be the authority name_short
    # VERBATIM. The previous display-only "morphogenesis" -> "morphog." abbreviation made
    # P56 ("Blood vessel morphogenesis") differ literally from the authority name_short, so
    # it is removed. name_short is taken exactly as written in program_names.tsv.
    conf = str(r["confidence"]).strip()
    star = " *" if conf == "brain-weak" else ""
    curated = f"{P} {short}{star}"
    rows.append({"program": program, "P": P, "curated_name": curated})
    v17_name = v17_map.get(P, "")
    changed = "YES" if v17_name != curated else "no"
    remap_rows.append({
        "P": P,
        "v17_wrong_name": v17_name,
        "v18_correct_name": curated,
        "confidence": conf,
        "changed": changed,
    })

out = pd.DataFrame(rows)
out.to_csv(f"{FIG}/program_names_curated_forfig_v18.tsv", sep="\t", index=False)
print("WROTE", f"{FIG}/program_names_curated_forfig_v18.tsv", "n=", len(out))

remap = pd.DataFrame(remap_rows)
remap.to_csv(f"{FIG}/program_name_remap_v18.tsv", sep="\t", index=False)
print("WROTE", f"{FIG}/program_name_remap_v18.tsv")

# --- verify the four flagged programs ---
print("\n=== FLAGGED P15/P34/P51/P53 (v17 -> v18) ===")
for P in ["P15", "P34", "P51", "P53"]:
    row = remap[remap["P"] == P].iloc[0]
    print(f"{P}: v17='{row['v17_wrong_name']}'  ->  v18='{row['v18_correct_name']}'  ({row['confidence']})")

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

# --- consistency spot-check vs program_names.tsv (= figure_release names) ---
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
