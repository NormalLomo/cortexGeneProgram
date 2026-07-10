#!/usr/bin/env python
"""Pass2 namefig: propagate new names into Fig A overview heatmap label source.
1. program_meta.csv: rebuild plabel = "P{program} {name_short}" from the UPDATED
   program_names.tsv for ALL programs (keeps single source of truth; only the 5
   renamed programs actually change). lineage/NES/fdr columns untouched.
2. avoid_names.tsv p8_short -> short new name (drives any reader; figA panel title
   is hardcoded separately and already patched in R).
"""
import pandas as pd

D="CORTEX_PROGRAM_ROOT/results/crossregion_v1"
META=f"{D}/markcorr_v2/figures/markcorr_overview/figdata/program_meta.csv"
AVN ="CORTEX_PROGRAM_ROOT/scripts/figmarkcorr_spatial/export_R/avoid_names.tsv"

names=pd.read_csv(f"{D}/program_names.tsv",sep="\t")
pid2short={int(r.program):r.name_short for r in names.itertuples()}

# --- program_meta.csv ---
meta=pd.read_csv(META)
old=dict(zip(meta["program"].astype(int), meta["plabel"]))
meta["plabel"]=meta["program"].astype(int).map(lambda n:f"P{n} {pid2short.get(n,'')}")
# also refresh name_short col so it stays consistent
meta["name_short"]=meta["program"].astype(int).map(lambda n:pid2short.get(n, ''))
meta.to_csv(META,index=False)
changed=[(n,old[n],f"P{n} {pid2short[n]}") for n in old if old[n]!=f"P{n} {pid2short.get(n,'')}"]
print("program_meta.csv plabel changed rows:")
for n,o,nw in changed: print(f"  P{n}: '{o}' -> '{nw}'")

# --- avoid_names.tsv p8_short ---
avn=pd.read_csv(AVN,sep="\t")
avn.loc[avn["key"]=="p8_short","value"]=pid2short.get(8,'')
avn.to_csv(AVN,sep="\t",index=False)
print("avoid_names.tsv p8_short ->", pid2short.get(8))
