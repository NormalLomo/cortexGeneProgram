#!/usr/bin/env python
"""
Figure-prep derived data for xregion_auroc figures.
1. Partner-network edge tables for the 3 hits (P24,P31,P36):
   for each hit, top co-activation partners in EACH region (so the network can be
   shown side-by-side per region to visualize rewiring). Edge weight = corr; we keep
   the union of top-8 partners across regions, with per-region corr.
2. Hit spatial fields: pull spatial_bin_program_score for hit programs, join meta
   (x,y,region,chip), restrict to well-covered chips, NO per-bin normalization
   (use the score as-is; report raw score column). One chip per chosen region.
"""
import numpy as np, pandas as pd, glob, os
import pyarrow.parquet as pq

B = "CORTEX_PROGRAM_ROOT/results/crossregion_v1"
OUT = f"{B}/xregion_auroc"
HITS = [24, 31, 36]
TOPN = 8

# ---- 1. partner networks per region (co-activation) ----
regions = sorted([os.path.basename(f).replace("m2b_coact_corr_","").replace(".tsv","")
                  for f in glob.glob(f"{OUT}/m2b_coact_corr_*.tsv")])
C = {r: pd.read_csv(f"{OUT}/m2b_coact_corr_{r}.tsv", sep="\t", index_col=0) for r in regions}
names = pd.read_csv(f"{B}/program_names.tsv", sep="\t").set_index("program")["name_short"].to_dict()

edge_rows = []
for h in HITS:
    # union of top-N partners across regions
    partner_union = set()
    for r in regions:
        v = C[r].loc[h].drop(str(h))
        partner_union |= set(v.sort_values(ascending=False).head(TOPN).index.astype(int))
    for r in regions:
        for q in partner_union:
            edge_rows.append((h, int(q), r, float(C[r].loc[h, str(q)]),
                              names.get(h,f"P{h}"), names.get(q,f"P{q}")))
ed = pd.DataFrame(edge_rows, columns=["hit","partner","region","corr","hit_name","partner_name"])
ed.to_csv(f"{OUT}/fig_partner_network_edges.tsv", sep="\t", index=False)
print(f"[prep] partner edges: {len(ed)} rows, {ed['partner'].nunique()} unique partners over {len(regions)} regions")

# also a compact "is partner in region" flag table (top-8 membership) for turnover viz
mem_rows=[]
for h in HITS:
    for r in regions:
        tops=set(C[r].loc[h].drop(str(h)).sort_values(ascending=False).head(TOPN).index.astype(int))
        for q in tops:
            mem_rows.append((h,int(q),r,names.get(q,f"P{q}")))
pd.DataFrame(mem_rows,columns=["hit","partner","region","partner_name"]).to_csv(
    f"{OUT}/fig_partner_topk_membership.tsv",sep="\t",index=False)

# ---- 2. hit spatial fields ----
# choose one chip per target region with good coverage; restrict programs to hits.
meta = pd.read_parquet(f"{B}/spatial_bin50_meta.parquet", columns=["bin","chip","x","y","region","majorDomain"])
# but program score parquet is the bin-level (358k). check its bin ids align with meta 'bin'
sc_cols = ["bin"] + [f"program_{h}" for h in HITS]
score = pd.read_parquet(f"{B}/spatial_bin_program_score.parquet", columns=sc_cols)
m = score.merge(meta, on="bin", how="inner")
print(f"[prep] spatial score joined rows={len(m)} regions={sorted(m['region'].dropna().unique())}")
# pick, per hit's target region preference, the chip with most bins
# choose 3 well-covered regions spanning AP for the field panel: DLPFC, M1, V1 (good chips)
chosen_regions = ["DLPFC","M1","V1","ITG"]
m = m[m["region"].isin(chosen_regions)]
# one chip per region = the chip with most bins
keep_chips = m.groupby(["region","chip"]).size().reset_index(name="n").sort_values("n",ascending=False)
sel = keep_chips.groupby("region").first().reset_index()[["region","chip"]]
print("[prep] selected chips:\n", sel.to_string(index=False))
m = m.merge(sel, on=["region","chip"], how="inner")
m.to_csv(f"{OUT}/fig_hit_spatial_fields.tsv", sep="\t", index=False)
print(f"[prep] hit spatial field bins={len(m)} (no per-bin normalization; raw program score)")
print("[prep] DONE")
