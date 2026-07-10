"""Finalize rep chip + exemplar programs, export per-bin tables for panels a,b,c,
and hexbin/pair data for panel f, WM-distance rings for panel h."""
import pyarrow.parquet as pq
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

base = "CORTEX_PROGRAM_ROOT/results/crossregion_v1/"
work = "CORTEX_PROGRAM_ROOT/scripts/fig2/"
PROGS = [f"program_{i}" for i in range(1,61)]
CELLTYPES = ["AST","CHANDELIER","ENDO","ET","L2-L3 IT LINC00507","L3-L4 IT RORB",
 "L4-L5 IT RORB","L6 CAR3","L6 CT","L6 IT","L6B","LAMP5","MICRO","NDNF","NP",
 "OLIGO","OPC","PAX6","PVALB","SST","VIP","VLMC"]
LAYER_ORDER=["ARACHNOID","L1","L2","L3","L4","L5","L6","WM"]

# L6 specific program: from per-layer global, find prog with max at L6 among those NOT WM-peaking
g = pd.read_csv(work+"prog_x_layer_global.tsv", sep="\t").set_index("majorDomain")[PROGS]
l6score = g.loc["L6"] - g.drop("L6").max()
print("Top L6-enriched programs:\n", l6score.sort_values(ascending=False).head(5).to_string())
l6prog = l6score.idxmax()

# Final exemplar programs for panel b (6): span layers + WM + vascular + inhibitory
EXEMPLARS = {
 "program_5":"L2/3 (program 5)",
 "program_35":"L4 (program 35)",
 l6prog:f"L6 ({l6prog.replace('program_','program ')})",
 "program_37":"OLIGO / WM (program 37)",
 "program_56":"ENDO / vascular (program 56)",
 "program_29":"Inhibitory (program 29)",
}
print("EXEMPLARS:", EXEMPLARS)

meta = pd.read_parquet(work+"_meta_cache.parquet").set_index("bin")
# rep chip: among 8-layer chips, biggest with balanced layers; pick top
cov = meta.groupby("chip")["majorDomain"].nunique()
full = cov[cov==8].index
sizes = meta[meta.chip.isin(full)].groupby("chip").size().sort_values(ascending=False)
REP = sizes.index[0]
print("REP CHIP:", REP, "n=", sizes.iloc[0])
with open(work+"_choices.txt","w") as fh:
    fh.write(f"REP={REP}\n")
    for k,v in EXEMPLARS.items(): fh.write(f"EX\t{k}\t{v}\n")
    fh.write(f"L6PROG={l6prog}\n")

# ---- per-bin for rep chip: need program scores + rctd for rep chip bins ----
rep_bins = meta.index[meta.chip==REP]
rep_meta = meta.loc[rep_bins, ["chip","x","y","region","majorDomain"]].reset_index()
rep_meta.to_csv(work+"repchip_meta.tsv", sep="\t", index=False)
print("rep meta rows:", len(rep_meta))

# pull program scores for rep chip (filter by bin membership via streaming)
repset = set(rep_bins)
ex_cols = list(EXEMPLARS.keys())
chunks=[]
pf=pq.ParquetFile(base+"spatial_bin50_program_score.parquet")
for b in pf.iter_batches(batch_size=400000, columns=["bin"]+ex_cols):
    d=b.to_pandas()
    d=d[d["bin"].isin(repset)]
    if len(d): chunks.append(d)
repprog=pd.concat(chunks)
repprog.to_csv(work+"repchip_progscores.tsv", sep="\t", index=False)
print("rep prog rows:", len(repprog))

# rctd weights for rep chip (for panel c scatterpie)
chunks=[]
pf=pq.ParquetFile(base+"spatial_bin50_rctd_weights.parquet")
for b in pf.iter_batches(batch_size=400000, columns=["bin"]+CELLTYPES+["rctd_pass_mask"]):
    d=b.to_pandas()
    d=d[d["bin"].isin(repset)]
    if len(d): chunks.append(d)
reprctd=pd.concat(chunks)
reprctd.to_csv(work+"repchip_rctd.tsv", sep="\t", index=False)
print("rep rctd rows:", len(reprctd))
