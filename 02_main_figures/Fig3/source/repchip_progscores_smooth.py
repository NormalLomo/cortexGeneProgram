#!/usr/bin/env python
"""DISPLAY-LAYER spatial smoothing for fig2 panel b (rep-chip program fields).

Reads the rep-chip per-bin program scores (repchip_progscores.tsv, built on SCT)
and the rep-chip meta (bin,chip,x,y,...). Per CHIP, builds a k-nearest-neighbour
graph in (x,y) and replaces each bin's program-z with the MEDIAN of its k nearest
bins (self included). Median (not mean) kills single-bin outlier spikes and is
consistent with the project's robustness rule.

Writes repchip_progscores_smooth.tsv (SAME schema as the raw tsv: bin + the 6
panel-b programs). Does NOT modify the raw per-bin parquet nor repchip_progscores.tsv.
"""
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

WORK = "CORTEX_PROGRAM_ROOT/scripts/fig2/"
K = 25  # ~5x5 bin50 neighbourhood (~250 um radius); preserves laminae, kills spikes

B_PROGS = ["program_34", "program_24", "program_15",
           "program_37", "program_56", "program_40"]

meta = pd.read_csv(WORK + "repchip_meta.tsv", sep="\t",
                   usecols=["bin", "chip", "x", "y"])
prog = pd.read_csv(WORK + "repchip_progscores.tsv", sep="\t")
prog["bin"] = prog["bin"].astype(str)
meta["bin"] = meta["bin"].astype(str)
df = meta.merge(prog, on="bin", how="inner")
print("merged bins:", len(df), "chips:", df["chip"].nunique(), flush=True)

out_parts = []
for chip, g in df.groupby("chip", sort=False):
    g = g.reset_index(drop=True)
    xy = g[["x", "y"]].to_numpy(float)
    k = min(K, len(g))
    nn = NearestNeighbors(n_neighbors=k).fit(xy)
    _, idx = nn.kneighbors(xy)            # (n, k) incl self
    sm = g[["bin"]].copy()
    for p in B_PROGS:
        v = g[p].to_numpy(float)
        sm[p] = np.median(v[idx], axis=1)  # MEDIAN over k nearest bins
    out_parts.append(sm)
    print(f"  chip {chip}: {len(g)} bins, k={k}", flush=True)

out = pd.concat(out_parts)[["bin"] + B_PROGS]
out.to_csv(WORK + "repchip_progscores_smooth.tsv", sep="\t", index=False)
print("wrote repchip_progscores_smooth.tsv rows", len(out),
      "K=", K, "method=KNN-median(x,y)")
