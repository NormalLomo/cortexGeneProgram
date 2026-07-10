#!/usr/bin/env python
"""Companion to panel i: export program-z feature columns for the bins that are
in the spatial EXPRESSION+Harmony UMAP embedding, so the R panel can colour the
feature insets without needing the R 'arrow' package (not installed).

Reads the new harmony embedding CSV (spatial_umap_expr_harmony.csv), filters the
SCT program-score parquet to those bins, writes a tiny CSV with P37/P34/P7.
"""
import pyarrow.parquet as pq
import pandas as pd

PROJ = "CORTEX_PROGRAM_ROOT"
EMB  = f"{PROJ}/figures/fig2/_intermediate/spatial_umap_expr_harmony.csv"
SCT  = f"{PROJ}/results/crossregion_v1/spatial_bin50_program_score_SCT.parquet"
OUT  = f"{PROJ}/figures/fig2/_intermediate/spatial_umap_feat.csv"
FEAT = ["program_37", "program_34", "program_7"]

emb_bins = set(pd.read_csv(EMB, usecols=["bin"])["bin"].astype(str))
print("embedding bins:", len(emb_bins), flush=True)

pf = pq.ParquetFile(SCT)
parts = []
for b in pf.iter_batches(batch_size=400000, columns=["bin"] + FEAT):
    d = b.to_pandas()
    d["bin"] = d["bin"].astype(str)
    d = d[d["bin"].isin(emb_bins)]
    if len(d):
        parts.append(d)
out = pd.concat(parts)
out.to_csv(OUT, index=False)
print("wrote", OUT, "rows", len(out), "cols", list(out.columns))
