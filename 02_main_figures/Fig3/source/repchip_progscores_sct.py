#!/usr/bin/env python
"""Re-export rep-chip per-bin program scores for panel b from the SCT parquet.
The previous repchip_progscores.tsv was built on the raw-tpm parquet; panel b
is now rendered on SCT. Rep-chip bins are fixed (from repchip_meta.tsv) and the
panel-b program set is fixed (P34/P24/P15/P37/P56/P40). Per-bin values (NOT an
aggregate) -> taken directly from SCT, no median. Backs up nothing (caller backs
up the .tsv).
"""
import pyarrow.parquet as pq
import pandas as pd

base = "CORTEX_PROGRAM_ROOT/results/crossregion_v1/"
work = "CORTEX_PROGRAM_ROOT/scripts/fig2/"
B_PROGS = ["program_34", "program_24", "program_15",
           "program_37", "program_56", "program_40"]

repmeta = pd.read_csv(work + "repchip_meta.tsv", sep="\t", usecols=["bin"])
repset = set(repmeta["bin"].astype(str))
print("rep-chip bins:", len(repset), flush=True)

pf = pq.ParquetFile(base + "spatial_bin50_program_score_SCT.parquet")
chunks = []
for b in pf.iter_batches(batch_size=400000, columns=["bin"] + B_PROGS):
    d = b.to_pandas()
    d["bin"] = d["bin"].astype(str)
    d = d[d["bin"].isin(repset)]
    if len(d):
        chunks.append(d)
out = pd.concat(chunks)[["bin"] + B_PROGS]
out.to_csv(work + "repchip_progscores.tsv", sep="\t", index=False)
print("wrote repchip_progscores.tsv (SCT) rows", len(out), "cols", list(out.columns))
