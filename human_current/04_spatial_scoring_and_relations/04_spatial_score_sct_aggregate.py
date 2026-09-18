#!/usr/bin/env python
"""
AGGREGATE the 44 per-chip SCT-score parquets -> final scored parquet.

Reads _SCT_score_perchip/*.parquet (produced by 04_spatial_score_sct.py).
Streaming / off-RAM (never holds all ~5.66M bins x 60 dense at once):
  PASS 1: stream each per-chip parquet, accumulate per-program sum / sumsq / count
          -> global per-program mu/sd ACROSS ALL BINS (NEVER per-bin).
  PASS 2: re-stream, z = (raw - mu) / sd, write z + meta to
            spatial_bin50_program_score_SCT.parquet

Meta carried: bin, majorDomain, domain, region, bin_total_umi, x, y.

DOES NOT touch spatial_bin50_program_score.parquet (raw.tpm baseline / before-fix).
"""
import os, glob, gc, time
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

t0 = time.time()
def log(*a): print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)

OUTDIR  = (__import__("os").environ["CORTEX_PROGRAM_ROOT"] + "/results/crossregion_v1")
TMPDIR  = os.path.join(OUTDIR, "_SCT_score_perchip")
PROG    = [f"program_{i}" for i in range(1, 61)]
META    = ["bin", "majorDomain", "domain", "region", "bin_total_umi", "x", "y"]
P_Z     = os.path.join(OUTDIR, "spatial_bin50_program_score_SCT.parquet")

parquets = sorted(glob.glob(os.path.join(TMPDIR, "*.parquet")))
log(f"found {len(parquets)} per-chip parquets")
assert len(parquets) == 44, f"expected 44 per-chip parquets, found {len(parquets)}"

# ---- PASS 1: global per-program mean/std over ALL bins ----
log("PASS 1: streaming global per-program mean/std over ALL bins ...")
ssum = np.zeros(60, dtype=np.float64)
ssq  = np.zeros(60, dtype=np.float64)
ncnt = 0
for p in parquets:
    t = pq.read_table(p, columns=PROG)
    mat = np.column_stack([t.column(c).to_numpy(zero_copy_only=False) for c in PROG]).astype(np.float64)
    ssum += mat.sum(axis=0)
    ssq  += (mat * mat).sum(axis=0)
    ncnt += mat.shape[0]
    del t, mat; gc.collect()
mu  = ssum / ncnt
var = ssq / ncnt - mu * mu
var[var < 0] = 0.0
sd  = np.sqrt(var)
sd[sd == 0] = 1.0
log(f"PASS 1 done: total bins={ncnt}; mu[:3]={np.round(mu[:3],4)}; sd[:3]={np.round(sd[:3],4)}")

# ---- PASS 2: write z + meta incrementally ----
log("PASS 2: writing z parquet ...")
wz = None
nwr = 0
mu32, sd32 = mu.astype(np.float32), sd.astype(np.float32)
for p in parquets:
    t = pq.read_table(p)
    raw = np.column_stack([t.column(c).to_numpy(zero_copy_only=False) for c in PROG]).astype(np.float32)
    z   = ((raw - mu32) / sd32).astype(np.float32)
    meta_arrs = {c: t.column(c) for c in META}
    zt = pa.table({**meta_arrs, **{c: pa.array(z[:, j]) for j, c in enumerate(PROG)}})
    if wz is None:
        wz = pq.ParquetWriter(P_Z + ".tmp", zt.schema)
    wz.write_table(zt)
    nwr += raw.shape[0]
    del t, raw, z, zt; gc.collect()
wz.close()
os.replace(P_Z + ".tmp", P_Z)
assert nwr == ncnt
log(f"PASS 2 done: wrote {nwr} rows to {P_Z}")
log("AGGREGATE_DONE")
