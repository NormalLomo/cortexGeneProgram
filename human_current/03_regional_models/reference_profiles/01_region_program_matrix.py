#!/usr/bin/env python
"""
T1: region x program aggregation matrix (foundational table for cross-region figures).
Builds region x program means, per-program z-scores across regions, region x subclass x
program long table, cell counts, and a cached per-cell joined parquet for T2/T3 reuse.

Paths are supplied through the documented environment variables. Run with Python 3.
"""
import sys
import pandas as pd
import numpy as np

CS = (__import__("os").environ["CORTEX_PROGRAM_ROOT"] + "/results/cnmf_snrna_joint_full1M_v1/snrna_joint_full1M_v1_k60_cell_scores.tsv")
OBS = (__import__("os").environ["CORTEX_PROGRAM_ROOT"] + "/inputs/snRNA_1M_obs.csv")
OUT = (__import__("os").environ["CORTEX_PROGRAM_ROOT"] + "/results/crossregion_v1")

# ---------------------------------------------------------------------------
# 1. Load
# ---------------------------------------------------------------------------
print("[1] loading cell_scores ...", flush=True)
cs = pd.read_csv(CS, sep="\t", index_col=0)
# program columns come in as strings "1".."60"; keep their original labels
progs = list(cs.columns)
print("    cell_scores shape:", cs.shape, "| program cols:", progs[:3], "...", progs[-2:])

print("[1] loading obs ...", flush=True)
obs = pd.read_csv(OBS, index_col=0)
print("    obs shape:", obs.shape)
print("    obs columns:", list(obs.columns))

# ---------------------------------------------------------------------------
# 2. Verify barcode join
# ---------------------------------------------------------------------------
ov = cs.index.isin(obs.index).sum()
pct = 100.0 * ov / len(cs)
print(f"[2] barcode overlap (raw): {ov} / {len(cs)} = {pct:.4f}%", flush=True)
print("    cs index sample :", list(cs.index[:3]))
print("    obs index sample:", list(obs.index[:3]))
if pct < 90.0:
    print("[2] BLOCKED: overlap < 90%, cannot recover ~1M cells.")
    sys.exit(2)
print("[2] join verified, no normalization needed.")

# ---------------------------------------------------------------------------
# 3. Join region + subclass (inner)
# ---------------------------------------------------------------------------
df = cs.join(obs[["region", "subclass"]], how="inner")
print("[3] joined per-cell table shape:", df.shape, flush=True)
assert df["region"].notna().all(), "NaN found in region after join"
n_region_na = df["region"].isna().sum()
n_subclass_na = df["subclass"].isna().sum()
print(f"    region NaN: {n_region_na} | subclass NaN: {n_subclass_na}")

# ---------------------------------------------------------------------------
# 4. Aggregations
# ---------------------------------------------------------------------------
# region x program mean
rp = df.groupby("region")[progs].mean()
rp.to_csv(f"{OUT}/region_program_mean.tsv", sep="\t")
print("[4] region_program_mean.tsv:", rp.shape, "| NaN:", int(rp.isna().sum().sum()), flush=True)

# per-program z-score across regions (mean/std over the region means)
z = (rp - rp.mean(axis=0)) / rp.std(axis=0)
z.to_csv(f"{OUT}/program_region_zscore.tsv", sep="\t")
print("[4] program_region_zscore.tsv:", z.shape, "| NaN:", int(z.isna().sum().sum()), flush=True)

# region x subclass x program LONG
rsp = (
    df.groupby(["region", "subclass"])[progs]
    .mean()
    .reset_index()
    .melt(id_vars=["region", "subclass"], var_name="program", value_name="mean")
)
rsp.to_csv(f"{OUT}/region_subclass_program_mean.tsv", sep="\t", index=False)
print("[4] region_subclass_program_mean.tsv (long):", rsp.shape, flush=True)

# counts: per region + per (region, subclass)
cnt_region = df.groupby("region").size().rename("n_cells").reset_index()
cnt_region.insert(1, "subclass", "__ALL__")
cnt_rs = (
    df.groupby(["region", "subclass"]).size().rename("n_cells").reset_index()
)
counts = pd.concat([cnt_region, cnt_rs], ignore_index=True)
counts.to_csv(f"{OUT}/region_cell_counts.tsv", sep="\t", index=False)
print("[4] region_cell_counts.tsv:", counts.shape, flush=True)

# ---------------------------------------------------------------------------
# 5. Cache per-cell joined table
# ---------------------------------------------------------------------------
cache_path = f"{OUT}/cell_program_region_subclass.parquet"
saved = None
try:
    df.to_parquet(cache_path)
    saved = cache_path
    print("[5] cached parquet:", cache_path, flush=True)
except Exception as e:
    print("[5] parquet failed (", repr(e), "), falling back to tsv.gz", flush=True)
    cache_path = f"{OUT}/cell_program_region_subclass.tsv.gz"
    df.to_csv(cache_path, sep="\t", compression="gzip")
    saved = cache_path
    print("[5] cached tsv.gz:", cache_path, flush=True)

# ---------------------------------------------------------------------------
# 6. Validation summary
# ---------------------------------------------------------------------------
print("\n========== VALIDATION SUMMARY ==========")
print("regions (n=%d):" % rp.shape[0])
per_region = df.groupby("region").size().sort_values(ascending=False)
for r, n in per_region.items():
    print(f"    {r}\t{n}")
print("region_program_mean shape:", rp.shape, "(expect 14 x 60)")
print("cached table:", saved, "| shape:", df.shape, "(expect ~1.036M x 62)")
print("any NaN in region_program_mean:", bool(rp.isna().any().any()))
print("========================================")
