"""Core aggregation: program x majorDomain (global + per-chip), pick exemplars & rep chip.
Writes small TSVs that R panels consume. Heavy parquet read done once here.

AGGREGATION STATISTIC = MEDIAN (user decision 2026-05-27): all DISPLAYED
per-domain (majorDomain/layer) program summaries use the MEDIAN across bins,
not the mean -- robust to residual per-bin UMI-inflation outliers in the SCT
scores. Per-bin spatial scatter panels are unaffected (handled elsewhere).

NOTE on column naming: the per-chip long output keeps the column name `mean_z`
and the global matrix keeps the same wide layout, ONLY so the existing R panel
scripts (fig2_svg_panels.R / 07_panels_R.R panels d,e,g + 05_panel_fgh_data.py)
consume the files unchanged. The VALUES in those columns are now MEDIANS.
Median cannot be computed from streaming sums, so we accumulate per-group
float32 chunks and take np.median per group at the end (peak RAM ~ full data
held once; box has ~2.9 GB free, data ~1.4 GB -> fits)."""
import pyarrow.parquet as pq
import pyarrow as pa
import numpy as np
import pandas as pd

base = "CORTEX_PROGRAM_ROOT/results/crossregion_v1/"
work = "CORTEX_PROGRAM_ROOT/scripts/fig2/"

LAYER_ORDER = ["ARACHNOID","L1","L2","L3","L4","L5","L6","WM"]  # pia -> WM
PROGS = [f"program_{i}" for i in range(1,61)]

meta = pd.read_parquet(work + "_meta_cache.parquet", columns=["bin","chip","x","y","region","majorDomain"])
meta = meta.set_index("bin")

# read program scores in chunks, accumulate per-group float32 chunks for MEDIAN
pf = pq.ParquetFile(base + "spatial_bin50_program_score_SCT.parquet")
from collections import defaultdict
# accumulator dicts keyed by (chip, layer) -> list of (n_i x 60) float32 chunks
chunks_cl = defaultdict(list)   # per (chip, majorDomain)
chunks_g  = defaultdict(list)   # per majorDomain (global)

n_done = 0
for batch in pf.iter_batches(batch_size=400000, columns=["bin"]+PROGS):
    df = batch.to_pandas()
    df = df.set_index("bin")
    sub = meta.loc[df.index, ["chip","majorDomain"]]
    mat = df[PROGS].to_numpy(dtype=np.float32)
    chips = sub["chip"].to_numpy()
    layers = sub["majorDomain"].to_numpy()
    # group row indices by (chip, layer) and by layer, slice mat, stash chunk
    key_cl = pd.Series(list(zip(chips, layers)))
    for k, idx in key_cl.groupby(key_cl).groups.items():
        chunks_cl[k].append(mat[np.asarray(idx)])
    ser_ly = pd.Series(layers)
    for ly, idx in ser_ly.groupby(ser_ly).groups.items():
        chunks_g[ly].append(mat[np.asarray(idx)])
    n_done += len(df)
    print("processed", n_done, flush=True)

# ----- global program x layer MEDIAN -----
rows=[]
for ly in LAYER_ORDER:
    if ly in chunks_g and len(chunks_g[ly]):
        arr = np.concatenate(chunks_g[ly], axis=0)
        med = np.median(arr, axis=0)
        rows.append([ly]+list(med))
        del arr
gdf = pd.DataFrame(rows, columns=["majorDomain"]+PROGS)
gdf.to_csv(work+"prog_x_layer_global.tsv", sep="\t", index=False)
print("wrote global matrix (MEDIAN)", gdf.shape)

# ----- per-chip program x layer MEDIAN (long format; column kept as mean_z) -----
recs=[]
for (ch,ly), clist in chunks_cl.items():
    arr = np.concatenate(clist, axis=0)
    n = arr.shape[0]
    if n>0:
        med = np.median(arr, axis=0)
        for i,p in enumerate(PROGS):
            recs.append((ch,ly,p,med[i],n))
    del arr
pcl = pd.DataFrame(recs, columns=["chip","majorDomain","program","mean_z","n"])
pcl.to_csv(work+"prog_x_layer_per_chip.tsv", sep="\t", index=False)
print("wrote per-chip long (MEDIAN, column 'mean_z' holds median)", pcl.shape)

# pick exemplar programs from global matrix: most layer-specific
G = gdf.set_index("majorDomain")[PROGS]
# for each program, which layer has max median, and the contrast (max-min)
exemplar=[]
for p in PROGS:
    col = G[p]
    exemplar.append((p, col.idxmax(), col.max()-col.min(), col.max()))
ex = pd.DataFrame(exemplar, columns=["program","peak_layer","contrast","peak_z"]).sort_values("contrast", ascending=False)
ex.to_csv(work+"program_layer_specificity.tsv", sep="\t", index=False)
print("\nTop layer-specific programs (by MEDIAN contrast):\n", ex.head(20).to_string())
print("\nProgram 37 / 26 / 56 peak:\n", ex[ex.program.isin(["program_37","program_26","program_56"])].to_string())
