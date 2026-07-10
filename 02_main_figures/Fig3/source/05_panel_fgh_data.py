"""Derived data for panels:
 f) program vs RCTD weight hexbin pairs (3 pairs), subsample across all bins
 g) per-chip program x layer profiles -> cross-chip correlation per program
 h) WM-distance ring profiles for WM programs (rep chip + 1 more), KDTree on OLIGO seeds
"""
import pyarrow.parquet as pq
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

base = "CORTEX_PROGRAM_ROOT/results/crossregion_v1/"
work = "CORTEX_PROGRAM_ROOT/scripts/fig2/"
PROGS=[f"program_{i}" for i in range(1,61)]
CELLTYPES=["AST","CHANDELIER","ENDO","ET","L2-L3 IT LINC00507","L3-L4 IT RORB",
 "L4-L5 IT RORB","L6 CAR3","L6 CT","L6 IT","L6B","LAMP5","MICRO","NDNF","NP",
 "OLIGO","OPC","PAX6","PVALB","SST","VIP","VLMC"]
LAYER_ORDER=["ARACHNOID","L1","L2","L3","L4","L5","L6","WM"]
rng=np.random.default_rng(7)

meta=pd.read_parquet(work+"_meta_cache.parquet").set_index("bin")

# ---------- panel f: 3 matched pairs across ALL bins, subsample 120k ----------
PAIRS=[("program_37","OLIGO"),("program_56","ENDO"),("program_5","L2-L3 IT LINC00507")]
pcols=sorted({p for p,_ in PAIRS})
ccols=sorted({c for _,c in PAIRS})
pf_p=pq.ParquetFile(base+"spatial_bin50_program_score_SCT.parquet")
pf_r=pq.ParquetFile(base+"spatial_bin50_rctd_weights.parquet")
itp=pf_p.iter_batches(batch_size=400000, columns=["bin"]+pcols)
itr=pf_r.iter_batches(batch_size=400000, columns=["bin"]+ccols+["rctd_pass_mask"])
parts=[]
for bp,br in zip(itp,itr):
    dp=bp.to_pandas().set_index("bin"); dr=br.to_pandas().set_index("bin")
    m=dr["rctd_pass_mask"].to_numpy().astype(bool)
    dd=pd.concat([dp[pcols], dr[ccols]], axis=1)[m]
    # subsample ~1.1% per chunk
    k=max(1,int(len(dd)*0.012))
    idx=rng.choice(len(dd), size=min(k,len(dd)), replace=False)
    parts.append(dd.iloc[idx])
fdat=pd.concat(parts)
print("panel f subsample n=", len(fdat))
fdat.reset_index().to_csv(work+"panelf_pairs.tsv", sep="\t", index=False)

# ---------- panel g: per-chip program x layer -> cross-chip corr per program ----
pcl=pd.read_csv(work+"prog_x_layer_per_chip.tsv", sep="\t")
# require a layer present with >=50 bins to be reliable
pcl=pcl[pcl["n"]>=50]
# build for each program: chips x layers matrix, corr across chips of layer-profile
g_recs=[]
chips=sorted(pcl.chip.unique())
for p in PROGS:
    sub=pcl[pcl.program==p]
    mat=sub.pivot(index="chip", columns="majorDomain", values="mean_z")
    mat=mat.reindex(columns=[l for l in LAYER_ORDER if l in mat.columns])
    mat=mat.dropna()
    if len(mat)<5:
        continue
    # mean profile
    meanprof=mat.mean(0)
    # each chip's corr to mean profile (reproducibility)
    for ch,row in mat.iterrows():
        if row.std()>0 and meanprof.std()>0:
            r=np.corrcoef(row.values, meanprof.values)[0,1]
            g_recs.append((p,ch,r))
gdf=pd.DataFrame(g_recs, columns=["program","chip","corr_to_mean"])
gdf.to_csv(work+"panelg_reproducibility.tsv", sep="\t", index=False)
# also a program x program cross-correlation summary not needed; per-program median
gsum=gdf.groupby("program")["corr_to_mean"].agg(["median","mean","std","count"]).reset_index()
gsum.to_csv(work+"panelg_summary.tsv", sep="\t", index=False)
print("panel g programs:", gsum.shape[0], "median repro overall:", gdf.corr_to_mean.median())

# ---------- panel h: WM-distance rings, rep chip + 1 more ----------
# load choices for rep chip
REP=open(work+"_choices.txt").read().split("REP=")[1].split("\n")[0].strip()
CHIPS_H=[REP]
# add a second large 8-layer chip
cov=meta.groupby("chip")["majorDomain"].nunique()
sizes=meta[meta.chip.isin(cov[cov==8].index)].groupby("chip").size().sort_values(ascending=False)
CHIPS_H=[sizes.index[0], sizes.index[1]]
print("panel h chips:", CHIPS_H)
WMPROGS=["program_37","program_26","program_5","program_56"]  # WM + a non-WM control(L2) + vascular

# need x,y + OLIGO weight + WMPROGS for these chips
binset=set(meta.index[meta.chip.isin(CHIPS_H)])
# gather rctd OLIGO
rc=[]
pf=pq.ParquetFile(base+"spatial_bin50_rctd_weights.parquet")
for b in pf.iter_batches(batch_size=400000, columns=["bin","OLIGO","rctd_pass_mask"]):
    d=b.to_pandas(); d=d[d["bin"].isin(binset)]
    if len(d): rc.append(d)
rc=pd.concat(rc).set_index("bin")
pc=[]
# SCT scores (repoint: was raw-tpm spatial_bin50_program_score.parquet; F2 is now
# rendered on SCT throughout -- this read was missed in the Stage-1 repoint).
pf=pq.ParquetFile(base+"spatial_bin50_program_score_SCT.parquet")
for b in pf.iter_batches(batch_size=400000, columns=["bin"]+WMPROGS):
    d=b.to_pandas(); d=d[d["bin"].isin(binset)]
    if len(d): pc.append(d)
pc=pd.concat(pc).set_index("bin")
mh=meta.loc[list(binset),["chip","x","y"]]
H=mh.join(rc[["OLIGO","rctd_pass_mask"]]).join(pc[WMPROGS])
H=H[H["rctd_pass_mask"].astype(bool)]

ring_recs=[]
for ch in CHIPS_H:
    sub=H[H.chip==ch].copy()
    pts=sub[["x","y"]].to_numpy()
    seeds=pts[sub["OLIGO"].to_numpy()>=0.5]  # WM seeds = high OLIGO
    if len(seeds)<20:
        seeds=pts[sub["OLIGO"].to_numpy()>=sub["OLIGO"].quantile(0.9)]
    tree=cKDTree(seeds)
    d,_=tree.query(pts, k=1)
    sub["dist_um"]=d  # x,y already in um (50um grid)
    # rings of 50um width up to 1000um
    bins=np.arange(0,1050,50)
    sub["ring"]=pd.cut(sub["dist_um"], bins=bins, labels=bins[:-1]+25, include_lowest=True)
    # MEDIAN per WM-distance ring (displayed per-distance aggregate; robust to
    # residual per-bin UMI-inflation outliers). Was .mean().
    gg=sub.groupby("ring", observed=True)[WMPROGS].median().reset_index()
    gg["chip"]=ch
    ring_recs.append(gg)
rings=pd.concat(ring_recs)
rings.to_csv(work+"panelh_rings.tsv", sep="\t", index=False)
print("panel h rings rows:", len(rings))
print(rings.head(12).to_string())
