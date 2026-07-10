#!/usr/bin/env python
"""Fig4 data prep + heavy compute.
Produces:
  - panel_c_partition.tsv  (compositional vs cell-autonomous variance partition)
  - panel_f_umap.tsv       (UMAP of L3-L4 IT RORB cells on 60 program scores, + region + p14)
  - panel_g_bootstrap.tsv  (bootstrap within-subclass region eta2 for top driver subclasses)
All read-only on source data.
"""
import os, sys, numpy as np, pandas as pd
np.random.seed(0)

RES = "CORTEX_PROGRAM_ROOT/results/crossregion_v1"
OUT = RES  # write intermediate TSVs here
PARQ = os.path.join(RES, "cell_program_region_subclass.parquet")
PROGS = [str(i) for i in range(1, 61)]

print("[load] parquet ...", flush=True)
df = pd.read_parquet(PARQ)
print("  shape", df.shape, flush=True)

progvar = pd.read_csv(os.path.join(RES, "program_variability.tsv"), sep="\t")

# ---------- helper: eta2 of region effect on a value vector within a group ----------
def region_eta2(values, region):
    """one-way ANOVA eta2 of `region` on `values`."""
    grand = values.mean()
    ss_tot = ((values - grand) ** 2).sum()
    if ss_tot <= 0:
        return 0.0
    ss_between = 0.0
    for r, idx in region.groupby(region).groups.items():
        v = values.loc[idx]
        ss_between += len(v) * (v.mean() - grand) ** 2
    return ss_between / ss_tot

# =====================================================================
# PANEL C: variance partition for top ~15 region-variable programs
#   global region eta2  = between-region SS / total SS (ignoring subclass)
#   within-subclass region eta2 (cell-autonomous, pooled) =
#       sum over subclasses of (between-region SS within subclass) / total SS
#   compositional component = global - within  (variance explained by
#       differential subclass composition across regions), clamped >=0
# =====================================================================
print("[panel C] variance partition ...", flush=True)
top_progs = progvar.sort_values("eta2_region", ascending=False).head(15)["program"].astype(str).tolist()
rows = []
region_all = df["region"]
for p in top_progs:
    v = df[p]
    grand = v.mean()
    ss_tot = ((v - grand) ** 2).sum()
    # global between-region SS
    ss_between_region = 0.0
    for r, idx in region_all.groupby(region_all).groups.items():
        vr = v.loc[idx]
        ss_between_region += len(vr) * (vr.mean() - grand) ** 2
    eta_global = ss_between_region / ss_tot
    # within-subclass region SS pooled
    ss_within_subclass_region = 0.0
    for sc, idx_sc in df.groupby("subclass").groups.items():
        vsc = v.loc[idx_sc]
        rsc = region_all.loc[idx_sc]
        gmean_sc = vsc.mean()
        for r, idx_r in rsc.groupby(rsc).groups.items():
            vr = vsc.loc[idx_r]
            ss_within_subclass_region += len(vr) * (vr.mean() - gmean_sc) ** 2
    eta_within = ss_within_subclass_region / ss_tot      # cell-autonomous
    eta_comp = max(eta_global - eta_within, 0.0)          # compositional
    rows.append(dict(program=p, eta_global=eta_global,
                     cell_autonomous=eta_within, compositional=eta_comp))
    print(f"  p{p}: global={eta_global:.3f} within={eta_within:.3f} comp={eta_comp:.3f}", flush=True)
pd.DataFrame(rows).to_csv(os.path.join(OUT, "panel_c_partition.tsv"), sep="\t", index=False)

# =====================================================================
# PANEL F: UMAP of L3-L4 IT RORB cells on the 60 program scores
# =====================================================================
print("[panel F] UMAP for L3-L4 IT RORB ...", flush=True)
SC = "L3-L4 IT RORB"
sub = df[df["subclass"] == SC].copy()
print("  n cells", len(sub), flush=True)
X = sub[PROGS].values.astype(np.float32)
# scale per program
from sklearn.preprocessing import StandardScaler
Xs = StandardScaler().fit_transform(X)
import scanpy as sc
ad = sc.AnnData(Xs)
ad.obs["region"] = sub["region"].values
ad.obs["p14"] = sub["14"].values
sc.pp.pca(ad, n_comps=30)
sc.pp.neighbors(ad, n_neighbors=15, n_pcs=30)
sc.tl.umap(ad, random_state=0)
um = ad.obsm["X_umap"]
outf = pd.DataFrame({"UMAP1": um[:, 0], "UMAP2": um[:, 1],
                     "region": sub["region"].values, "p14": sub["14"].values})
# subsample for plotting if huge
if len(outf) > 60000:
    outf = outf.sample(60000, random_state=0)
outf.to_csv(os.path.join(OUT, "panel_f_umap.tsv"), sep="\t", index=False)
print("  wrote", len(outf), "umap pts", flush=True)

# =====================================================================
# PANEL G: bootstrap within-subclass region eta2 for top driver subclasses
# =====================================================================
print("[panel G] bootstrap eta2 ...", flush=True)
TOP_SC = ["L6 CT", "L6 IT", "ET", "NP", "L6B", "AST", "OPC", "L3-L4 IT RORB"]
N_BOOT = 20
boot_rows = []
for sc_name in TOP_SC:
    sub_g = df[df["subclass"] == sc_name]
    n = len(sub_g)
    vals = sub_g[PROGS].values
    reg = sub_g["region"].values
    # only programs that are significant-ish for this subclass: take top 10 by full-data eta2
    full_eta = {}
    reg_ser = pd.Series(reg)
    for j, p in enumerate(PROGS):
        full_eta[p] = region_eta2(pd.Series(vals[:, j]), reg_ser)
    top10 = sorted(full_eta, key=full_eta.get, reverse=True)[:10]
    idx_top = [PROGS.index(p) for p in top10]
    for b in range(N_BOOT):
        bidx = np.random.randint(0, n, n)
        vb = vals[bidx][:, idx_top]
        rb = pd.Series(reg[bidx])
        for k, p in enumerate(top10):
            e = region_eta2(pd.Series(vb[:, k]), rb)
            boot_rows.append(dict(subclass=sc_name, program=p, boot=b, eta2=e,
                                  full_eta2=full_eta[p]))
    print(f"  {sc_name}: top10 progs done", flush=True)
pd.DataFrame(boot_rows).to_csv(os.path.join(OUT, "panel_g_bootstrap.tsv"), sep="\t", index=False)

print("[DONE] all prep written", flush=True)
