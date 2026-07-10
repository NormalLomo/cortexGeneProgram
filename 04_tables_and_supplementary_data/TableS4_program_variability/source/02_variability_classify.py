#!/usr/bin/env python
"""T2: cross-region program variability/stability classification + region-axis gradient.

Backbone: global cNMF K=60 programs ('1'..'60') over 14 cortical regions.
- One-way ANOVA (F, p) across 14 regions per program on per-cell scores.
- eta2 = SS_between / SS_total ; CV = std/mean of 14 region means.
- BH-FDR; class = "variable" if (FDR<0.05 & eta2>0.05) else "stable".
- Region axis = PC1 of region(14) x program(60) z-score matrix; per-program
  Pearson r / linear slope / p of 14 region means vs PC1 region coordinate.
"""
import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.stats.multitest as mt
from sklearn.decomposition import PCA

OUT = "CORTEX_PROGRAM_ROOT/results/crossregion_v1"

df = pd.read_parquet(f"{OUT}/cell_program_region_subclass.parquet")
rp = pd.read_csv(f"{OUT}/region_program_mean.tsv", sep="\t", index_col=0)
z = pd.read_csv(f"{OUT}/program_region_zscore.tsv", sep="\t", index_col=0)

# program columns as they appear in the region-mean matrix (strings '1'..'60')
progs = [str(c) for c in rp.columns]
regs = list(rp.index)
print(f"[info] n_programs={len(progs)} n_regions={len(regs)} n_cells={len(df)}")
print(f"[info] regions={regs}")

# ---- 1+3+4: ANOVA, eta2, FDR, class ----
rows = []
for p in progs:
    # per-region cell-level scores, NaN-dropped per group
    groups = [df.loc[df.region == r, p].dropna().values for r in regs]
    F, pv = stats.f_oneway(*groups)
    allvals = np.concatenate(groups)
    gm = allvals.mean()
    ssb = sum(len(g) * (g.mean() - gm) ** 2 for g in groups)
    sst = ((allvals - gm) ** 2).sum()
    eta2 = ssb / sst if sst > 0 else np.nan
    rows.append((p, F, pv, eta2))

v = pd.DataFrame(rows, columns=["program", "F", "p", "eta2_region"])
v["fdr"] = mt.multipletests(v["p"].values, method="fdr_bh")[1]
# CV across the 14 region means
cv = (rp.std(axis=0, ddof=1) / rp.mean(axis=0))
cv.index = [str(i) for i in cv.index]
v["cv"] = cv.reindex(v["program"]).values
v["class"] = np.where((v["fdr"] < 0.05) & (v["eta2_region"] > 0.05), "variable", "stable")
v.to_csv(f"{OUT}/program_variability.tsv", sep="\t", index=False)

# ---- 5: region axis via PCA on z-score matrix (14 regions x 60 programs) ----
zmat = z.loc[regs, [str(c) for c in z.columns]] if all(str(c) in [str(x) for x in z.columns] for c in progs) else z.loc[regs]
# align z columns to programs (z columns are same '1'..'60')
zmat = z.loc[regs]
zmat.columns = [str(c) for c in zmat.columns]
zmat = zmat[progs]
pca = PCA(n_components=4).fit(zmat.values)
coords = pca.transform(zmat.values)  # regions x PC
evr = pca.explained_variance_ratio_
print(f"[info] PCA explained var ratio (PC1-4): {np.round(evr[:4], 4)}")

rc = pd.DataFrame({"region": regs, "PC1": coords[:, 0], "PC2": coords[:, 1]})
rc.to_csv(f"{OUT}/region_pc_coords.tsv", sep="\t", index=False)

pc1 = coords[:, 0]
g = []
for p in progs:
    y = rp.loc[regs, p].values
    r, pp = stats.pearsonr(y, pc1)
    sl = np.polyfit(pc1, y, 1)[0]
    g.append((p, sl, r, pp))
grad = pd.DataFrame(g, columns=["program", "axis_slope", "axis_r", "axis_p"])
grad.to_csv(f"{OUT}/program_gradient.tsv", sep="\t", index=False)

# combined region-axis file (also requested as one optional file)
rc.to_csv(f"{OUT}/region_axis_gradient.tsv", sep="\t", index=False)

# ---- VALIDATION / REPORT ----
nvar = int((v["class"] == "variable").sum())
nsta = int((v["class"] == "stable").sum())
print(f"\n=== VALIDATION ===")
print(f"program_variability.tsv rows={len(v)} (expect 60)")
print(f"region_pc_coords.tsv rows={len(rc)} (expect 14)")
print(f"program_gradient.tsv rows={len(grad)} (expect 60)")
print(f"class counts: variable={nvar} stable={nsta}")

print(f"\n=== TOP-10 most-variable programs by eta2 ===")
top = v.sort_values("eta2_region", ascending=False).head(10)
for _, rw in top.iterrows():
    print(f"  prog {rw['program']:>3}  eta2={rw['eta2_region']:.4f}  F={rw['F']:.1f}  fdr={rw['fdr']:.2e}  class={rw['class']}")

order = rc.sort_values("PC1")["region"].tolist()
print(f"\n=== PC1 region order (low->high) ===\n  {order}")
print(f"PC1 var explained: {evr[0]:.3f} ; PC2: {evr[1]:.3f}")
print("\nDONE")
