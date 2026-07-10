#!/usr/bin/env python
"""T3: within-subclass region ANOVA (eta2).

For each (subclass, program) pair, run one-way ANOVA of the program score
across cortical regions present in that subclass (region group included only
if it has >=20 cells in that subclass; require >=3 such groups else skip).
Records F, p, eta2 = SS_between / SS_total, n_regions_used, n_cells.
BH-FDR across all (subclass, program) tests. Aggregates per subclass to rank
which cell types drive regional program differences.
"""
import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.stats.multitest as mt

OUT = "CORTEX_PROGRAM_ROOT/results/crossregion_v1"

df = pd.read_parquet(f"{OUT}/cell_program_region_subclass.parquet")
progs = [c for c in df.columns if c not in ("region", "subclass")]

subs = sorted(df.subclass.unique())
counts = {s: int((df.subclass == s).sum()) for s in subs}
print("=== subclasses (n=%d) + cell counts ===" % len(subs))
for s in subs:
    print(f"  {s}\t{counts[s]}")

rows = []
skipped = []
for s in subs:
    d = df[df.subclass == s]
    rc = d.region.value_counts()
    use = rc[rc >= 20].index.tolist()
    if len(use) < 3:
        skipped.append((s, len(use)))
        continue
    for p in progs:
        groups = [d.loc[d.region == r, p].values for r in use]
        F, pv = stats.f_oneway(*groups)
        gm = d[p].mean()
        ssb = sum(len(g) * (g.mean() - gm) ** 2 for g in groups)
        sst = ((d[p] - gm) ** 2).sum()
        eta2 = (ssb / sst) if sst > 0 else np.nan
        rows.append((s, p, F, pv, eta2, len(use), len(d)))

if skipped:
    print("=== skipped subclasses (<3 region groups with >=20 cells) ===")
    for s, n in skipped:
        print(f"  {s}\t(usable regions={n})")

w = pd.DataFrame(rows, columns=["subclass", "program", "F", "p", "eta2", "n_regions", "n_cells"])
w["fdr"] = mt.multipletests(w.p.fillna(1), method="fdr_bh")[1]
w = w[["subclass", "program", "F", "p", "fdr", "eta2", "n_regions", "n_cells"]]
w.to_csv(f"{OUT}/within_subclass_region_eta2.tsv", sep="\t", index=False)

agg = w.groupby("subclass").agg(
    median_eta2=("eta2", "median"),
    mean_eta2=("eta2", "mean"),
    n_sig_programs=("fdr", lambda x: ((x < 0.05) & (w.loc[x.index, "eta2"] > 0.05)).sum()),
    total_cells=("n_cells", "first"),
).sort_values("median_eta2", ascending=False)
agg.to_csv(f"{OUT}/subclass_driver_rank.tsv", sep="\t")

print("=== eta2 table rows: %d ; driver rank rows: %d ===" % (len(w), len(agg)))
print("=== top-10 driver subclasses by median eta2 ===")
print(agg.head(10).to_string())
print("=== top-10 (subclass,program) by eta2 ===")
top = w.sort_values("eta2", ascending=False).head(10)
print(top.to_string(index=False))
print("=== outputs written ===")
print(f"  {OUT}/within_subclass_region_eta2.tsv")
print(f"  {OUT}/subclass_driver_rank.tsv")
