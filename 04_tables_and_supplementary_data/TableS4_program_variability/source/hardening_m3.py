#!/usr/bin/env python
"""TEST M3 -- cohort-robust region-variable program set.

Original F3 (02_variability_classify.py): per program, one-way ANOVA of per-cell
score across ALL 14 regions; class="variable" if FDR<0.05 & eta2>0.05. -> 14 variable
programs: [6,1,14,19,18,57,35,9,4,8,10,52,37,3].

CONFOUND: region x cohort is severely confounded -- 7 regions are us-only
(FPPFC/ITG/PoCG/S1E/SMG/SPL/VLPFC), the other 7 are dual-cohort
(ACC/AG/DLPFC/M1/S1/STG/V1); edlein never appears alone. So a region ANOVA over all
14 regions partly measures the us-vs-edlein batch axis.

Hardening (three cohort-controlled re-derivations of the region-variable axis):
 (S7) DUAL-COHORT-ONLY: restrict to the 7 regions present in BOTH cohorts; re-run the
      one-way region ANOVA + eta2 there. Here cohort is (nearly) balanced across regions
      so a surviving region effect is not a pure batch artifact.
 (ANCOVA) cohort-as-covariate: additive linear model usage ~ C(subclass) + C(batch) +
      C(region) on all cells; region PARTIAL eta2 = SS(region | subclass,batch) /
      (SS(region|.) + RSS_full). This already exists in program_cohort_eta2.tsv over 14
      regions; here we ALSO recompute it restricted to the 7 dual-cohort regions (the
      clean estimate). survive if region_partial_eta2 > 0.05 (no FDR for partial-eta2,
      so we pair it with the S7 ANOVA FDR).
 (DONOR) donor-aware: of the 7 dual-cohort regions, fit a MIXED-ish check -- per program
      one-way ANOVA on DONOR-LEVEL mean usage (10 donors collapsed to donor x region
      means) across region, AND report region eta2 on donor-region means. This guards
      against pseudoreplication (1M cells inflate any F).

SURVIVE criterion (primary) = S7 region ANOVA FDR<0.05 AND eta2(7-region)>0.05 AND
dual-cohort region_partial_eta2>0.05. Report which of the 14 survive vs drop.

Outputs -> results/crossregion_v1/hardening/M3_cohort_robust.tsv (+ console verdict).
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.stats.multitest as mt

ROOT = "CORTEX_PROGRAM_ROOT"
RES  = f"{ROOT}/results/crossregion_v1"
OUT  = f"{RES}/hardening"
OBS  = f"{ROOT}/inputs/snRNA_1M_obs.csv"
SCORES = f"{ROOT}/results/cnmf_snrna_joint_full1M_v1/snrna_joint_full1M_v1_k60_cell_scores.tsv"
os.makedirs(OUT, exist_ok=True)

DUAL = ["ACC", "AG", "DLPFC", "M1", "S1", "STG", "V1"]

# original variable set
v0 = pd.read_csv(f"{RES}/program_variability.tsv", sep="\t")
VAR14 = v0[v0["class"] == "variable"]["program"].astype(int).tolist()
print("=== TEST M3: cohort-robust region-variable set ===")
print(f"  original 14 region-variable programs: {sorted(VAR14)}\n")

# scores + obs (join by index, NOT position!)
print("  loading scores + obs (index-join) ...", flush=True)
S = pd.read_csv(SCORES, sep="\t", index_col=0)
S.columns = [int(c) for c in S.columns]
obs = pd.read_csv(OBS, index_col=0, usecols=["Unnamed: 0", "batch", "region", "subclass", "donor"])
assert set(obs.index) == set(S.index)
obs = obs.reindex(S.index)
for c in ["batch", "region", "subclass", "donor"]:
    obs[c] = obs[c].astype(str)
PROGS = list(S.columns)

# ---- (S7) one-way region ANOVA on the 7 dual-cohort regions ----
print("  (S7) 7-dual-region one-way ANOVA ...", flush=True)
m7 = obs["region"].isin(DUAL).values
S7 = S.loc[m7]; reg7 = obs.loc[m7, "region"].values
rows = []
for p in PROGS:
    y = S7[p].values
    groups = [y[reg7 == r] for r in DUAL]
    F, pv = stats.f_oneway(*groups)
    gm = y.mean()
    ssb = sum(len(g) * (g.mean() - gm) ** 2 for g in groups)
    sst = ((y - gm) ** 2).sum()
    eta2 = ssb / sst if sst > 0 else np.nan
    rows.append((p, F, pv, eta2))
s7 = pd.DataFrame(rows, columns=["program", "F7", "p7", "eta2_7region"])
s7["fdr7"] = mt.multipletests(s7["p7"].fillna(1).values, method="fdr_bh")[1]

# ---- (ANCOVA) region partial-eta2 on the 7 dual regions: usage~subclass+batch+region ----
print("  (ANCOVA) dual-region partial-eta2 (usage~subclass+batch+region) ...", flush=True)
def design(factors):
    n = len(factors[0]); cols = [np.ones(n)]
    for f in factors:
        cats = pd.Categorical(f); codes = cats.codes; K = len(cats.categories)
        for k in range(1, K):
            cols.append((codes == k).astype(float))
    return np.column_stack(cols)
def rss_all(X, Y):
    beta = np.linalg.solve(X.T @ X, X.T @ Y)
    r = Y - X @ beta
    return (r * r).sum(axis=0)
sub7 = obs.loc[m7, "subclass"]; bat7 = obs.loc[m7, "batch"]; regser7 = obs.loc[m7, "region"]
Y = S7.values.astype(np.float64)
Y = (Y - Y.mean(0)) / (Y.std(0) + 1e-12)
X_sb  = design([sub7, bat7])               # subclass + batch
X_sbr = design([sub7, bat7, regser7])      # + region
rss_sb  = rss_all(X_sb, Y)
rss_sbr = rss_all(X_sbr, Y)
ss_reg = rss_sb - rss_sbr
region_partial_eta2_7 = ss_reg / (ss_reg + rss_sbr)
pe = dict(zip(PROGS, region_partial_eta2_7))

# ---- (DONOR) donor x region means -> region eta2 on donor-level means (7 dual regions) ----
print("  (DONOR) donor-level region ANOVA (7 dual regions) ...", flush=True)
don7 = obs.loc[m7, "donor"].values
dl = pd.DataFrame({"donor": don7, "region": reg7})
donor_rows = []
for p in PROGS:
    dl_p = dl.copy(); dl_p["y"] = S7[p].values
    dm = dl_p.groupby(["donor", "region"])["y"].mean().reset_index()
    groups = [dm.loc[dm.region == r, "y"].values for r in DUAL]
    groups = [g for g in groups if len(g) >= 2]
    if len(groups) < 3:
        donor_rows.append((p, np.nan, np.nan, np.nan)); continue
    F, pv = stats.f_oneway(*groups)
    allv = np.concatenate(groups); gm = allv.mean()
    ssb = sum(len(g) * (g.mean() - gm) ** 2 for g in groups)
    sst = ((allv - gm) ** 2).sum()
    eta2 = ssb / sst if sst > 0 else np.nan
    donor_rows.append((p, F, pv, eta2))
dd = pd.DataFrame(donor_rows, columns=["program", "F_donor", "p_donor", "eta2_donor"])
dd["fdr_donor"] = mt.multipletests(dd["p_donor"].fillna(1).values, method="fdr_bh")[1]

# ---- assemble + verdict ----
tab = (s7.merge(dd, on="program")
         .assign(region_partial_eta2_7=lambda d: d["program"].map(pe))
         .assign(orig_variable=lambda d: d["program"].isin(VAR14)))
# original eta2 (14-region one-way)
tab = tab.merge(v0[["program", "eta2_region", "fdr"]].rename(
        columns={"eta2_region": "eta2_14region_orig", "fdr": "fdr_14region_orig"}), on="program")
tab["survive_S7"] = (tab["fdr7"] < 0.05) & (tab["eta2_7region"] > 0.05)
tab["survive_ANCOVA"] = tab["region_partial_eta2_7"] > 0.05
tab["survive_donor"] = (tab["fdr_donor"] < 0.05) & (tab["eta2_donor"] > 0.05)
tab["survive_primary"] = tab["survive_S7"] & tab["survive_ANCOVA"]

tab = tab.sort_values(["orig_variable", "eta2_7region"], ascending=[False, False])
tab.to_csv(f"{OUT}/M3_cohort_robust.tsv", sep="\t", index=False)

v14 = tab[tab["orig_variable"]].copy()
n_surv = int(v14["survive_primary"].sum())
n_surv_donor = int((v14["survive_primary"] & v14["survive_donor"]).sum())
print(f"\n  of 14 original region-variable programs:")
print(f"    survive PRIMARY (7-region FDR<0.05 & eta2>0.05 & dual-cohort partial-eta2>0.05): {n_surv}/14")
print(f"    survive PRIMARY + donor-level FDR<0.05 & eta2>0.05: {n_surv_donor}/14")
nm = pd.read_csv(f"{RES}/program_names.tsv", sep="\t").set_index("program")["name_short"].to_dict()
print("\n  per-program (14 originals):")
print("   prog  name                         eta2_14  eta2_7  partialEta2_7  donorEta2  fdr7      survive")
for _, r in v14.iterrows():
    p = int(r["program"])
    print("   P%-3d %-28s %.3f    %.3f   %.3f          %s      %.1e   %s"
          % (p, str(nm.get(p, ""))[:28], r["eta2_14region_orig"], r["eta2_7region"],
             r["region_partial_eta2_7"],
             ("%.3f" % r["eta2_donor"]) if pd.notna(r["eta2_donor"]) else " NA ",
             r["fdr7"], "YES" if r["survive_primary"] else "drop"))
dropped = sorted(v14.loc[~v14["survive_primary"], "program"].astype(int).tolist())
kept = sorted(v14.loc[v14["survive_primary"], "program"].astype(int).tolist())
print(f"\n  SURVIVING set: {kept}")
print(f"  DROPPED set:   {dropped}")
shrink = 1 - n_surv / 14
verdict = ("ROBUST (>=11/14 survive)" if n_surv >= 11 else
           ("SHRINKS MODERATELY (%d/14 survive, %.0f%% lost)" % (n_surv, 100*shrink) if n_surv >= 7 else
            "SHRINKS SUBSTANTIALLY (only %d/14 survive, %.0f%% lost)" % (n_surv, 100*shrink)))
print(f"\n  VERDICT F3 region-variable axis: {verdict}")
print("[write]", f"{OUT}/M3_cohort_robust.tsv")
print("\nDONE M3")
