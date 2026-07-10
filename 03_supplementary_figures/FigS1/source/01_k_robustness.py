#!/usr/bin/env python
"""
Extended Data Fig 3 : cNMF K-robustness.

Show the cross-region variability + cell-driver conclusions (canonical K=60) hold
across K in {40,50,65,70}.

Per K (sequential, memory-safe):
  - join cell_scores (1.036M x K) with obs[region,subclass]
  - region x program mean  -> per-program cross-region variability:
        eta2_region = SS_between/SS_total via one-way ANOVA over 14 regions (script02 logic)
        CV          = std/mean of 14 region means
        class       = "variable" if (FDR<0.05 & eta2>0.05) else "stable"
  - within-subclass region eta2 (cell-driver, script03 logic) -> per-subclass median_eta2 driver rank
Cross-K program matching to K60 by gene-spectra (factor_loadings) Pearson corr, Hungarian best-match.
Concordance: (a) Spearman of eta2 ranking K vs K60 on matched programs;
             (b) eta2_Kx vs eta2_K60 scatter; (c) top-driver class recurrence; (d) top-var program laminar-IT-like?

Outputs -> results/crossregion_v1/k_robustness/ ; figure -> figures/extended/ed_fig3_k_robustness.{pdf,png}
"""
import os, json, gc
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import linear_sum_assignment
import statsmodels.stats.multitest as mt

BASE = "CORTEX_PROGRAM_ROOT"
CNMF = f"{BASE}/results/cnmf_snrna_joint_full1M_v1"
CRX  = f"{BASE}/results/crossregion_v1"
OUTD = f"{CRX}/k_robustness"
os.makedirs(OUTD, exist_ok=True)
OBS  = f"{BASE}/inputs/snRNA_1M_obs.csv"

CS = {
    40: f"{CNMF}/snrna_joint_full1M_v1_k40_cell_scores.tsv",
    50: f"{CNMF}/snrna_joint_full1M_v1_k50_cell_scores.tsv",
    60: f"{CNMF}/snrna_joint_full1M_v1_k60_cell_scores.tsv",
    65: f"{CNMF}/snrna_joint_full1M_v1_k65_k65_cell_scores.tsv",
    70: f"{CNMF}/snrna_joint_full1M_v1_k70_cell_scores.tsv",
}
FL = {
    40: f"{CNMF}/snrna_joint_full1M_v1_k40_factor_loadings.tsv",
    50: f"{CNMF}/snrna_joint_full1M_v1_k50_factor_loadings.tsv",
    60: f"{CNMF}/snrna_joint_full1M_v1_k60_factor_loadings.tsv",
    65: f"{CNMF}/snrna_joint_full1M_v1_k65_k65_factor_loadings.tsv",
    70: f"{CNMF}/snrna_joint_full1M_v1_k70_factor_loadings.tsv",
}
ALTK = [40, 50, 65, 70]

# deep excitatory + glia driver classes (canonical conclusion: these are top drivers)
DEEP_EX = {"L6 CT", "L6 IT", "ET", "NP", "L6B", "L6 CAR3"}
GLIA    = {"AST", "OLIGO", "OPC", "MICRO"}
DRIVER_CLASSES = DEEP_EX | GLIA
# laminar IT / RORB classes (canonical: top-variable program is one of these)
LAMINAR_IT = {"L3-L4 IT RORB", "L4-L5 IT RORB", "L2-L3 IT LINC00507", "L6 IT", "L6 CAR3"}
RORB_IT    = {"L3-L4 IT RORB", "L4-L5 IT RORB"}

print("[load] obs region/subclass ...", flush=True)
obs = pd.read_csv(OBS, index_col=0)[["region", "subclass"]]
print(f"  obs: {obs.shape}", flush=True)


def per_program_variability(df, progs, regs):
    """script02 logic: ANOVA eta2 across regions + CV + class."""
    rp = df.groupby("region")[progs].mean()
    rows = []
    for p in progs:
        groups = [df.loc[df.region == r, p].dropna().values for r in regs]
        F, pv = stats.f_oneway(*groups)
        allv = np.concatenate(groups)
        gm = allv.mean()
        ssb = sum(len(g) * (g.mean() - gm) ** 2 for g in groups)
        sst = ((allv - gm) ** 2).sum()
        eta2 = ssb / sst if sst > 0 else np.nan
        rows.append((p, F, pv, eta2))
    v = pd.DataFrame(rows, columns=["program", "F", "p", "eta2_region"])
    v["fdr"] = mt.multipletests(v["p"].values, method="fdr_bh")[1]
    cv = (rp.std(axis=0, ddof=1) / rp.mean(axis=0))
    cv.index = [str(i) for i in cv.index]
    v["cv"] = cv.reindex(v["program"]).values
    v["class"] = np.where((v["fdr"] < 0.05) & (v["eta2_region"] > 0.05), "variable", "stable")
    return v, rp


def within_subclass_eta2(df, progs):
    """script03 logic: within-subclass region eta2; per-subclass median driver rank."""
    subs = sorted(df.subclass.unique())
    rows = []
    for s in subs:
        d = df[df.subclass == s]
        rc = d.region.value_counts()
        use = rc[rc >= 20].index.tolist()
        if len(use) < 3:
            continue
        for p in progs:
            groups = [d.loc[d.region == r, p].values for r in use]
            F, pv = stats.f_oneway(*groups)
            gm = d[p].mean()
            ssb = sum(len(g) * (g.mean() - gm) ** 2 for g in groups)
            sst = ((d[p] - gm) ** 2).sum()
            eta2 = (ssb / sst) if sst > 0 else np.nan
            rows.append((s, p, F, pv, eta2, len(use), len(d)))
    w = pd.DataFrame(rows, columns=["subclass", "program", "F", "p", "eta2", "n_regions", "n_cells"])
    w["fdr"] = mt.multipletests(w.p.fillna(1), method="fdr_bh")[1]
    agg = w.groupby("subclass").agg(
        median_eta2=("eta2", "median"),
        mean_eta2=("eta2", "mean"),
        n_sig_programs=("fdr", lambda x: ((x < 0.05) & (w.loc[x.index, "eta2"] > 0.05)).sum()),
        total_cells=("n_cells", "first"),
    ).sort_values("median_eta2", ascending=False)
    return w, agg


def program_subclass_assoc(df, progs):
    """argmax-subclass association per program = subclass with highest mean usage (z over subclasses)."""
    sp = df.groupby("subclass")[progs].mean()             # subclass x program
    spz = (sp - sp.mean(axis=0)) / sp.std(axis=0)         # z over subclasses per program
    top = spz.idxmax(axis=0)                               # program -> subclass (most-enriched)
    return top  # Series index=program


# ---------------------------------------------------------------------------
# main per-K loop
# ---------------------------------------------------------------------------
results = {}   # K -> dict(var=DataFrame, agg=DataFrame, assoc=Series, n_cells=int)
for K in [60] + ALTK:
    print(f"\n===== K={K} =====", flush=True)
    cs = pd.read_csv(CS[K], sep="\t", index_col=0)
    progs = list(cs.columns)
    print(f"  cell_scores {cs.shape} progs[{progs[0]}..{progs[-1]}]", flush=True)
    df = cs.join(obs, how="inner")
    del cs; gc.collect()
    ov = df["region"].notna().sum()
    print(f"  joined {df.shape} ; region notna={ov} ({100*ov/len(df):.3f}%)", flush=True)
    regs = sorted(df.region.unique())
    v, rp = per_program_variability(df, progs, regs)
    w, agg = within_subclass_eta2(df, progs)
    assoc = program_subclass_assoc(df, progs)
    results[K] = dict(var=v, agg=agg, assoc=assoc, n_cells=int(len(df)), n_prog=len(progs))
    v.to_csv(f"{OUTD}/program_variability_k{K}.tsv", sep="\t", index=False)
    agg.to_csv(f"{OUTD}/subclass_driver_rank_k{K}.tsv", sep="\t")
    assoc.rename("top_subclass").to_csv(f"{OUTD}/program_top_subclass_k{K}.tsv", sep="\t")
    nvar = int((v["class"] == "variable").sum())
    print(f"  variable={nvar}/{len(v)} ; top-eta2 prog={v.sort_values('eta2_region',ascending=False).iloc[0]['program']} "
          f"eta2={v['eta2_region'].max():.4f} ; top driver={agg.index[0]} med_eta2={agg['median_eta2'].iloc[0]:.4f}", flush=True)
    del df, w; gc.collect()

# ---------------------------------------------------------------------------
# cross-K program matching to K60 via gene-spectra Pearson corr (Hungarian)
# ---------------------------------------------------------------------------
print("\n===== cross-K matching (gene-spectra Pearson, Hungarian) =====", flush=True)
fl60 = pd.read_csv(FL[60], sep="\t", index_col=0)   # 60 x genes
match = {}   # K -> DataFrame[prog_K, match_k60, corr]
for K in ALTK:
    flk = pd.read_csv(FL[K], sep="\t", index_col=0)  # K x genes
    genes = fl60.columns.intersection(flk.columns)
    A = flk[genes].values    # K x G
    B = fl60[genes].values   # 60 x G
    # pearson corr matrix K x 60
    Az = (A - A.mean(1, keepdims=True)) / A.std(1, keepdims=True)
    Bz = (B - B.mean(1, keepdims=True)) / B.std(1, keepdims=True)
    C = (Az @ Bz.T) / Az.shape[1]   # K x 60 correlation
    ri, ci = linear_sum_assignment(-C)   # maximize corr
    rows = []
    for i, j in zip(ri, ci):
        rows.append((str(flk.index[i]), str(fl60.index[j]), float(C[i, j])))
    m = pd.DataFrame(rows, columns=["prog_K", "match_k60", "corr"]).sort_values("corr", ascending=False)
    m.to_csv(f"{OUTD}/match_k{K}_to_k60.tsv", sep="\t", index=False)
    match[K] = m
    print(f"  K={K}: matched {len(m)} progs ; corr median={m['corr'].median():.3f} "
          f"min={m['corr'].min():.3f} ; n(corr>0.5)={int((m['corr']>0.5).sum())}", flush=True)
    del flk; gc.collect()

# ---------------------------------------------------------------------------
# concordance metrics
# ---------------------------------------------------------------------------
print("\n===== concordance metrics =====", flush=True)
v60 = results[60]["var"].set_index("program")
eta60 = v60["eta2_region"]

conc_rows = []
scatter = {}   # K -> DataFrame[prog_K, eta_K, prog_60, eta_60, corr]
for K in ALTK:
    m = match[K]
    vk = results[K]["var"].set_index("program")
    sub = []
    for _, r in m.iterrows():
        pk, p60, cc = r["prog_K"], r["match_k60"], r["corr"]
        if pk in vk.index and p60 in eta60.index:
            sub.append((pk, vk.loc[pk, "eta2_region"], p60, eta60.loc[p60], cc))
    sc = pd.DataFrame(sub, columns=["prog_K", "eta_K", "prog_60", "eta_60", "corr"])
    scatter[K] = sc
    # Spearman over all matched pairs
    rho_all, p_all = stats.spearmanr(sc["eta_K"], sc["eta_60"])
    # confident-match subset (corr>0.5)
    scc = sc[sc["corr"] > 0.5]
    rho_c, p_c = stats.spearmanr(scc["eta_K"], scc["eta_60"]) if len(scc) >= 4 else (np.nan, np.nan)
    # Pearson too
    pear, _ = stats.pearsonr(sc["eta_K"], sc["eta_60"])
    conc_rows.append((K, len(sc), rho_all, p_all, len(scc), rho_c, p_c, pear, m["corr"].median()))
    print(f"  K={K}: Spearman(all,n={len(sc)})={rho_all:.3f} (p={p_all:.1e}) | "
          f"Spearman(corr>0.5,n={len(scc)})={rho_c:.3f} | Pearson={pear:.3f}", flush=True)

conc = pd.DataFrame(conc_rows, columns=["K", "n_matched", "spearman_all", "p_all",
                                        "n_conf", "spearman_conf", "p_conf", "pearson_all", "match_corr_median"])
conc.to_csv(f"{OUTD}/concordance_summary.tsv", sep="\t", index=False)

# ---- top-driver class recurrence (fraction of top-N drivers that are deep-EX/glia) ----
print("\n----- top-driver class recurrence -----", flush=True)
TOPN = 6
drv_rows = []
for K in [60] + ALTK:
    agg = results[K]["agg"]
    topN = list(agg.index[:TOPN])
    n_deep = sum(s in DEEP_EX for s in topN)
    n_glia = sum(s in GLIA for s in topN)
    n_drv = sum(s in DRIVER_CLASSES for s in topN)
    drv_rows.append((K, TOPN, n_deep, n_glia, n_drv, n_drv / TOPN, ";".join(topN)))
    print(f"  K={K}: topN={topN} -> deepEX={n_deep} glia={n_glia} driver_frac={n_drv/TOPN:.2f}", flush=True)
drv = pd.DataFrame(drv_rows, columns=["K", "topN", "n_deepEX", "n_glia", "n_driver", "driver_frac", "top_subclasses"])
drv.to_csv(f"{OUTD}/topdriver_recurrence.tsv", sep="\t", index=False)

# ---- is the top-variable program laminar-IT/RORB-like at each K? ----
print("\n----- top-variable program identity -----", flush=True)
tv_rows = []
for K in [60] + ALTK:
    v = results[K]["var"]; assoc = results[K]["assoc"]
    tvp = v.sort_values("eta2_region", ascending=False).iloc[0]["program"]
    tv_eta = v["eta2_region"].max()
    sc_assoc = assoc.get(tvp, "NA")
    is_laminar = sc_assoc in LAMINAR_IT
    is_rorb = sc_assoc in RORB_IT
    # also report top-3 var programs' subclasses
    top3 = v.sort_values("eta2_region", ascending=False).head(3)["program"].tolist()
    top3_sc = [assoc.get(p, "NA") for p in top3]
    tv_rows.append((K, tvp, round(float(tv_eta), 4), sc_assoc, is_laminar, is_rorb, ";".join(f"{p}:{s}" for p, s in zip(top3, top3_sc))))
    print(f"  K={K}: top-var prog={tvp} eta2={tv_eta:.4f} subclass={sc_assoc} laminarIT={is_laminar} RORB={is_rorb} | top3={top3_sc}", flush=True)
tv = pd.DataFrame(tv_rows, columns=["K", "top_var_prog", "top_var_eta2", "assoc_subclass", "is_laminar_IT", "is_RORB", "top3_prog_subclass"])
tv.to_csv(f"{OUTD}/top_variable_program_identity.tsv", sep="\t", index=False)

# save scatter tables + a compact summary json
for K in ALTK:
    scatter[K].to_csv(f"{OUTD}/scatter_eta_k{K}_vs_k60.tsv", sep="\t", index=False)

summary = {
    "n_cells": results[60]["n_cells"],
    "K_list": [60] + ALTK,
    "n_programs": {str(K): results[K]["n_prog"] for K in [60] + ALTK},
    "concordance": conc.to_dict(orient="records"),
    "topdriver_recurrence": drv.to_dict(orient="records"),
    "top_variable_identity": tv.to_dict(orient="records"),
}
with open(f"{OUTD}/k_robustness_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print("\n===== DONE compute. Tables in", OUTD, "=====", flush=True)
print(conc.to_string(index=False))
print(drv[["K", "n_deepEX", "n_glia", "driver_frac"]].to_string(index=False))
print(tv[["K", "top_var_prog", "top_var_eta2", "assoc_subclass", "is_laminar_IT"]].to_string(index=False))
