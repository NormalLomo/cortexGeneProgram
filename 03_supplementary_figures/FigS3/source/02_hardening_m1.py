#!/usr/bin/env python
"""TEST M1 -- circularity of program->disease enrichment.

Replicates fig8_analysis.py enrichment exactly, then runs TWO hardening variants
for the two exemplar pairs P40->Alzheimer's Disease (AD) and P1->Autism Spectrum
Disorder (ASD):

 (A) LEAVE-OUT / conditional enrichment: remove the program-DEFINING overlap genes
     (the top-N loading genes of the program that are ALSO in the disease set) from
     BOTH the program top-list AND the disease set, then recompute hypergeometric
     OR + p on the residual sets against the residual background.
        - "removing the defining genes" => drop the observed overlap genes g0 = topN(P) & D.
        - residual program list   = topN(P) \\ g0                (still N-|g0| genes)
          (we re-take the next loading genes so the draw size n stays = N, the
           "leave-out then refill" variant; ALSO report strict "drop, no refill")
        - residual disease set     = D \\ g0
        - residual background      = BG \\ g0
     If the enrichment was ONLY the defining genes, residual OR -> ~1, p -> ns.

 (B) SPECIFICITY-MATCHED background: instead of all 18,742 genes, the null background
     = genes matched to the program's loading distribution (so the test asks: among
     genes with SIMILAR loading magnitude in this program, is the disease set still
     over-represented in the very top?). We bin all 18,742 genes by this program's
     loading decile; for each disease gene present, sample size-matched non-disease
     genes from the same decile to form a matched background of equal disease-prior.
     Operationalized as: restrict BG to {disease genes} U {loading-decile-matched
     controls}, recompute OR/p. Repeat 200x, report median OR + empirical p.

FDR: we recompute the full 60 x n_disease BH-FDR grid for the leave-out variant so the
two exemplar q-values are comparable to the original pipeline.

Outputs -> results/crossregion_v1/hardening/M1_circularity.tsv (+ console verdict).
"""
import os, json
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, fisher_exact
from statsmodels.stats.multitest import multipletests

ROOT = "CORTEX_PROGRAM_ROOT"
LOAD = f"{ROOT}/results/cnmf_snrna_joint_full1M_v1/snrna_joint_full1M_v1_k60_factor_loadings.tsv"
BB   = f"{ROOT}/data/brainbase/disease_gene_associations.txt"
OUT  = f"{ROOT}/results/crossregion_v1/hardening"
os.makedirs(OUT, exist_ok=True)
N = 150
MIN_GENES = 10
FDR_SIG = 0.05
RNG = np.random.default_rng(0)

load = pd.read_csv(LOAD, sep="\t", index_col=0)
load.index = [int(i) for i in load.index]
GENES = list(load.columns)
BG = set(GENES); N_BG = len(GENES)
assert N_BG == 18742, N_BG

bb = pd.read_csv(BB, sep="\t"); bb.columns = [c.strip() for c in bb.columns]
gcol = "Gene symbol" if "Gene symbol" in bb.columns else "Gene"
mrna = bb[bb["Type"] == "mRNA"].copy()
mrna[gcol] = mrna[gcol].astype(str).str.strip()
disease_sets = {}
for dis, sub in mrna.groupby("Disease"):
    g = set(sub[gcol]) & BG
    if len(g) >= MIN_GENES:
        disease_sets[dis] = g
diseases = sorted(disease_sets, key=lambda d: -len(disease_sets[d]))

def topN(p, n=N, exclude=frozenset()):
    row = load.loc[p].drop(labels=[g for g in exclude if g in load.columns], errors="ignore")
    return list(row.sort_values(ascending=False).head(n).index)

def hyper(top_set, dis_set, n_bg):
    K = len(dis_set); n = len(top_set); x = len(set(top_set) & dis_set)
    p = hypergeom.sf(x - 1, n_bg, K, n) if x > 0 else 1.0
    a, b, c, d = x, n - x, K - x, n_bg - n - K + x
    orr, _ = fisher_exact([[a, b], [c, d]], alternative="greater")
    return p, orr, x

# ---------- ORIGINAL grid (to reproduce baseline q) ----------
def full_grid(topfn):
    rows = []
    for p in load.index:
        ts = topfn(p)
        for dis in diseases:
            pv, orr, ov = hyper(ts, disease_sets[dis], N_BG)
            rows.append(dict(program=p, disease=dis, pval=pv, odds_ratio=orr, overlap=ov))
    df = pd.DataFrame(rows)
    df["fdr"] = multipletests(df["pval"].values, method="fdr_bh")[1]
    return df

orig = full_grid(lambda p: topN(p))

EX = [(40, "Alzheimer's Disease"), (1, "Autism Spectrum Disorder")]

def get(df, p, d):
    r = df[(df.program == p) & (df.disease == d)].iloc[0]
    return r.odds_ratio, r.pval, r.fdr, int(r.overlap)

results = []
print("=== TEST M1: circularity / leave-out enrichment ===\n")
for p, dis in EX:
    D = disease_sets[dis]
    ts0 = set(topN(p))
    g0 = sorted(ts0 & D)                       # the program-DEFINING disease genes
    or0, p0, q0, ov0 = get(orig, p, dis)
    print(f"--- P{p} -> {dis} ---")
    print(f"  ORIGINAL: overlap={ov0}/{len(D)}  OR={or0:.2f}  p={p0:.2e}  FDR={q0:.3f}")
    print(f"  defining (top-{N} loading genes in {dis} set), n={len(g0)}: {g0[:25]}{'...' if len(g0)>25 else ''}")

    # (A1) leave-out, REFILL: drop g0, refill top-N from remaining genes; disease & bg lose g0
    ts_refill = set(topN(p, exclude=set(g0)))
    Dres = D - set(g0)
    BGres = N_BG - len(g0)
    pA1, orA1, ovA1 = hyper(ts_refill, Dres, BGres)

    # (A2) leave-out, STRICT (no refill): program list = ts0 minus g0 (size N-|g0|)
    ts_strict = ts0 - set(g0)
    pA2, orA2, ovA2 = hyper(ts_strict, Dres, BGres)

    # FDR for A1 across full grid (refill variant)
    gridA1 = full_grid(lambda pp: topN(pp, exclude=(set(topN(pp)) & disease_sets[dis])) if pp == p else topN(pp))
    # NOTE: only program p's row is leave-out; others unchanged -> conservative q for p
    qA1 = gridA1[(gridA1.program == p) & (gridA1.disease == dis)].iloc[0].fdr

    print(f"  (A) LEAVE-OUT refill : overlap={ovA1}  OR={orA1:.2f}  p={pA1:.2e}  FDR(grid)={qA1:.3f}")
    print(f"  (A) LEAVE-OUT strict : overlap={ovA2}  OR={orA2:.2f}  p={pA2:.2e}  (n_draw={len(ts_strict)})")

    # (B) specificity-matched background via loading deciles of THIS program
    lo = load.loc[p]
    deciles = pd.qcut(lo.rank(method="first"), 10, labels=False)   # 0..9 by loading rank
    dec_of = dict(zip(lo.index, deciles))
    dgenes_in = [g for g in D if g in dec_of]
    nondisease_by_dec = {k: [g for g in lo.index if dec_of[g] == k and g not in D] for k in range(10)}
    ts_for_match = ts0                                              # observed top-N (with defining genes)
    ors_m, ps_m = [], []
    NB = 200
    for _ in range(NB):
        ctrl = []
        for g in dgenes_in:
            pool = nondisease_by_dec[dec_of[g]]
            if pool:
                ctrl.append(RNG.choice(pool))
        bg_match = set(dgenes_in) | set(ctrl)
        nbg = len(bg_match)
        # within this matched background, draw = top-N genes that fall in bg_match
        draw = [g for g in ts_for_match if g in bg_match]
        pv, orr, ov = hyper(draw, set(dgenes_in), nbg)
        ors_m.append(orr); ps_m.append(pv)
    orB = float(np.median(ors_m)); pB = float(np.median(ps_m))
    pB_emp = float(np.mean([pp >= 0.05 for pp in ps_m]))            # frac of matched draws ns
    print(f"  (B) MATCHED-bg (loading-decile, {NB}x): median OR={orB:.2f}  median p={pB:.2e}  "
          f"frac(p>=0.05)={pB_emp:.2f}\n")

    surv_A = (qA1 < FDR_SIG) and (orA1 > 1.5)
    surv_B = (pB < FDR_SIG) and (orB > 1.5)
    verdict = "SURVIVES" if (surv_A and surv_B) else ("PARTIAL" if (surv_A or surv_B) else "COLLAPSES")
    print(f"  VERDICT P{p}->{dis}: {verdict}  (leave-out {'pass' if surv_A else 'fail'}; matched-bg {'pass' if surv_B else 'fail'})\n")

    results.append(dict(program=p, disease=dis, n_disease_genes=len(D),
                        orig_overlap=ov0, orig_OR=or0, orig_p=p0, orig_FDR=q0,
                        n_defining_genes=len(g0),
                        leaveout_refill_overlap=ovA1, leaveout_refill_OR=orA1, leaveout_refill_p=pA1, leaveout_refill_FDR=qA1,
                        leaveout_strict_overlap=ovA2, leaveout_strict_OR=orA2, leaveout_strict_p=pA2,
                        matched_bg_median_OR=orB, matched_bg_median_p=pB, matched_bg_frac_ns=pB_emp,
                        verdict=verdict))

res = pd.DataFrame(results)
res.to_csv(f"{OUT}/M1_circularity.tsv", sep="\t", index=False)
print("[write]", f"{OUT}/M1_circularity.tsv")
print("\nDONE M1")
