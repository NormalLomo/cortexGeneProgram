#!/usr/bin/env python
"""
Metric 2a: Spatial co-localization neighborhood conservation AUROC (CORROBORATION).
Uses markcorr final progprog_byarea_median_iqr.tsv (program x program g(r) per area).
Restricted to the 50-100 um ring band, averaged. Only areas with >=3 chips (9 areas).
Same background-rank AUROC framework as metric 2b, but unit = spatial co-localization.

Output (xregion_auroc/):
  m2a_spatial_coloc_<AREA>.tsv
  m2a_spatial_neighborhood_conservation_per_program.tsv
  m2a_spatial_neigh_pairwise_auroc_matrix.tsv
  m2a_spatial_neigh_program_region_self_auroc.tsv
"""
import itertools
import numpy as np
import pandas as pd
from scipy.stats import rankdata

BASE = "CORTEX_PROGRAM_ROOT"
OUT  = f"{BASE}/results/crossregion_v1/xregion_auroc"
SRC  = f"{BASE}/results/crossregion_v1/markcorr/final/progprog_byarea_median_iqr.tsv"
K = 10
RING_LO, RING_HI = 50, 100   # ring band for neighborhood

def main():
    d = pd.read_csv(SRC, sep="\t")
    d["Ai"] = d["A"].str.replace("program_","").astype(int)
    d["Bi"] = d["B"].str.replace("program_","").astype(int)
    chips = d.groupby("area")["n_chips_area"].first()
    areas = sorted(chips[chips>=3].index)
    print(f"[m2a] areas>=3chips: {areas}", flush=True)
    progint = list(range(1,61)); P=60
    pmap = {p:i for i,p in enumerate(progint)}

    band = d[(d["ring_um"]>=RING_LO)&(d["ring_um"]<=RING_HI)]
    # per area: program x program coloc = mean log2_median_g over the band
    C = {}
    for a in areas:
        da = band[band["area"]==a]
        M = np.full((P,P), np.nan)
        agg = da.groupby(["Ai","Bi"])["log2_median_g"].mean()
        for (ai,bi), v in agg.items():
            if ai in pmap and bi in pmap:
                M[pmap[ai], pmap[bi]] = v
        # symmetrize (g(r) cross is symmetric in expectation); fill diagonal high
        M = np.where(np.isnan(M), M.T, M)
        # remaining nan -> 0 (no signal)
        M = np.nan_to_num(M, nan=0.0)
        np.fill_diagonal(M, np.nanmax(M)+1e-6)
        C[a] = M
        pd.DataFrame(M, index=progint, columns=progint).to_csv(f"{OUT}/m2a_spatial_coloc_{a}.tsv", sep="\t")
        print(f"[m2a] {a} coloc done", flush=True)

    # top-k partners per area
    topk = {}
    for a in areas:
        Ca=C[a]; dd={}
        for pi in range(P):
            v=Ca[pi].copy(); v[pi]=-np.inf
            dd[pi]=set(np.argsort(-v)[:K].tolist())
        topk[a]=dd

    self_rows=[]; ridx={r:i for i,r in enumerate(areas)}; Rn=len(areas)
    ps=np.zeros((Rn,Rn)); pn=np.zeros((Rn,Rn))
    for rA,rB in itertools.permutations(areas,2):
        CB=C[rB]
        for pi in range(P):
            partners=topk[rA][pi]
            scores=CB[pi].copy()
            mask=np.ones(P,bool); mask[pi]=False
            qs=np.where(mask)[0]
            lab=np.array([1 if q in partners else 0 for q in qs])
            sc=scores[qs]; npos=lab.sum(); nneg=len(lab)-npos
            if npos==0 or nneg==0: auc=np.nan
            else:
                order=rankdata(sc)
                auc=(order[lab==1].sum()-npos*(npos+1)/2)/(npos*nneg)
            self_rows.append((progint[pi],rA,rB,auc))
            if np.isfinite(auc):
                ps[ridx[rA],ridx[rB]]+=auc; pn[ridx[rA],ridx[rB]]+=1

    sd=pd.DataFrame(self_rows,columns=["program","regionA","regionB","self_auroc"])
    sd.to_csv(f"{OUT}/m2a_spatial_neigh_program_region_self_auroc.tsv",sep="\t",index=False)
    pp=sd.groupby("program")["self_auroc"].agg(["mean","count"]).reset_index()
    pp.columns=["program","spatial_neigh_cons_auroc","n_pairs"]
    pp=pp.sort_values("spatial_neigh_cons_auroc",ascending=False)
    pp.to_csv(f"{OUT}/m2a_spatial_neighborhood_conservation_per_program.tsv",sep="\t",index=False)
    with np.errstate(invalid="ignore"):
        pm=ps/pn
    pd.DataFrame(pm,index=areas,columns=areas).to_csv(f"{OUT}/m2a_spatial_neigh_pairwise_auroc_matrix.tsv",sep="\t")
    print("[m2a] DONE", flush=True)
    print(pp.head(6).to_string(index=False)); print("...")
    print(pp.tail(6).to_string(index=False))
    print(f"[m2a] global mean spatial neigh conservation AUROC = {pp['spatial_neigh_cons_auroc'].mean():.4f} (9 areas only, low-coverage caveat)")

if __name__=="__main__":
    main()
