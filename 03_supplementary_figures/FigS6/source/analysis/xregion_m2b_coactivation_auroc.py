#!/usr/bin/env python
"""
Metric 2b: Co-activation neighborhood conservation AUROC (within-human, cross-region).

From cell-level program scores (1.036M cells x 60 programs) with region label,
compute per-region program x program co-activation (Spearman across cells).
For each program p, its "neighborhood" = top-k (k=10) co-activation partners in a
reference region. Conservation = how well that neighborhood (ranking of partners)
is preserved in another region, scored as a background-rank AUROC:

  For ordered region pair (rA,rB), program p:
    partners ranked by |corr| or corr in rA (top-k by corr, exclude self).
    In rB, look at corr(p, q) for all q != p; AUROC = P( a true top-k partner
    ranks above a non-partner ), i.e. AUC of distinguishing the k partners from
    the (P-1-k) non-partners using rB co-activation strengths.
  neighborhood conservation of p = mean over all ordered region pairs.

Co-activation uses Spearman on cell-level scores within each region.

Input : cell_program_region_subclass.parquet  (cols '1'..'60', 'region', 'subclass')
Output (xregion_auroc/):
  m2b_coact_corr_<REGION>.tsv            per-region 60x60 spearman co-activation
  m2b_neighborhood_conservation_per_program.tsv   program, neigh_cons_auroc, n_pairs
  m2b_neigh_pairwise_auroc_matrix.tsv    14x14 region-pair mean self-AUROC
  m2b_neigh_program_region_self_auroc.tsv  program x region-pair AUROC (long)
  m2b_partner_turnover.tsv               program, mean Jaccard of top-k partner sets across regions
"""
import sys, itertools
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import rankdata

BASE = "CORTEX_PROGRAM_ROOT"
OUT  = f"{BASE}/results/crossregion_v1/xregion_auroc"
PARQ = f"{BASE}/results/crossregion_v1/cell_program_region_subclass.parquet"
K = 10

def spearman_corr_matrix(X):
    """X: n x P cell-level scores. Return P x P spearman corr via column ranks."""
    n, P = X.shape
    R = np.empty_like(X, dtype=float)
    for j in range(P):
        R[:, j] = rankdata(X[:, j])
    R -= R.mean(0)
    cov = R.T @ R
    d = np.sqrt(np.diag(cov))
    denom = np.outer(d, d)
    with np.errstate(invalid="ignore", divide="ignore"):
        C = cov / denom
    C[~np.isfinite(C)] = 0.0
    np.fill_diagonal(C, 1.0)
    return C

def main():
    progs = [str(i) for i in range(1, 61)]
    print("[m2b] reading parquet (region + 60 program cols)...", flush=True)
    tbl = pq.read_table(PARQ, columns=progs + ["region"])
    df = tbl.to_pandas()
    print(f"[m2b] cells={len(df)} regions={df['region'].nunique()}", flush=True)
    regions = sorted(df["region"].unique())
    P = 60
    progint = list(range(1, 61))

    # per-region co-activation matrix
    C = {}
    for r in regions:
        Xr = df.loc[df["region"] == r, progs].to_numpy(float)
        C[r] = spearman_corr_matrix(Xr)
        pd.DataFrame(C[r], index=progint, columns=progint).to_csv(
            f"{OUT}/m2b_coact_corr_{r}.tsv", sep="\t")
        print(f"[m2b] {r}: n={Xr.shape[0]} corr done", flush=True)

    # top-k partner sets per region per program (by corr, exclude self)
    topk = {}  # (r) -> dict prog_idx -> set of partner idx
    for r in regions:
        Cr = C[r]; d = {}
        for pi in range(P):
            v = Cr[pi].copy(); v[pi] = -np.inf
            d[pi] = set(np.argsort(-v)[:K].tolist())
        topk[r] = d

    # neighborhood conservation AUROC
    self_rows = []
    ridx = {r:i for i,r in enumerate(regions)}
    Rn = len(regions)
    pair_sum = np.zeros((Rn, Rn)); pair_n = np.zeros((Rn, Rn))
    for rA, rB in itertools.permutations(regions, 2):
        CB = C[rB]
        for pi in range(P):
            partners = topk[rA][pi]           # k partners defined in rA
            # in rB, score all q != pi by corr(pi,q); AUC partners vs non-partners
            scores = CB[pi].copy()
            mask = np.ones(P, bool); mask[pi] = False
            qs = np.where(mask)[0]
            lab = np.array([1 if q in partners else 0 for q in qs])
            sc  = scores[qs]
            npos = lab.sum(); nneg = len(lab) - npos
            if npos == 0 or nneg == 0:
                auc = np.nan
            else:
                order = rankdata(sc)  # average ranks
                auc = (order[lab==1].sum() - npos*(npos+1)/2) / (npos*nneg)
            self_rows.append((progint[pi], rA, rB, auc))
            if np.isfinite(auc):
                pair_sum[ridx[rA], ridx[rB]] += auc
                pair_n[ridx[rA], ridx[rB]] += 1

    selfdf = pd.DataFrame(self_rows, columns=["program","regionA","regionB","self_auroc"])
    selfdf.to_csv(f"{OUT}/m2b_neigh_program_region_self_auroc.tsv", sep="\t", index=False)

    perprog = selfdf.groupby("program")["self_auroc"].agg(["mean","count"]).reset_index()
    perprog.columns = ["program","neigh_cons_auroc","n_pairs"]
    perprog = perprog.sort_values("neigh_cons_auroc", ascending=False)
    perprog.to_csv(f"{OUT}/m2b_neighborhood_conservation_per_program.tsv", sep="\t", index=False)

    with np.errstate(invalid="ignore"):
        pairmat = pair_sum / pair_n
    pd.DataFrame(pairmat, index=regions, columns=regions).to_csv(
        f"{OUT}/m2b_neigh_pairwise_auroc_matrix.tsv", sep="\t")

    # partner turnover: mean Jaccard of top-k partner sets across all region pairs (lower=more turnover)
    turn_rows = []
    for pi in range(P):
        js = []
        for rA, rB in itertools.combinations(regions, 2):
            a = topk[rA][pi]; b = topk[rB][pi]
            j = len(a & b) / len(a | b) if (a | b) else np.nan
            js.append(j)
        turn_rows.append((progint[pi], float(np.nanmean(js)), 1.0 - float(np.nanmean(js))))
    pd.DataFrame(turn_rows, columns=["program","mean_partner_jaccard","partner_turnover"]).to_csv(
        f"{OUT}/m2b_partner_turnover.tsv", sep="\t", index=False)

    print("[m2b] DONE", flush=True)
    print(perprog.head(8).to_string(index=False))
    print("...")
    print(perprog.tail(6).to_string(index=False))
    print(f"[m2b] global mean neighborhood conservation AUROC = {perprog['neigh_cons_auroc'].mean():.4f}")

if __name__ == "__main__":
    main()
