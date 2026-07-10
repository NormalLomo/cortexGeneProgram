#!/usr/bin/env python
"""
Metric 1: Expression conservation AUROC (within-human, cross-14-region).
Unit = program x subclass composition profile (21/22-subclass axis).

For each (region, program): build the subclass-composition SHAPE vector over
subclasses (NOT per-cell/per-bin library normalization). We use z-score across
the subclass axis as the shape (composition de-meaned within program x region).

Similarity = Spearman across the subclass axis.
Background = for a fixed region pair (rA, rB), the distribution of similarity of
program p in rA against ALL 60 programs in rB.
AUROC = rank of the TRUE self-program (p in rA vs p in rB) among the 60 (per
fixed region pair, both directions). Program conservation score = mean over all
ordered region pairs (bidirectional) of that program's self-vs-others AUROC.

Inputs:
  region_subclass_program_mean.tsv  (region, subclass, program, mean)  long
Outputs (xregion_auroc/):
  m1_expr_conservation_per_program.tsv     program, expr_cons_auroc, n_pairs
  m1_expr_pairwise_auroc_matrix.tsv        14x14 region-pair mean self-AUROC (overall)
  m1_expr_program_region_self_auroc.tsv    program x region-pair self-AUROC (long, for distance curve)
"""
import sys, itertools
import numpy as np
import pandas as pd
from scipy.stats import rankdata

BASE = "CORTEX_PROGRAM_ROOT"
OUT  = f"{BASE}/results/crossregion_v1/xregion_auroc"
SRC  = f"{BASE}/results/crossregion_v1/region_subclass_program_mean.tsv"

def zshape(v):
    # composition shape over subclass axis: z-score (de-mean + unit sd).
    v = np.asarray(v, float)
    sd = v.std()
    if sd == 0 or not np.isfinite(sd):
        return v - v.mean()
    return (v - v.mean()) / sd

def spearman_vec(a, B):
    """spearman of vector a (len S) vs each column of B (S x P). returns len P."""
    ar = rankdata(a)
    ar = ar - ar.mean()
    out = np.empty(B.shape[1])
    for j in range(B.shape[1]):
        br = rankdata(B[:, j]); br = br - br.mean()
        denom = np.sqrt((ar**2).sum() * (br**2).sum())
        out[j] = (ar*br).sum()/denom if denom > 0 else 0.0
    return out

def main():
    df = pd.read_csv(SRC, sep="\t")
    regions = sorted(df["region"].unique())
    programs = sorted(df["program"].unique())
    subclasses = sorted(df["subclass"].unique())
    P = len(programs); S = len(subclasses); R = len(regions)
    print(f"[m1] regions={R} programs={P} subclasses={S}", flush=True)

    # build per-region matrix  M[region] : S x P  of subclass-composition shapes
    pidx = {p:i for i,p in enumerate(programs)}
    sidx = {s:i for i,s in enumerate(subclasses)}
    Mraw = {r: np.zeros((S, P)) for r in regions}
    for r, sub, prog, m in df[["region","subclass","program","mean"]].itertuples(index=False):
        Mraw[r][sidx[sub], pidx[prog]] = m
    # shape: z over subclass axis per (region, program) column
    M = {}
    for r in regions:
        Z = np.apply_along_axis(zshape, 0, Mraw[r])  # column-wise over subclass axis
        M[r] = Z

    # For each ordered region pair (rA,rB), each program p:
    #   sims = spearman( shape(p,rA) , shape(q,rB) for all q )  -> len P
    #   auroc(p) = P( true self rank ): fraction of OTHER programs with sim < sim_self
    #   (ties counted 0.5). This is the rank-based AUROC of self vs background.
    self_rows = []  # program, rA, rB, auroc
    pair_self_sum = np.zeros((R, R)); pair_self_n = np.zeros((R, R))
    ridx = {r:i for i,r in enumerate(regions)}

    for rA, rB in itertools.permutations(regions, 2):
        A = M[rA]; Bm = M[rB]
        for pi in range(P):
            sims = spearman_vec(A[:, pi], Bm)        # len P
            s_self = sims[pi]
            others = np.delete(sims, pi)
            # AUROC = P(self > other): #(self>other) + 0.5*#(ties)  / (P-1)
            gt = np.sum(s_self > others)
            eq = np.sum(s_self == others)
            auroc = (gt + 0.5*eq) / (P-1)
            self_rows.append((programs[pi], rA, rB, auroc))
            pair_self_sum[ridx[rA], ridx[rB]] += auroc
            pair_self_n[ridx[rA], ridx[rB]] += 1

    selfdf = pd.DataFrame(self_rows, columns=["program","regionA","regionB","self_auroc"])
    selfdf.to_csv(f"{OUT}/m1_expr_program_region_self_auroc.tsv", sep="\t", index=False)

    # per-program conservation = mean over all ordered pairs
    perprog = selfdf.groupby("program")["self_auroc"].agg(["mean","count"]).reset_index()
    perprog.columns = ["program","expr_cons_auroc","n_pairs"]
    perprog = perprog.sort_values("expr_cons_auroc", ascending=False)
    perprog.to_csv(f"{OUT}/m1_expr_conservation_per_program.tsv", sep="\t", index=False)

    # region-pair mean self-AUROC matrix (averaged over programs)
    with np.errstate(invalid="ignore"):
        pairmat = pair_self_sum / pair_self_n
    pm = pd.DataFrame(pairmat, index=regions, columns=regions)
    pm.to_csv(f"{OUT}/m1_expr_pairwise_auroc_matrix.tsv", sep="\t")

    print("[m1] DONE", flush=True)
    print(perprog.head(10).to_string(index=False))
    print("...")
    print(perprog.tail(5).to_string(index=False))
    print(f"[m1] global mean expr conservation AUROC = {perprog['expr_cons_auroc'].mean():.4f}")

if __name__ == "__main__":
    main()
