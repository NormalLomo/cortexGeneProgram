#!/usr/bin/env python3
"""K=60 redundancy statistics for the 'light defense' of K.
Reads program x program activity_pearson + gene_cosine (60x60), computes:
  - # of truly redundant pairs (activity r>0.8 AND gene_cosine>0.8)
  - which programs are involved; connected-component clusters at that threshold
  - how many of 60 are distinct vs fall in redundant clusters
Honest disclosure only -- we keep K=60.
"""
import sys
import numpy as np
import pandas as pd

BASE = "CORTEX_PROGRAM_ROOT/results/crossregion_v1/markcorr/similarity"
AP = f"{BASE}/activity_pearson.csv"
GC = f"{BASE}/gene_cosine.csv"

ap = pd.read_csv(AP, index_col=0)
gc = pd.read_csv(GC, index_col=0)

# sanity: same labels / order
assert list(ap.index) == list(ap.columns), "activity_pearson not square-labelled"
assert list(gc.index) == list(gc.columns), "gene_cosine not square-labelled"
assert list(ap.index) == list(gc.index), "label order mismatch between matrices"

labels = list(ap.index)            # 'P{n} <name>'
n = len(labels)
print(f"# programs (n) = {n}")
assert n == 60, f"expected 60 programs, got {n}"

A = ap.values.astype(float)
G = gc.values.astype(float)

# distributions (context for honesty about inflated activity baseline)
iu = np.triu_indices(n, k=1)
ap_off = A[iu]
gc_off = G[iu]
print("\n=== off-diagonal distribution (1770 unordered pairs) ===")
print(f"activity_pearson: min {ap_off.min():.3f}  median {np.median(ap_off):.3f}  "
      f"mean {ap_off.mean():.3f}  max {ap_off.max():.3f}  "
      f">0.8: {(ap_off>0.8).sum()}  >0.9: {(ap_off>0.9).sum()}")
print(f"gene_cosine     : min {gc_off.min():.3f}  median {np.median(gc_off):.3f}  "
      f"mean {gc_off.mean():.3f}  max {gc_off.max():.3f}  "
      f">0.8: {(gc_off>0.8).sum()}  >0.9: {(gc_off>0.9).sum()}")

def report(thr_a, thr_g, tag):
    print(f"\n=== REDUNDANCY at activity_pearson > {thr_a} AND gene_cosine > {thr_g}  [{tag}] ===")
    mask = (A > thr_a) & (G > thr_g)
    np.fill_diagonal(mask, False)
    # redundant pairs (unordered)
    pu = np.triu(mask, k=1)
    pair_idx = np.argwhere(pu)
    n_pairs = len(pair_idx)
    print(f"redundant pairs = {n_pairs}")
    for i, j in pair_idx:
        print(f"   {labels[i]:42s} <-> {labels[j]:42s}  "
              f"r={A[i,j]:.3f}  cos={G[i,j]:.3f}")
    # connected components over the redundancy graph (undirected)
    # adjacency = symmetric mask
    adj = mask | mask.T
    visited = [False]*n
    comps = []
    for s in range(n):
        if visited[s]:
            continue
        # BFS
        stack = [s]; comp = []
        visited[s] = True
        while stack:
            u = stack.pop(); comp.append(u)
            for v in range(n):
                if adj[u, v] and not visited[v]:
                    visited[v] = True; stack.append(v)
        comps.append(comp)
    # clusters of size >=2 are redundant clusters; singletons are distinct
    redundant_clusters = [c for c in comps if len(c) >= 2]
    singletons = [c for c in comps if len(c) == 1]
    progs_in_redundant = sorted([u for c in redundant_clusters for u in c])
    n_in_red = len(progs_in_redundant)
    n_distinct_singletons = len(singletons)
    n_clusters = len(redundant_clusters)
    # "distinct count" = singletons + (one representative per cluster)
    n_distinct_axes = n_distinct_singletons + n_clusters
    print(f"\n# redundant clusters (size>=2): {n_clusters}")
    for ci, c in enumerate(sorted(redundant_clusters, key=len, reverse=True), 1):
        pids = [labels[u].split()[0] for u in sorted(c)]
        names = "; ".join(labels[u] for u in sorted(c))
        print(f"   cluster {ci} (size {len(c)}): {{{', '.join(pids)}}}")
        print(f"       -> {names}")
    print(f"\n# programs falling in redundant clusters: {n_in_red}/{n} "
          f"({[labels[u].split()[0] for u in progs_in_redundant]})")
    print(f"# singleton (distinct) programs           : {n_distinct_singletons}/{n}")
    print(f"# distinct biological axes (singletons + 1/cluster) = {n_distinct_axes}")
    return dict(thr_a=thr_a, thr_g=thr_g, n_pairs=n_pairs,
                n_clusters=n_clusters, n_in_red=n_in_red,
                n_singletons=n_distinct_singletons, n_distinct_axes=n_distinct_axes,
                clusters=[sorted([labels[u].split()[0] for u in c]) for c in redundant_clusters],
                progs_in_red=[labels[u].split()[0] for u in progs_in_redundant])

# known clusters claimed in task: oligo/myelin {P13,26,36,37,45}, IT-neuropil {P8,24,57}, vascular {P56,59}
main = report(0.8, 0.8, "primary threshold")
# sensitivity bands
report(0.9, 0.9, "stricter")
report(0.85, 0.85, "intermediate")

# check the task's named known clusters at primary threshold
print("\n=== check task-named candidate clusters at r>0.8 & cos>0.8 ===")
def pidx(p):
    for i, l in enumerate(labels):
        if l.split()[0] == p:
            return i
    return None
for name, members in [("oligo/myelin", ["P13","P26","P36","P37","P45"]),
                      ("IT-neuropil", ["P8","P24","P57"]),
                      ("vascular", ["P56","P59"])]:
    idxs = [pidx(p) for p in members]
    print(f"\n  {name} {members}:")
    for a in range(len(idxs)):
        for b in range(a+1, len(idxs)):
            i, j = idxs[a], idxs[b]
            if i is None or j is None:
                print(f"    {members[a]}-{members[b]}: MISSING label")
                continue
            red = (A[i,j]>0.8) and (G[i,j]>0.8)
            print(f"    {members[a]}-{members[b]}: r={A[i,j]:.3f} cos={G[i,j]:.3f} "
                  f"{'REDUNDANT' if red else ''}")

print("\n=== DONE ===")
